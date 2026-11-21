"""
Credentials Manager for the Enzyme Design Pipeline.
Resolves credentials across environment variables, .env files, user profile YAMLs,
and dry-run mock providers.
"""

from __future__ import annotations

import os
import re
from pathlib import Path
from typing import Any, Dict, List, Optional

from hackathon_extensions.config.config_schema import BackendChoice, CredentialsConfig

# Try importing pyyaml, else provide lightweight line-based YAML/JSON parser
try:
    import yaml
    HAS_YAML = True
except ImportError:
    HAS_YAML = False


def _simple_yaml_load(text: str) -> Dict[str, Any]:
    """Fallback parser for basic flat or nested YAML key-value structures without external dependencies."""
    result: Dict[str, Any] = {}
    current_section: Optional[str] = None

    for raw_line in text.splitlines():
        line = raw_line.strip()
        if not line or line.startswith("#"):
            continue

        # Check section header (e.g. "eth_euler:")
        if line.endswith(":") and not ":" in line[:-1]:
            current_section = line[:-1].strip()
            result[current_section] = {}
            continue

        if ":" in line:
            parts = line.split(":", 1)
            key = parts[0].strip()
            raw_val = parts[1].strip()
            # Strip inline comment if any
            if " #" in raw_val:
                raw_val = raw_val.split(" #", 1)[0].strip()
            elif "\t#" in raw_val:
                raw_val = raw_val.split("\t#", 1)[0].strip()

            # Strip quotes
            if (raw_val.startswith('"') and raw_val.endswith('"')) or (raw_val.startswith("'") and raw_val.endswith("'")):
                val: Any = raw_val[1:-1]
            elif raw_val.lower() == "true":
                val = True
            elif raw_val.lower() == "false":
                val = False
            elif raw_val.lower() in ("null", "none", "~", ""):
                val = None
            elif raw_val.isdigit():
                val = int(raw_val)
            else:
                try:
                    val = float(raw_val)
                except ValueError:
                    val = raw_val

            if current_section and raw_line.startswith(("  ", "\t")):
                result[current_section][key] = val
            else:
                current_section = None
                result[key] = val

    return result


def parse_yaml_file(file_path: Path) -> Dict[str, Any]:
    """Parse YAML file using PyYAML or internal fallback."""
    if not file_path.exists():
        return {}
    content = file_path.read_text(encoding="utf-8")
    if HAS_YAML:
        data = yaml.safe_load(content)
        return data if isinstance(data, dict) else {}
    return _simple_yaml_load(content)


class CredentialsManager:
    """
    Manages loading and validating authentication credentials across backends.
    Priority order:
    1. User profile YAML (`hackathon_extensions/config/user_profiles/<name>.yaml`)
    2. Shared credentials YAML (`hackathon_extensions/config/credentials.yaml`)
    3. Project root `.env` / system environment variables
    4. Mock credentials (when `dry_run=True`)
    """

    def __init__(self, workspace_root: Optional[Path] = None):
        if workspace_root is None:
            # Anchor relative to this file's location: .../hackathon_extensions/config/ -> repo root
            self.workspace_root = Path(__file__).resolve().parent.parent.parent
        else:
            self.workspace_root = Path(workspace_root).resolve()

        self.config_dir = self.workspace_root / "hackathon_extensions" / "config"
        self.profiles_dir = self.config_dir / "user_profiles"
        self._load_dotenv_if_present()

    def _load_dotenv_if_present(self) -> None:
        """Loads key-value pairs from .env file into os.environ if not already present."""
        env_paths = [
            self.workspace_root / ".env",
            self.config_dir / ".env",
        ]
        for p in env_paths:
            if p.exists():
                try:
                    lines = p.read_text(encoding="utf-8").splitlines()
                    for line in lines:
                        line = line.strip()
                        if not line or line.startswith("#") or "=" not in line:
                            continue
                        k, v = line.split("=", 1)
                        k = k.strip()
                        v = v.strip().strip("'\"")
                        if k not in os.environ:
                            os.environ[k] = v
                except Exception:
                    pass

    def get_credentials(
        self,
        user_profile: Optional[str] = None,
        dry_run: bool = False
    ) -> CredentialsConfig:
        """
        Resolve credentials merging profile, yaml, environment, and dry_run defaults.
        """
        merged_data: Dict[str, Any] = {}

        # 1. Load from credentials.yaml if present
        cred_yaml = self.config_dir / "credentials.yaml"
        if cred_yaml.exists():
            merged_data.update(parse_yaml_file(cred_yaml))

        # 2. Load from user profile if specified
        if user_profile:
            profile_path = self._resolve_profile_path(user_profile)
            if profile_path and profile_path.exists():
                profile_data = parse_yaml_file(profile_path)
                merged_data.update(profile_data)

        # 3. Overlay environment variables
        env_map = {
            "tamarind_api_key": os.getenv("TAMARIND_API_KEY", "04b8a253-a5cc-4ea6-9bfa-5ffc7c60d52a"),
            "modal_token_id": os.getenv("MODAL_TOKEN_ID"),
            "modal_token_secret": os.getenv("MODAL_TOKEN_SECRET"),
            "eth_euler_host": os.getenv("ETH_EULER_HOST"),
            "eth_euler_user": os.getenv("ETH_EULER_USER"),
            "eth_euler_port": os.getenv("ETH_EULER_PORT"),
            "eth_euler_key_filename": os.getenv("ETH_EULER_KEY") or os.getenv("ETH_EULER_KEY_FILENAME"),
            "eth_euler_password": os.getenv("ETH_EULER_PASSWORD"),
            "eth_euler_slurm_partition": os.getenv("ETH_EULER_PARTITION") or os.getenv("ETH_EULER_SLURM_PARTITION"),
            "eth_euler_scratch_dir": os.getenv("ETH_EULER_SCRATCH") or os.getenv("ETH_EULER_SCRATCH_DIR"),
        }
        for k, v in env_map.items():
            if v is not None and v != "":
                if k == "eth_euler_port":
                    try:
                        merged_data[k] = int(v)
                    except ValueError:
                        pass
                else:
                    merged_data[k] = v

        # Construct CredentialsConfig
        return CredentialsConfig(**{k: v for k, v in merged_data.items() if v is not None})

    def _resolve_profile_path(self, user_profile: str) -> Optional[Path]:
        """Finds user profile file by name or path."""
        p = Path(user_profile)
        if p.is_file():
            return p
        if not user_profile.endswith((".yaml", ".yml", ".json")):
            for ext in (".yaml", ".yml", ".json"):
                cand = self.profiles_dir / f"{user_profile}{ext}"
                if cand.is_file():
                    return cand
        cand = self.profiles_dir / user_profile
        if cand.is_file():
            return cand
        return None

    @staticmethod
    def validate_backend_credentials(
        creds: CredentialsConfig,
        backend: BackendChoice,
        dry_run: bool = False
    ) -> bool:
        """
        Validates if credentials are functional for the chosen backend.
        Returns True if valid, False otherwise.
        """
        if dry_run:
            return True
        if backend == BackendChoice.LOCAL:
            return True
        elif backend == BackendChoice.TAMARIND:
            return bool(creds.tamarind_api_key)
        elif backend == BackendChoice.MODAL:
            return bool(creds.modal_token_id and creds.modal_token_secret)
        elif backend == BackendChoice.ETH_CLUSTER:
            has_user = bool(creds.eth_euler_user)
            has_auth = bool(creds.eth_euler_key_filename or creds.eth_euler_password)
            return has_user and has_auth
        return False
