"""
Configuration package for the Enzyme Design Pipeline.
Exports schemas, credentials management, and pipeline loader utilities.
"""

from __future__ import annotations

from pathlib import Path
from typing import Optional

from hackathon_extensions.config.config_schema import (
    BackendChoice,
    ComputeRouteConfig,
    CredentialsConfig,
    DesignMode,
    GateCutoffsConfig,
    PipelineConfig,
    ThinkingBudget,
)
from hackathon_extensions.config.credentials_manager import (
    CredentialsManager,
    parse_yaml_file,
)


def load_pipeline_config(
    config_path: Optional[str] = None,
    user_profile: Optional[str] = None,
    dry_run: Optional[bool] = None,
) -> PipelineConfig:
    """
    Convenience function to load and instantiate the pipeline configuration.
    Merges YAML configuration file with credentials and CLI overrides.
    """
    default_yaml = Path(__file__).parent / "default_pipeline.yaml"
    target_path = Path(config_path) if config_path else default_yaml

    config_data = parse_yaml_file(target_path) if target_path.exists() else {}
    if dry_run is not None:
        config_data["dry_run"] = dry_run

    cfg = PipelineConfig.from_dict(config_data)

    # Attach resolved credentials
    creds_mgr = CredentialsManager()
    cfg.credentials = creds_mgr.get_credentials(
        user_profile=user_profile,
        dry_run=cfg.dry_run
    )

    return cfg


__all__ = [
    "BackendChoice",
    "ThinkingBudget",
    "DesignMode",
    "CredentialsConfig",
    "ComputeRouteConfig",
    "GateCutoffsConfig",
    "PipelineConfig",
    "CredentialsManager",
    "parse_yaml_file",
    "load_pipeline_config",
]
