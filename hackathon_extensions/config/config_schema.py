"""
Configuration schema for the Enzyme Design Pipeline.
Provides strongly-typed Pydantic models for credentials, compute routing,
gate cutoffs, thinking budgets, and full pipeline configuration.
"""

from __future__ import annotations

import json
from enum import Enum
from typing import Any, Dict, Optional, Union

# Robust Pydantic import with standard-library fallback
try:
    from pydantic import BaseModel, Field
    PYDANTIC_AVAILABLE = True
except ImportError:
    PYDANTIC_AVAILABLE = False
    # Fallback minimal BaseModel for environments without pydantic installed
    class FieldShim:
        def __init__(self, default=None, default_factory=None, description=""):
            self.default = default
            self.default_factory = default_factory
            self.description = description

    def Field(default=..., default_factory=None, description=""):
        return FieldShim(default=default, default_factory=default_factory, description=description)

    class BaseModel:
        def __init__(self, **kwargs):
            # Assign class defaults first
            for k, v in self.__class__.__dict__.items():
                if not k.startswith("_") and not callable(v):
                    if isinstance(v, FieldShim):
                        val = v.default_factory() if v.default_factory is not None else v.default
                    else:
                        val = v
                    setattr(self, k, val)
            # Assign passed kwargs
            for k, v in kwargs.items():
                setattr(self, k, v)

        def dict(self) -> Dict[str, Any]:
            return self.model_dump()

        def model_dump(self) -> Dict[str, Any]:
            result = {}
            for k, v in self.__dict__.items():
                if k.startswith("_"):
                    continue
                if isinstance(v, BaseModel):
                    result[k] = v.model_dump()
                elif isinstance(v, Enum):
                    result[k] = v.value
                elif isinstance(v, list):
                    result[k] = [item.model_dump() if isinstance(item, BaseModel) else item for item in v]
                elif isinstance(v, dict):
                    result[k] = {
                        dk: dv.model_dump() if isinstance(dv, BaseModel) else dv
                        for dk, dv in v.items()
                    }
                else:
                    result[k] = v
            return result

        def model_dump_json(self, indent: int = 2) -> str:
            return json.dumps(self.model_dump(), indent=indent)

        def json(self, **kwargs) -> str:
            return self.model_dump_json()

        @classmethod
        def model_validate(cls, data: Dict[str, Any]):
            return cls(**data)


class BackendChoice(str, Enum):
    """Available compute backend execution environments."""
    LOCAL = "local"
    TAMARIND = "tamarind"
    MODAL = "modal"
    ETH_CLUSTER = "eth_cluster"


class ThinkingBudget(str, Enum):
    """Prescribed thinking/reasoning effort budgets for LLM guidance."""
    LOW = "low"
    MEDIUM = "medium"
    HIGH = "high"


class DesignMode(str, Enum):
    """Operational design pipeline mode."""
    FAST_SCREEN = "fast_screen"
    BALANCED = "balanced"
    DEEP_RIGOR = "deep_rigor"


class CredentialsConfig(BaseModel):
    """Credentials and connection parameters across all backends."""
    tamarind_api_key: Optional[str] = Field(default=None, description="Tamarind Bio REST API Key")
    modal_token_id: Optional[str] = Field(default=None, description="Modal AI Token ID")
    modal_token_secret: Optional[str] = Field(default=None, description="Modal AI Token Secret")
    
    # ETH Euler SSH & Slurm parameters
    eth_euler_host: str = Field(default="euler.ethz.ch", description="ETH Euler login hostname")
    eth_euler_user: Optional[str] = Field(default=None, description="ETH Euler username")
    eth_euler_port: int = Field(default=22, description="SSH port for Euler")
    eth_euler_key_filename: Optional[str] = Field(default=None, description="Path to SSH private key")
    eth_euler_password: Optional[str] = Field(default=None, description="SSH password if key not used")
    eth_euler_slurm_partition: str = Field(default="gpu.2080", description="Slurm partition (RTX 2080)")
    eth_euler_scratch_dir: str = Field(default="/cluster/scratch", description="Shared scratch directory")


class ComputeRouteConfig(BaseModel):
    """
    Explicit compute backend routing for each pipeline stage and gate.
    Defaults represent the optimal speed, cost, and VRAM layout.
    """
    stage1_rfd4: BackendChoice = Field(
        default=BackendChoice.TAMARIND,
        description="Generative diffusion backbone (RFD4-Proteina). Requires >= 24GB VRAM."
    )
    stage1_rfd3: BackendChoice = Field(
        default=BackendChoice.TAMARIND,
        description="De novo structural scaffolding (RFD3). High VRAM."
    )
    stage1_ligandmpnn: BackendChoice = Field(
        default=BackendChoice.ETH_CLUSTER,
        description="Sequence design with ligand context (LigandMPNN). 8GB VRAM compatible."
    )
    gate1_cofolding: BackendChoice = Field(
        default=BackendChoice.TAMARIND,
        description="All-atom co-folding verification (Boltz-2, AF3, Protenix-v1). High VRAM."
    )
    gate2_placer: BackendChoice = Field(
        default=BackendChoice.MODAL,
        description="PLACER ligand docking and sidechain conformational ensemble (100 samples)."
    )
    gate3_thermo: BackendChoice = Field(
        default=BackendChoice.ETH_CLUSTER,
        description="Thermodynamic stability prediction (ThermoMPNN/Rosetta). Low VRAM."
    )
    gate3_esm3: BackendChoice = Field(
        default=BackendChoice.MODAL,
        description="Zero-shot sequence fitness scoring (ESM3 / SaProt). Serverless high-memory GPU."
    )
    gate4_openmm: BackendChoice = Field(
        default=BackendChoice.ETH_CLUSTER,
        description="Molecular dynamics relaxation and NAC trajectory analysis. Low-to-medium VRAM."
    )
    gate4_fieldtools: BackendChoice = Field(
        default=BackendChoice.ETH_CLUSTER,
        description="Electrostatic field projection and preorganization analysis. Low VRAM / CPU."
    )
    gate4_orca: BackendChoice = Field(
        default=BackendChoice.ETH_CLUSTER,
        description="QM/MM activation barrier estimation with ORCA. Slurm CPU/GPU cluster."
    )
    stage3_aizymes: BackendChoice = Field(
        default=BackendChoice.ETH_CLUSTER,
        description="AIzymes directed evolution simulation & multi-mutant epistasis screening."
    )


class GateCutoffsConfig(BaseModel):
    """
    Quantitative physical and biophysical thresholds required for candidate survival.
    Strictly enforced by Stage 2 verification gates.
    """
    # Gate 1: Structure & Co-folding
    rmsd_cutoff_angstrom: float = Field(
        default=1.2,
        description="Maximum backbone RMSD between generated model and catalytic motif"
    )
    attack_distance_cutoff_angstrom: float = Field(
        default=3.3,
        description="Maximum catalytic nucleophile-to-substrate electrophile attack distance"
    )
    boltz2_affinity_prob_cutoff: float = Field(
        default=0.72,
        description="Minimum Boltz-2 ligand affinity probability score"
    )
    plddt_cutoff: float = Field(
        default=70.0,
        description="Minimum average pLDDT confidence across the catalytic core"
    )
    ipae_cutoff: float = Field(
        default=10.0,
        description="Maximum interface predicted aligned error (iPAE)"
    )

    # Gate 2: PLACER Catalytic Preorganization & Docking
    placer_ensemble_size: int = Field(
        default=100,
        description="Number of conformational sidechain/ligand ensemble samples"
    )
    preorg_fraction_cutoff: float = Field(
        default=0.30,
        description="Minimum fraction of ensemble poses meeting the ideal catalytic geometry"
    )

    # Gate 3: Thermodynamic & Zero-shot Fitness
    max_ddg_destabilization: float = Field(
        default=1.5,
        description="Maximum allowable delta delta G of folding destabilization (kcal/mol)"
    )
    min_esm3_log_likelihood: float = Field(
        default=-2.5,
        description="Minimum ESM3 pseudo-log-likelihood score for evolutionary tolerance"
    )

    # Gate 4: Electric Field & QM/MM NAC Persistence
    min_nac_persistence_fraction: float = Field(
        default=0.25,
        description="Minimum fraction of MD trajectory in near-attack conformation (NAC)"
    )
    min_electric_field_projection_mv_per_cm: float = Field(
        default=10.0,
        description="Minimum oriented catalytic electric field along TS dipole (MV/cm)"
    )
    require_dG_barrier_less_than_water: bool = Field(
        default=True,
        description="Require QM/MM activation barrier lower than the uncatalyzed water reference"
    )


class PipelineConfig(BaseModel):
    """
    Master configuration combining compute routes, gate cutoffs, budgets,
    execution mode, and backend credentials.
    """
    compute_routes: ComputeRouteConfig = Field(
        default_factory=ComputeRouteConfig,
        description="Routing table for stages and gates"
    )
    gate_cutoffs: GateCutoffsConfig = Field(
        default_factory=GateCutoffsConfig,
        description="Quantitative quality cutoffs"
    )
    thinking_budget: ThinkingBudget = Field(
        default=ThinkingBudget.MEDIUM,
        description="LLM reasoning effort allocation"
    )
    design_mode: DesignMode = Field(
        default=DesignMode.BALANCED,
        description="Pipeline depth preset: fast_screen, balanced, or deep_rigor"
    )
    dry_run: bool = Field(
        default=False,
        description="When True, run pre-flight validation check at $0 cost without submitting live GPU jobs"
    )
    credentials: Optional[CredentialsConfig] = Field(
        default=None,
        description="Loaded authentication credentials for backends"
    )

    def to_dict(self) -> Dict[str, Any]:
        """Convert configuration to dictionary."""
        return self.model_dump()

    @classmethod
    def from_dict(cls, data: Dict[str, Any]) -> PipelineConfig:
        """Create PipelineConfig from dictionary, handling nested models."""
        routes_data = data.get("compute_routes", {})
        cutoffs_data = data.get("gate_cutoffs", {})
        creds_data = data.get("credentials")
        
        if isinstance(routes_data, dict):
            clean_routes = {}
            for k, v in routes_data.items():
                if isinstance(v, str):
                    try:
                        clean_routes[k] = BackendChoice(v.strip().lower())
                    except ValueError:
                        clean_routes[k] = v
                else:
                    clean_routes[k] = v
            routes = ComputeRouteConfig(**clean_routes)
        else:
            routes = routes_data
        cutoffs = (
            GateCutoffsConfig(**cutoffs_data)
            if isinstance(cutoffs_data, dict)
            else cutoffs_data
        )
        creds = (
            CredentialsConfig(**creds_data)
            if isinstance(creds_data, dict) and creds_data
            else creds_data
        )

        return cls(
            compute_routes=routes or ComputeRouteConfig(),
            gate_cutoffs=cutoffs or GateCutoffsConfig(),
            thinking_budget=ThinkingBudget(data.get("thinking_budget", "medium")),
            design_mode=DesignMode(data.get("design_mode", "balanced")),
            dry_run=bool(data.get("dry_run", True)),
            credentials=creds
        )
