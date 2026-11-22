"""Stage 1 Generative De Novo Enzyme Design Module.

Exposes primary generative design adapters and execution runner:
- RFD4ProteinaAdapter: 32D force-field-guided joint sequence & all-atom structure codesign.
- RFD3Adapter: Small-molecule & catalytic motif scaffolding with automated LigandMPNN sequence design.
- LigandMPNNAdapter: Active-site-aware sequence design with catalytic residue protection.
- Stage1Runner: Standalone orchestrator with VRAM safety checks and manifest output.
"""

from hackathon_extensions.stage1_design.base_engine import (
    BaseDesignEngine,
    CandidateDesign,
)
from hackathon_extensions.stage1_design.ligand_mpnn_adapter import (
    LigandMPNNAdapter,
    LigandMPNNConfig,
)
from hackathon_extensions.stage1_design.rfd3_adapter import (
    RFD3Adapter,
    RFD3Config,
)
from hackathon_extensions.stage1_design.rfd4_proteina_adapter import (
    THINKING_BUDGET_PRESETS,
    RFD4LatentGuidanceConfig,
    RFD4ProteinaAdapter,
    RFD4ThinkingBudgetPreset,
)
from hackathon_extensions.stage1_design.runner import (
    Stage1Runner,
)

__all__ = [
    "BaseDesignEngine",
    "CandidateDesign",
    "RFD4ProteinaAdapter",
    "RFD4LatentGuidanceConfig",
    "RFD4ThinkingBudgetPreset",
    "THINKING_BUDGET_PRESETS",
    "RFD3Adapter",
    "RFD3Config",
    "LigandMPNNAdapter",
    "LigandMPNNConfig",
    "Stage1Runner",
]
