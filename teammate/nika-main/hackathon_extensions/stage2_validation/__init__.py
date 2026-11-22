"""Stage 2 Validation Package.

Implements the 4-gate multi-scale biophysical and computational validation pipeline:
- Gate 1: Co-folding & Binding Affinity (Boltz-2 / AlphaFold3 / Protenix-v1)
- Gate 2: Active-Site Preorganization (PLACER 100 ensemble samples)
- Gate 3: Stability & Sequence Fitness (ThermoMPNN & ESM3 / SaProt)
- Gate 4: Dynamics & Transition-State Energetics (OpenMM MD, FieldTools, QM/MM via ORCA)
"""

from hackathon_extensions.stage2_validation.scorecard import DesignScorecard
from hackathon_extensions.stage2_validation.gate1_cofolding import Gate1CofoldingEvaluator
from hackathon_extensions.stage2_validation.gate2_preorganization import Gate2PreorganizationEvaluator
from hackathon_extensions.stage2_validation.gate3_stability import Gate3StabilityEvaluator
from hackathon_extensions.stage2_validation.gate4_dynamics_qm import Gate4DynamicsQMEvaluator
from hackathon_extensions.stage2_validation.pipeline import Stage2ValidationPipeline, PipelineConfig

__all__ = [
    "DesignScorecard",
    "Gate1CofoldingEvaluator",
    "Gate2PreorganizationEvaluator",
    "Gate3StabilityEvaluator",
    "Gate4DynamicsQMEvaluator",
    "Stage2ValidationPipeline",
    "PipelineConfig",
]
