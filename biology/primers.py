"""Domain-oriented biology modules for the suite.

This is the start of the architectural refactor. The legacy `structure.py`
remains the source of truth for now, while these modules provide a cleaner
import path for future modularization.
"""

from structure import (
    AMPDesigner,
    ARGScanner,
    CodonOptimizer,
    CRISPRAnimator,
    CRISPRDesigner,
    PhylogeneticsBuilder,
    PlasmidDesigner,
    PrimerDesigner,
    RestrictionMapper,
    RNAStructurePredictor,
    VaccineEpitopeDesigner,
    siRNADesigner,
)

__all__ = [
    "AMPDesigner",
    "ARGScanner",
    "CodonOptimizer",
    "CRISPRAnimator",
    "CRISPRDesigner",
    "PhylogeneticsBuilder",
    "PlasmidDesigner",
    "PrimerDesigner",
    "RestrictionMapper",
    "RNAStructurePredictor",
    "VaccineEpitopeDesigner",
    "siRNADesigner",
]
