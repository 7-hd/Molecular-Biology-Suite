from rich.console import Console
from rich.panel import Panel
from rich.prompt import Prompt, Confirm
from rich.table import Table
from rich.align import Align
from rich import box

from utilities import get_sequence_from_user, truncate, CRISPR_SYSTEMS, CODON_TABLES
from biology import (
    PrimerDesigner, CRISPRDesigner, CRISPRAnimator, RestrictionMapper,
    CodonOptimizer, PlasmidDesigner, siRNADesigner, RNAStructurePredictor,
    PhylogeneticsBuilder, VaccineEpitopeDesigner, AMPDesigner, ARGScanner,
)
