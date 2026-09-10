"""polyfind: efficient search for stable atomic arrangements in semi-crystalline polymers."""
from .polymers import PVDF, PE, Monomer, Polymer, copolymer, get_polymer
from .ris import RISModel

__all__ = ["PVDF", "PE", "Monomer", "Polymer", "copolymer", "get_polymer", "RISModel"]
__version__ = "0.1.0"
