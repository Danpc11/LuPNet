"""LuPNet - Lung Perfusion Network model.

Two symmetric binary trees (pulmonary arterial/venous and airway) that share
their N = 2^G terminal lung units. Units are lost (fibrotic or honeycomb) with a
hazard that depends on their tidal strain; after every step the vascular tree
remodels by a local shear set-point rule (shared target law of InFlow with law="shared"), the
breathing pattern is re-chosen by minimum mechanical power, and O2 exchange is
solved unit by unit.
"""
__version__ = "0.4.0"
from .params import Params
from .lung import Lung

__all__ = ["Params", "Lung", "__version__"]
