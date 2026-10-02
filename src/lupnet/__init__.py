"""LuPNet - Lung Perfusion-Ventilation Network prototype.

Two symmetric binary trees (pulmonary arterial/venous and airway) that share
their N = 2^G terminal lung units. Units are lost (fibrotic or honeycomb) with a
hazard that depends on their tidal strain; after every step the vascular tree
remodels by the local shear set-point rule of vascular-shear-setpoint, the
breathing pattern is re-chosen by minimum mechanical power, and O2 exchange is
solved unit by unit.
"""
from .params import Params
from .lung import Lung
