"""agvarpred_core: internal reusable feature-engineering library for AGVarPred."""

__version__ = "1.0.0"

from .af_source import AFSource, LocalGnomADSource, OnlineAFSource, NoAFSource, resolve_af_source
from .feature_generator import FeatureGenerator
from .feature_selector import FeatureSelector
from .preprocessing import clean_name, encode_vep_features

__all__ = [
    "AFSource",
    "LocalGnomADSource",
    "OnlineAFSource",
    "NoAFSource",
    "resolve_af_source",
    "FeatureGenerator",
    "FeatureSelector",
    "clean_name",
    "encode_vep_features",
]
