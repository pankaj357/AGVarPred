"""AGVarPred: pathogenicity prediction from germline variants."""

__version__ = "1.0.0"

from .predictor import AGVarPredAutoPredictor, AGVarPredPredictor

__all__ = ["AGVarPredPredictor", "AGVarPredAutoPredictor"]
