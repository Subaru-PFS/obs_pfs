"""Model registry — lookup table from MODEL string (FITS header) to model class."""

from __future__ import annotations

from .base import BlockFitResult, Model
from .polynomial import PolynomialModel

MODEL_REGISTRY: dict[str, type[Model]] = {}


def registerModel(modelClass: type[Model], *, overwrite: bool = False) -> None:
    """Register a model class under its ``modelName`` attribute."""
    name = modelClass.modelName  # type: ignore[attr-defined]
    if name in MODEL_REGISTRY and not overwrite:
        raise ValueError(
            f"model name {name!r} already registered; "
            f"pass overwrite=True to replace"
        )
    MODEL_REGISTRY[name] = modelClass


# Register the built-in model.
registerModel(PolynomialModel)

__all__ = [
    "MODEL_REGISTRY",
    "Model",
    "BlockFitResult",
    "PolynomialModel",
    "registerModel",
]
