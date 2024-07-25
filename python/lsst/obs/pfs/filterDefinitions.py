from lsst.obs.base import FilterDefinition, FilterDefinitionCollection

__all__ = ("pfsFilterDefinitions",)


pfsFilterDefinitions = FilterDefinitionCollection(
    FilterDefinition("b", "b", "PFS blue arm"),
    FilterDefinition("r", "r", "PFS red arm"),
    FilterDefinition("m", "m", "PFS medium-resolution red arm"),
    FilterDefinition("n", "n", "PFS near-infrared arm"),
)
