from lsst.obs.base import FilterDefinition, FilterDefinitionCollection

__all__ = ("pfsFilterDefinitions",)


pfsFilterDefinitions = FilterDefinitionCollection(
    FilterDefinition("b_PFS", "b", "PFS blue arm"),
    FilterDefinition("r_PFS", "r", "PFS red arm"),
    FilterDefinition("m_PFS", "m", "PFS medium-resolution red arm"),
    FilterDefinition("n_PFS", "n", "PFS near-infrared arm"),
)
