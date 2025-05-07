from pfs.datamodel import PfsConfig

__all__ = ("PfsConfigEditor",)


class PfsConfigEditor:
    """Class to edit PfsConfig objects at read-time."""

    def __call__(self, pfsConfig: PfsConfig) -> "PfsConfig":
        """Edit the PfsConfig object.

        Parameters
        ----------
        pfsConfig : `PfsConfig`
            Fiber configuration to edit (modified).

        Returns
        -------
        pfsConfig : `PfsConfig`
            Edited fiber configuration.
        """
        self.updateFiberStatus(pfsConfig)
        return pfsConfig

    def updateFiberStatus(self, pfsConfig: PfsConfig) -> None:
        """Update the fiber status for the PfsConfig object.

        Set ``BAD_PSF`` for certain fibers.

        Parameters
        ----------
        pfsConfig : `PfsConfig`
            Fiber configuration to edit.
        """
        pass
