from gfam.modula.storage import DiskStorageEngine, NotFoundError


class GFamDiskStorageEngine(DiskStorageEngine):
    """Disk storage engine for GFam that has an empty `store` method.  This is
    because `GFamCalculation` writes directly into the results file and always
    returns ``None``, so there is no need to store the results explicitly.
    """

    def get_filename(self, source_or_module):
        """Retrieves the filename corresponding to a data source."""
        try:
            module = self.module_manager.get(source_or_module)
            if hasattr(module, "filename"):
                return module.filename
            return self._get_module_result_filename(module)
        except Exception as ex:
            raise NotFoundError(source_or_module, source_or_module, ex)

    def store(self, module, result):
        """Empty, does nothing"""
        pass
