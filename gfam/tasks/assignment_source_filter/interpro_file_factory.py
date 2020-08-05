from typing import Optional

from gfam.interpro import InterPro


class InterproFileFactory:
    @classmethod
    def get_from_file(cls, interpro_file: Optional[str] = None) -> InterPro:
        if interpro_file is None:
            return InterPro()
        else:
            return InterPro.from_file(interpro_file)
