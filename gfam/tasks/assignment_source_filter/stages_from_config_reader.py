import re
from typing import List, Set, Union

from gfam.utils import complementerset


class StagesFromConfigReader:
    def __init__(self, parser):
        self.parser = parser
        self._stages_read: bool = False
        self._stages_from_config: List[Union[Set[str], complementerset]]

    def get_stages_from_config(self) -> List[Union[Set[str], complementerset]]:
        if not self._stages_read:
            self._stages_from_config = self._get_stages_from_config()
        return self._stages_from_config

    def _get_stages_from_config(self) -> List[Union[Set[str],
                                                    complementerset]]:
        """Turns to the configuration file specified at startup to
        fetch the data sources to be used in each stage of the algorithm.
        If there is no configuration file specified or it does not
        contain the corresponding keys, it will simply use a default
        stage setup which ignores HMMPanther and Gene3D in the first
        and second steps, but uses all sources in the third step.

        The method will be looking for configuration keys named like
        ``stages.1``, ``stages.2`` and so on in the ``analysis:iprscan_filter``
        section of the config file. The value of each such config key must
        be an expression consisting of assignment source names and the
        operators ``+`` and ``-``, with their usual meaning of addition
        and exclusion. The special source name ``ALL`` means all possible
        data sources, enabling us to write expressions like ``ALL-HMMPanther``
        (meaning all the sources except HMMPanther). Some examples:

        - ``HMMPanther`` means HMMPanther only.
        - ``ALL`` means all possible data sources.
        - ``HMMPanther+HMMPfam`` means HMMPanther or HMMPfam.
        - ``ALL-HMMPanther-Gene3D`` means all possible data sources but
          HMMPanther or Gene3D.
        - ``ALL+HMMPanther`` does not really make sense as you are extending
          all data sources with HMMPanther, so it is equivalent to ``ALL``.
          GFam will figure out what you meant anyway.
        """
        cfg = self.parser.config
        if cfg is None:
            spec = ["ALL-HMMPanther-Gene3D", "ALL-HMMPanther-Gene3D",
                    "ALL"]
        else:
            spec, idx = [], 1
            section = "analysis:iprscan_filter"
            while cfg.has_option(section, "stages.%d" % idx):
                spec.append(cfg.get(section, "stages.%d" % idx))
                idx += 1

        regexp = re.compile(r"([-+])?\s*([^-+]+)")
        result = []
        for item in spec:
            sources: Union[Set[str], complementerset] = set()
            for match in regexp.finditer(item):
                sign, source = match.groups()
                s_source: Union[Set[str], complementerset]
                if source == "ALL":
                    s_source = complementerset()
                else:
                    s_source = set([source.strip()])
                if sign == "-":
                    sources -= s_source
                else:
                    sources |= s_source
            result.append(sources)

        return result
