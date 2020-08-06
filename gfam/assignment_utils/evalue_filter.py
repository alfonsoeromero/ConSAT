class EValueFilter(object):
    """Given an `Assignment`, this filter tells whether the assignment's
    E-value is satisfactory to accept it.

    The filter supports different E-values for different data sources. By
    default, the E-value threshold is infinity for all data sources.
    """

    def __init__(self):
        self.default_e_value = float('inf')
        self.thresholds = {}

    def set_threshold(self, source, evalue):
        """Sets the E-value threshold for the given data source."""
        self.thresholds[source] = evalue

    def is_acceptable(self, assignment):
        """Checks whether the given assignment is acceptable.

        This method looks up the E-value threshold corresponding to the
        ``source`` of the assignment and returns ``True`` if the E-value
        of the assignment is less than the threshold, ``False``
        otherwise."""
        threshold = self.thresholds.get(assignment.source,
                                        self.default_e_value)
        return assignment.evalue <= threshold

    @classmethod
    def from_string(cls, description):
        """Constructs an E-value filter from a string description that
        can be used in command line arguments and configuration files.

        The string description is a semicolon-separated list of
        source-threshold pairs. For instance, the following is a valid
        description giving an E-value of 0.001 for HMMPfam sources,
        0.005 for HMMSmart sources and 0.007 for everything else::

            HMMPfam=0.001;HMMSmart=0.005;0.007

        The last entry denotes the default E-value; in particular, if a
        number is not preceded by a source name, it is assumed to be a
        default E-value. If no default E-value is given, infinity will
        be used.
        """
        result = cls()
        for part in description.split(";"):
            part = part.strip()
            if "=" in part:
                source, evalue = part.split("=", 1)
                evalue = float(evalue)
                result.set_threshold(source, evalue)
            else:
                result.default_e_value = float(part)
        return result
