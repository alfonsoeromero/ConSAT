from gfam.assignment_utils.overlap_type import OverlapType


class AssignmentOverlapChecker:
    """Static class that contains the central logic of determining
    whether an assignment can be added to a partially assigned
    `SequenceWithAssignments`.

    The class has a class variable named `max_overlap` which stores
    the maximum allowed overlap size. This is 20 by default.
    """
    #: The maximum allowed overlap size.
    max_overlap = 20

    #: This is the minimum number of residues that a domain with the same
    #: InterPro id can leave uncovered in the parent until it starts.
    #: Otherwise it is discarded and considered as the same domain. Therefore,
    #: if two domains with the same InterPro are inserted and the beginning of
    #: both has a difference between 0 and 20 residues, the insertion is
    #: labeled as a synonym insertion
    min_parent_inserted_size = 20

    #: Here we give a scale of priorities for the different overlap types.
    #: When checking an assignment against all previously added assignments
    #: we will first check if there is overlap of the first kind
    #: (`OverlapType.OVERLAP`) with any of the already added members.
    #: If so, this overlap type is returned. If not, we pass on to the next
    #: `OverlapType`, until all the types have been tested. This is done to
    #: avoid unwanted behaviour. For instance, an assignment C might be found
    #: as an insertion in a large assignment A, which also contains an
    #: assignment B already inserted, overlapping with C. If we checked only
    #: if there were an insertion, A would be accepted, as no further checks
    #: are done.
    priority = [OverlapType.OVERLAP,
                OverlapType.DIFFERENT,
                OverlapType.SYNONYM_INSERTION,
                OverlapType.INSERTION_DIFFERENT,
                OverlapType.INSERTION,
                OverlapType.NO_OVERLAP]

    @classmethod
    def check(cls, sequence, assignment, interpro):
        """Checks whether an `assignment` can be added to a partially
        assigned `sequence`. `sequence` must be an instance of
        `SequenceWithAssignments`, `assignment` must be an instance
        of `Assignment`. InterPro is an interpro data structure
        (we use it to check for similar IPR domains).

        The most strict type of overlap (following the order defined
        in "priority" is returned)
        """
        overlaps_found = set(cls.check_single(assignment, other_assignment,
                                              interpro)
                             for other_assignment in sequence.assignments)
        for overlap_type in cls.priority:
            if overlap_type in overlaps_found:
                return overlap_type
        return OverlapType.NO_OVERLAP

    @classmethod
    def similar(cls, ipr1, ipr2, interpro):
        if ipr1 == ipr2:
            return True
        else:
            return interpro.tree.get_most_remote_ancestor(ipr1) ==\
                interpro.tree.get_most_remote_ancestor(ipr2)

    @classmethod
    def check_single(cls, assignment, other_assignment, interpro=None):
        """Checks whether the given `assignment` overlaps with another
        assignment `other_assignment`. Returns one of the following:

        - `OverlapType.NO_OVERLAP`: there is no overlap between the
          two given assignments

        - `OverlapType.SYNONYM_INSERTION`: `assignment` is inserted into
          `other_assignment` or vice versa, but they share the same InterPro
           term, and very few amino acids in the parent (less than
           `min_parent_inserted_size`) are covered in the starting fragment.
           Therefore we might think it is not a real insertion, but two
           expressions of the same domain (belonging for example to two
           different signature databases).

        - `OverlapType.INSERTION`: there is a complete domain insertion in
          either direction

        - `OverlapType.INSERTION_DIFFERENT`: `assignment` is inserted into
          `other_assignment` or vice versa, but they have different data
          sources.

        - `OverlapType.DIFFERENT`: `other_assignment` overlaps with
        `assignment` partially, but they have different data sources

        - `OverlapType.OVERLAP`: `other_assignment` overlaps with `assignment`
          partially, they have the same data source, but the size of
          the overlap is larger than the maximum allowed overlap specified
          in `AssignmentOverlapChecker.max_overlap`.
        """
        start, end = assignment.start, assignment.end
        other_start, other_end = other_assignment.start, other_assignment.end

        if ((other_start <= start and other_end >= end)
                or (other_start >= start and other_end <= end)):
            # Domains are one inside the other: could be INSERTION,
            # INSERTION_DIFFERENT,
            # SYNONYM_INSERTION or OVERLAP
            similar_iprs = (cls.similar(assignment.interpro_id,
                                        other_assignment.interpro_id,
                                        interpro) or
                            assignment.domain == other_assignment.domain)

            if abs(other_start-start) + abs(other_end-end)\
               < cls.min_parent_inserted_size:
                if similar_iprs and\
                   (abs(other_start - start) < cls.min_parent_inserted_size):
                    return OverlapType.SYNONYM_INSERTION
                else:
                    return OverlapType.OVERLAP
            elif other_assignment.source == assignment.source:
                # This is a valid domain insertion, other_assignment is
                # inserted into assignment
                if not similar_iprs or\
                   (abs(other_start - start) >= cls.min_parent_inserted_size):
                    return OverlapType.INSERTION
                else:
                    return OverlapType.OVERLAP
            else:
                return OverlapType.INSERTION_DIFFERENT

        if (other_start <= start and start <= other_end <= end) or\
           (start <= other_start and other_start <= end <= other_end):
            # Domains could overlap (either A B A' B' or B A B' A')
            if other_assignment.source == assignment.source:
                # This is a partial overlap
                overlap_size = min(other_end - start, end - other_start) + 1
                if overlap_size > cls.max_overlap:
                    return OverlapType.OVERLAP
            else:
                return OverlapType.DIFFERENT

        return OverlapType.NO_OVERLAP
