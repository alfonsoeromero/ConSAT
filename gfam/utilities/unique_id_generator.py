class UniqueIdGenerator(object):
    """A dictionary-like class that can be used to assign unique integer IDs to
    names.

    Usage:

    >>> gen = UniqueIdGenerator()
    >>> gen["A"]
    0
    >>> gen["B"]
    1
    >>> gen["C"]
    2
    >>> gen["A"]      # Retrieving already existing ID
    0
    >>> len(gen)      # Number of already used IDs
    3
    """

    def __init__(self, id_generator=None):
        """Creates a new unique ID generator. `id_generator` specifies how do we
        assign new IDs to elements that do not have an ID yet. If it is `None`,
        elements will be assigned integer identifiers starting from 0. If it is
        an integer, elements will be assigned identifiers starting from the
        given
        integer. If it is an iterator or generator, its `next` method will be
        called every time a new ID is needed."""
        if id_generator is None:
            id_generator = 0
        if isinstance(id_generator, int):
            import itertools
            self._generator = itertools.count(id_generator)
        else:
            self._generator = id_generator
        self._ids = {}

    def __getitem__(self, item):
        """Retrieves the ID corresponding to `item`. Generates a new ID for
        `item` if it is the first time we request an ID for it."""
        try:
            return self._ids[item]
        except KeyError:
            self._ids[item] = next(self._generator)
            return self._ids[item]

    def __len__(self):
        """Retrieves the number of added elements
        in this `UniqueIDGenerator`"""
        return len(self._ids)

    def values(self):
        """Returns the list of items added so far. Items are ordered according
        to the standard sorting order of their keys, so the values will be
        exactly in the same order they were added if the ID generator generates
        IDs in ascending order. This hold, for instance, to numeric ID
        generators that assign integers starting from a given number."""
        return sorted(self._ids.keys(), key=self._ids.__getitem__)
