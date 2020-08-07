class Bidict:
    """Bidirectional many-to-many mapping.

    This class models many-to-many mappings that are used in some places in
    GFam. For instance, the mapping of GO term identifiers to InterPro
    domain identifiers is typically a many-to-many mapping (domains are
    annotated by multiple GO terms, and GO terms may belong to multiple
    domains).

    Being a general many-to-many mapping class, instances contain two
    member dictionaries: `left` (which maps items from one side of the
    relationship to *sets* of items from the other side of the relationship)
    and `right` (which contains the exact opposite). You should only query
    these dictionaries directly, manipulation should be done by the methods
    provided by the `Bidict` class to ensure that the two dicts are kept
    in sync.

    Example::

        >>> bd = Bidict()
        >>> bd.add_left("foo", "bar")
        >>> bd.add_left("foo", "baz")
        >>> bd.get_left("foo") == set(['bar', 'baz'])
        True
        >>> bd.get_right("bar")
        set(['foo'])
        >>> bd.get_right("baz") == set(['foo', 'frob'])
        True
        >>> bd.len_left()
        2
    """

    def __init__(self, items=None):
        object.__init__(self)

        self.left = {}
        self.right = {}

        if isinstance(items, Bidict):
            # Copying an existing Bidict
            self.left = dict(items.left)
            self.right = dict(items.right)
        elif isinstance(items, dict):
            # items assigns an element from the left dict to a
            # set of elements from the right dict
            for key, values in items.items():
                self.add_left_multi(key, values)
        elif items is not None:
            raise TypeError("items must be dict or Bidict, got %r" %
                            type(items))

    def add_left(self, v1, v2):
        """Adds a pair of items `v1` and `v2` to the mapping s.t. `v1` as a
        left item is mapped to `v2` as a right item."""
        try:
            self.left[v1].add(v2)
        except KeyError:
            self.left[v1] = set([v2])
        try:
            self.right[v2].add(v1)
        except KeyError:
            self.right[v2] = set([v1])

    def add_left_multi(self, v1, v2s):
        """Associates multiple items in `v2s` to `v1` when `v1` is interpreted
        as a left item"""
        try:
            self.left[v1].update(v2s)
        except KeyError:
            self.left[v1] = set(v2s)
        for v2 in v2s:
            try:
                self.right[v2].add(v1)
            except KeyError:
                self.right[v2] = set([v1])

    def get_left(self, v1, default=None):
        """Returns the items associated to `v1` when `v1` is looked up from the
        left dictionary. `default` will be returned if `v1` is not in the left
        dictionary."""
        return self.left.get(v1, default)

    def get_right(self, v1, default=None):
        """Returns the items associated to `v1` when `v1` is looked up from the
        right dictionary. `default` will be returned if `v2` is not in the
        right dictionary."""
        return self.right.get(v1, default)

    def len_left(self):
        """Returns the number of unique left items"""
        return len(self.left)

    def iteritems_left(self):
        """Iterates over the left dictionary"""
        return self.left.items()
