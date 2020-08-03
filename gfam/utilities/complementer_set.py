class ComplementerSet:
    """This object behaves more or less like a set, with one exception,
    the membership checking. For a `complementerset` object, you can
    define the elements which are *not* in the set, everything else is
    contained in it. The semantics of the operators are the same as for
    sets.

    Usage example::

        >>> s = complementerset()
        >>> "abc" in s
        True
        >>> s in s
        True
        >>> s -= set(["abc"])
        >>> s
        complementerset(['abc'])
        >>> "abc" in s
        False
    """

    __slots__ = ("_set", )

    def __init__(self, iterable=()):
        """Constructs a complementer set that contains everything except
        the members of the given iterable."""
        self._set = set(iterable)

    def difference_update(self, *args):
        """Removes all elements of another set from this set.
        """
        self._set.update(*args)

    def discard(self, member):
        """Removes an element from the complementer set if it is a member.
        """
        self._set.add(member)

    def iterexcluded(self):
        """Iterates over the items excluded from the complementer set.

        Example::

            >>> s = complementerset([5, 7, 4])
            >>> print sorted(list(s.iterexcluded()))
            [4, 5, 7]
        """
        return iter(self._set)

    def remove(self, member):
        """Removes an element from the complementer set; it must be a member.
        If the element is not a member, raises a ``KeyError``.

        Example::

            >>> s = complementerset()
            >>> s.remove(2)
            >>> print s
            complementerset([2])
            >>> s.remove(2)
            Traceback (most recent call last):
              File "<stdin>", line 4, in ?
            KeyError: 2
        """
        if member in self._set:
            raise KeyError(member)
        self._set.add(member)

    def __and__(self, other):
        """Example::

            >>> s = complementerset()
            >>> s2 = set([1,2,3])
            >>> (s & s2) == s2
            True
            >>> s = complementerset([3,4])
            >>> (s & s2) == set([1,2])
            True
            >>> (s & complementerset([1,2,3])) == complementerset([1,2,3,4])
            True
        """
        if isinstance(other, self.__class__):
            return ComplementerSet(self._set | other._set)

        self._ensure_set(other)
        return set(other) - self._set

    def __contains__(self, what):
        if hasattr(what, "__hash__"):
            return what not in self._set
        else:
            return True

    def __hash__(self):
        return id(self)

    def __eq__(self, other):
        return isinstance(other, self.__class__)\
            and self._set == other._set

    def __ge__(self, other):
        """Example::

            >>> s = complementerset()
            >>> s2 = complementerset()
            >>> s2 >= s
            True
            >>> s >= s2
            True
            >>> s >= complementerset([1,2])
            True
            >>> s >= set([1,2,3])
            True
            >>> s >= 2
            Traceback (most recent call last):
              File "<stdin>", line 1, in ?
            NotImplementedError
            >>> complementerset([1,2,3]) >= complementerset([1,2,3])
            True
            >>> complementerset([1,2]) >= complementerset([2,3,4])
            False
        """
        if isinstance(other, self.__class__):
            return self._set <= other._set

        self._ensure_set(other)
        return True

    def __gt__(self, other):
        """Example::

            >>> s = complementerset()
            >>> s2 = complementerset()
            >>> s2 > s
            False
            >>> s > s2
            False
            >>> s > complementerset([1,2])
            True
            >>> s > set([1,2,3])
            True
            >>> s > 2
            Traceback (most recent call last):
              File "<stdin>", line 1, in ?
            NotImplementedError
            >>> complementerset([1,2,3]) > complementerset([1,2,3])
            False
            >>> complementerset([1,2]) > complementerset([1,2,3])
            True
            >>> complementerset([1,2,3]) > complementerset([1,2])
            False
            >>> complementerset([1,2]) > complementerset([2,3,4])
            False
        """
        return self != other and self >= other

    def __iand__(self, other):
        """Example::

            >>> s = complementerset([3,4,5])
            >>> s &= s
            >>> s == complementerset([3, 4, 5])
            True
            >>> s &= complementerset([5,6])
            >>> s
            complementerset([3, 4, 5, 6])
            >>> s &= set([1,2,3])
            >>> s == set([1, 2])
            True
            >>> s = complementerset()
            >>> s &= 2
            Traceback (most recent call last):
              File "<stdin>", line 1, in ?
            NotImplementedError
        """
        if isinstance(other, self.__class__):
            self._set |= other._set
            return self

        self._ensure_set(other)
        return other - self._set

    def __ior__(self, other):
        """Example::

            >>> s = complementerset([3,4,5])
            >>> s |= s
            >>> s == complementerset([3, 4, 5])
            True
            >>> s |= set([1,2,3])
            >>> s == complementerset([4, 5])
            True
            >>> s |= complementerset([1, 4])
            >>> s
            complementerset([4])
            >>> s |= 2
            Traceback (most recent call last):
              File "<stdin>", line 1, in ?
            NotImplementedError
        """
        if isinstance(other, self.__class__):
            self._set &= other._set
            return self

        self._ensure_set(other)
        self._set -= other
        return self

    def __isub__(self, other):
        """Example::

            >>> s = complementerset()
            >>> s -= set([1,2,3])
            >>> s
            complementerset([1, 2, 3])
            >>> s -= 2
            Traceback (most recent call last):
              File "<stdin>", line 1, in ?
            NotImplementedError
            >>> s -= complementerset([1,2,4])
            >>> s
            set([4])
        """
        if isinstance(other, self.__class__):
            return other._set - self._set

        self._ensure_set(other)
        self._set |= other
        return self

    def __le__(self, other):
        """Example::

            >>> s = complementerset()
            >>> s2 = complementerset()
            >>> s2 <= s
            True
            >>> s <= s2
            True
            >>> s <= complementerset([1,2])
            False
            >>> s <= set([1,2,3])
            False
            >>> s <= 2
            Traceback (most recent call last):
              File "<stdin>", line 1, in ?
            NotImplementedError
            >>> complementerset([1,2,3]) <= complementerset([1,2,3])
            True
            >>> complementerset([1,2]) <= complementerset([2,3,4])
            False
        """
        if isinstance(other, self.__class__):
            return self._set >= other._set

        self._ensure_set(other)
        return False

    def __lt__(self, other):
        """Example::

            >>> s = complementerset()
            >>> s2 = complementerset()
            >>> s2 < s
            False
            >>> s < s2
            False
            >>> s < complementerset([1,2])
            False
            >>> complementerset([1,2]) < s
            True
            >>> s < set([1,2,3])
            False
            >>> s < 2
            Traceback (most recent call last):
              File "<stdin>", line 1, in ?
            NotImplementedError
            >>> complementerset([1,2,3]) < complementerset([1,2,3])
            False
            >>> complementerset([1,2]) < complementerset([1,2,3])
            False
            >>> complementerset([1,2,3]) < complementerset([1,2])
            True
            >>> complementerset([1,2]) < complementerset([2,3,4])
            False
        """
        return self != other and self <= other

    def __ne__(self, other):
        return not self == other

    def __or__(self, other):
        """Example::

            >>> s = complementerset([2,4,5])
            >>> s | set([1,2,3]) == complementerset([4, 5])
            True
            >>> s | complementerset([1,2])
            complementerset([2])
            >>> s | 2
            Traceback (most recent call last):
              File "<stdin>", line 1, in ?
            NotImplementedError
        """
        if isinstance(other, self.__class__):
            return ComplementerSet(self._set & other._set)

        self._ensure_set(other)
        return ComplementerSet(self._set - other)

    def __rand__(self, other):
        """Example::

            >>> s = complementerset([2,4,5])
            >>> set([1,2,3]) & s == set([1,3])
            True
            >>> complementerset([1,2,3]) & s == complementerset([1,2,3,4,5])
            True
            >>> 2 & s
            Traceback (most recent call last):
              File "<stdin>", line 1, in ?
            NotImplementedError
        """
        return self & other

    def __ror__(self, other):
        """Example::

            >>> s = complementerset([2,4,5])
            >>> set([1,2,3]) | s == complementerset([4, 5])
            True
            >>> complementerset([1, 2]) | s
            complementerset([2])
            >>> 2 | s
            Traceback (most recent call last):
              File "<stdin>", line 1, in ?
            NotImplementedError
        """
        return self | other

    def __repr__(self):
        return "%s(%r)" % (self.__class__.__name__, list(self._set))

    def __rsub__(self, other):
        """Example::

            >>> set([1,2,3]) - complementerset()
            set([])
            >>> (set([1,2,3]) - complementerset([1,2,4])) == set([1,2])
            True
            >>> 2 - complementerset()
            Traceback (most recent call last):
              File "<stdin>", line 1, in ?
            NotImplementedError
        """
        if isinstance(other, self.__class__):
            raise NotImplementedError

        self._ensure_set(other)
        return other.intersection(self._set)

    def __sub__(self, other):
        """Example::

            >>> complementerset() - set([1,2,3])
            complementerset([1, 2, 3])
            >>> complementerset([1,2,3]) - complementerset([2,4])
            set([4])
            >>> complementerset([3,4]) - set([1,2,3])
            complementerset([1, 2, 3, 4])
            >>> complementerset() - 2
            Traceback (most recent call last):
              File "<stdin>", line 1, in ?
            NotImplementedError
        """
        if isinstance(other, self.__class__):
            return other._set - self._set

        self._ensure_set(other)
        return ComplementerSet(self._set | other)

    @staticmethod
    def _ensure_set(obj):
        if not isinstance(obj, (set, frozenset, ComplementerSet)):
            raise NotImplementedError
