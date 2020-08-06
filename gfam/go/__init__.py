#!/usr/bin/env python
"""
A higher level Gene Ontology representation in Python
"""
from collections import deque

import gfam.go.obo

__author__ = "Tamas Nepusz"
__email__ = "tamas@cs.rhul.ac.uk"
__copyright__ = "Copyright (c) 2009, Tamas Nepusz"
__license__ = "MIT"
__version__ = "0.1"

__all__ = ["Tree", "Term"]


class Tree(object):
    """Class representing the GO tree. A GO tree contains many GO terms
    represented by `Term` objects.
    """

    def __init__(self):
        self.terms = {}
        self.aliases = {}

    def add(self, term):
        """Adds a `Term` to this GO tree"""
        self.terms[term.id] = term

    def add_alias(self, canonical, alias):
        """Adds an alias to the given canonical term in the GO tree"""
        self.aliases[alias] = canonical

    def lookup(self, identifier):
        """Looks up a `Term` in this tree by ID. Also cares about alternative
        IDs"""
        try:
            return self.terms[identifier]
        except KeyError:
            return self.terms[self.aliases[identifier]]

    def ensure_term(self, term_or_id):
        """Given a `Term` or a GO term ID, returns an object
        that is surely a `Term`"""
        if isinstance(term_or_id, Term):
            return term_or_id
        return self.lookup(term_or_id)

    def ancestors(self, *args):
        """Returns all the ancestors of a given `Term`
        (or multiple terms) in this tree. The result is a
        list of `Term` instances."""
        unprocessed_terms = deque(self.ensure_term(term_or_id)
                                  for term_or_id in args)
        result = set()
        while unprocessed_terms:
            term = unprocessed_terms.popleft()
            result.add(term)
            unprocessed_terms.extend(self.parents(term))
        return result

    def parents(self, term_or_id):
        """Returns the direct parents of a `Term` in this tree.
        `term_or_id` can be a GO term ID or a `Term`.
        The result is a list of `Term` instances."""
        term = self.ensure_term(term_or_id)
        parent_ids = term.tags.get("is_a", [])
        return [self.lookup(identifier.value) for identifier in parent_ids]

    def paths_to_root(self, *args):
        """Finds all the paths from a term (or multiple terms if multiple
        arguments are used) to the root."""
        terms = [self.ensure_term(term_or_id) for term_or_id in args]
        queue = deque([(term, [term]) for term in terms])
        while queue:
            first, path = queue.popleft()
            parents = self.parents(first)
            if not parents:
                yield path
            for parent in self.parents(first):
                queue.append((parent, path + [parent]))

    def to_igraph(self, rel="is_a"):
        """Returns an :mod:`igraph` graph representing this GO tree. This is
        handy if you happen to use igraph_.

        .. _igraph: http://igraph.sf.net
        """
        import igraph
        graph = igraph.Graph(n=len(self.terms), directed=True)
        graph.vs["id"] = self.terms.keys()
        graph.vs["name"] = [term.name for term in self.terms.values()]

        term_id_to_idx = dict(zip(self.terms.keys(),
                                  range(len(self.terms))))
        edgelist = []
        for identifier, term in self.terms.items():
            source = term_id_to_idx[identifier]
            for parent_id in term.tags.get(rel, []):
                target = term_id_to_idx.get(parent_id.value, None)
                if target is None:
                    continue
                edgelist.append((source, target))

        graph.add_edges(edgelist)
        return graph

    @classmethod
    def from_obo(cls, f_pointer):
        """Constructs a GO tree from an OBO file. `f_pointer` is a file pointer
        to the OBO file we want to use"""
        parser = gfam.go.obo.Parser(f_pointer)
        tree = cls()
        for stanza in parser:
            term = Term.from_stanza(stanza)
            tree.add(term)
            for alt_id in stanza.tags.get("alt_id", []):
                tree.add_alias(term.id, alt_id.value)
        return tree

    def __len__(self):
        return len(self.terms)


class Term(object):
    """Class representing a single GO term"""

    __slots__ = ("id", "name", "tags")

    def __init__(self, go_id, name="", tags=None):
        """Constructs a GO term with the given go_id, the given human-
        readable name and the given tags."""
        self.id = str(go_id)
        self.name = str(name)
        if tags:
            self.tags = dict(tags)
        else:
            self.tags = {}

    def __repr__(self):
        """String representation of a GO term"""
        return "%s(%r, %r, %r)" % (self.__class__.__name__,
                                   self.id, self.name, self.tags)

    def __str__(self):
        """Returns just the ID of the GO term"""
        return self.id

    @classmethod
    def from_stanza(cls, stanza):
        """Constructs a GO term from a stanza coming from an
        OBO file. `stanza` must be an instance of `gfam.go.obo.Stanza`.
        """
        identifier = stanza.tags["id"][0]
        name = stanza.tags.get("name", [id])[0]
        return cls(identifier, name, stanza.tags)


def test():
    from time import time

    start = time()
    tree = Tree.from_obo(open("gene_ontology.obo"))
    end = time()

    print("{} terms in GO tree parsed in {:.2f} seconds\n".format(
        len(tree), end-start))

    print("All paths from GO:0009651 to root:")
    for path in tree.paths_to_root("GO:0009651"):
        print(" -> ".join(node.name for node in path))
    print()

    print("Ancestors of GO:0009651:")
    print("\n".join(str(term) for term in tree.ancestors("GO:0009651")))
    print()


if __name__ == "__main__":
    import sys
    sys.exit(test())
