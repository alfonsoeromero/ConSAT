"""
Overrepresentation analysis of Gene Ontology terms
"""
import sys
from collections import Counter
from math import exp, log
from operator import itemgetter

from gfam.utilities.bidict import Bidict

__author__ = "Tamas Nepusz"
__email__ = "tamas@cs.rhul.ac.uk"
__copyright__ = "Copyright (c) 2010, Tamas Nepusz"
__license__ = "MIT"

__all__ = ["OverrepresentationAnalyser"]

try:
    from scipy.special import gammaln
except ImportError:
    def gammaln(n):
        """Logarithm of Euler's gamma function for discrete values."""
        if n < 1:
            return float('inf')
        if n < 3:
            return 0.0
        c = [76.18009172947146, -86.50532032941677,
             24.01409824083091, -1.231739572450155,
             0.001208650973866179, -0.5395239384953 * 0.00001]
        x, y = float(n), float(n)
        tm = x + 5.5
        tm -= (x + 0.5) * log(tm)
        se = 1.0000000000000190015
        for j in range(6):
            y += 1.0
            se += c[j] / y
        return -tm + log(2.5066282746310005 * se / x)


def memoize(function):
    cache = {}

    def decorated_function(*args):
        if args in cache:
            return cache[args]
        else:
            val = function(*args)
            cache[args] = val
            return val
    return decorated_function


def logchoose(n, k):
    """ Computes the logarithm of the binomial
    coefficient, using the row-symmetry property
    """
    if n-k < k:
        return _logchoose(n, n-k)
    else:
        return _logchoose(n, k)


@memoize
def _logchoose(n, k):
    """Calculates the logarithm of n-choose-k"""
    lgn1 = gammaln(n+1)
    lgk1 = gammaln(k+1)
    lgnk1 = gammaln(n-k+1)
    return lgn1 - (lgnk1 + lgk1)


def hypergeom_pmf(k, M, n, N):
    """Hypergeometric probability moment function"""
    tot, good = M, n
    bad = tot - good
    return exp(logchoose(good, k) + logchoose(bad, N-k) - logchoose(tot, N))


def hypergeom_sf(k, M, n, N):
    """Tail distribution of a hypergeometric distribution. This is
    used to calculate p-values in the overrepresentation analysis
    tests"""
    tot, good = M, n
    bad = tot - good
    den = logchoose(tot, N)
    result = sum([exp(logchoose(good, x) + logchoose(bad, N-x) - den)
                  for x in range(k, N+1)])
    return result


def binomial_sf(k, M, n, N):
    """Tail distribution of a binomial distribution, used to compute
    p-value in the overrepresentation analysis test.

        ``k``: frequency of the GO term on the sample
        ``M``: number of objects in the BD
        ``n``: number of objects mapped to that GO term in the BD
        ``N``: sample size

    Translated to the "binomial" language, the number of successes
    is ``k`` (of ``N`` trials). The prior probability of getting a
    hit is n/M. We compute the probability of getting at least k
    hits in a sample of size N
    """
    result = 0.
    p = float(n/M)
    logp = log(p)
    log1_min_p = log(1.0-p)
    for x in range(k, N+1):
        result += exp(logchoose(N, x) + float(x)*logp + float(N-x)*log1_min_p)
    return result


class Correction(object):
    none, off = 0, 0
    fdr = 1
    bonferroni = 2
    sidak = 3


class OverrepresentationAnalyser(object):
    """Performs overrepresentation analysis of Gene Ontology
    terms on sets of entities that are annotated by some
    GO terms."""

    def __init__(self, tree, mapping, confidence=0.05, min_count=5,
                 correction="fdr", test="hypergeometric"):
        """Initializes the overrepresentation analysis algorithm by associating
        it to a given Gene Ontology tree and a given mapping from entities to
        their respective GO terms.

        `tree` must be an instance of `gfam.go.Tree`. `mapping` must be
        a bidirectional dictionary object (`gfam.utilities.bidict.Bidict`) that
        maps entities to GO terms and vice versa. For `mapping`, if an entity
        is annotated by a GO term, it is not necessary to list all the
        ancestors of that GO term for that entity, this will be taken care of
        by the class itself which always works on the copy of the given
        mapping.

        `confidence` is the confidence level of the test.

        `min_count` specifies which GO terms are to be excluded from the
        overrepresentation analysis; if a GO term occurs less than `min_count`
        times in `mapping`, it will not be considered.

        `correction` specifies the multiple hypothesis testing correction to be
        used and it must be one of the following:

        - ``None``, ``"none"`` or ``"off"``: no correction

        - ``"fdr"``: controlling the false discovery rate according
          to the method of Benjamini and Hochberg

        - ``"bonferroni"``: Bonferroni correction of the family-wise
          error rate

        - ``"sidak"``: Sidak correction of the family-wise error rate

        `test` specifies whether the test performed is an hypergeometric one
        (without replacement) or a binomial (with replacement)
        """
        self.tree = tree
        self.mapping = self._propagate_go_term_ancestors(mapping)
        self.confidence = float(confidence)
        self.min_count = max(1, int(min_count))
        cor = correction.lower()
        if cor in Correction.__dict__:
            self.correction = Correction.__dict__[cor]

        if test == "hypergeometric":
            self.test_func = hypergeom_sf
        elif test == "binomial":
            self.test_func = binomial_sf
        else:
            print("Test must be either 'hypergeometric' or 'binomial'")
            sys.exit(-1)

    def _propagate_go_term_ancestors(self, mapping):
        """Given a mapping object which maps entities to GO terms, this
        method ensures that all the ancestor terms of each GO term
        appear at the respective entities. Returns a copy of the mapping
        which contains these modifications."""
        result = Bidict(mapping)
        for entity, terms in result.iteritems_left():
            ancestors = self.tree.ancestors(*terms)
            result.add_left_multi(entity, ancestors)
        return result

    def enrichment_p(self, term_or_id, count, group_size):
        """Calculates the enrichment p-score of the given GO term or ID
        (`term_or_id`) if it occurs `count` times in a group of size
        `group_size`.
        """
        term = self.tree.ensure_term(term_or_id)
        objs = self.mapping.right[term]
        return self.test_func(count, self.mapping.len_left(),
                              len(objs), group_size)

    def test_counts(self, counts, group_size):
        """Given a dict that maps Gene Ontology terms to their
        occurrence counts and the number of entities in the group
        from which these term counts originate, calculates a
        list of overrepresented Gene Ontology terms."""
        confidence = self.confidence
        correction = self.correction
        min_count = self.min_count

        # Determine how many tests will be performed
        if correction:
            num_tests = sum(1
                            for term, count in counts.items()
                            if len(self.mapping.right[term]) >= min_count)
            if num_tests == 0:
                return []

        # If we are doing Bonferroni or Sidak correction, adjust the
        # confidence level
        if correction == Correction.bonferroni:
            confidence /= num_tests
        elif correction == Correction.sidak:
            confidence = 1 - (1. - confidence) ** (1. / num_tests)

        # Do the testing
        result = [(term, self.enrichment_p(term, count, group_size))
                  for term, count in counts.items()
                  if len(self.mapping.right[term]) >= min_count]

        # Filter the results
        if correction == Correction.fdr:
            result.sort(key=itemgetter(1))
            num_tests = float(num_tests)
            for k in range(len(result)):
                if result[k][1] > confidence * ((k+1) / num_tests):
                    result = result[0:k]
                    break
        else:
            result = [item for item in result if item[1] <= confidence]
            result.sort(key=itemgetter(1))
            if correction == Correction.bonferroni:
                result = [(c, p * num_tests) for c, p in result]
            elif correction == Correction.sidak:
                result = [(c, 1 - (1 - p) ** num_tests) for c, p in result]

        return result

    def test_group(self, group):
        """Overrepresentation analysis of the given group of objects.
        `group` must be an iterable yielding objects that are in
        `self.mapping.left`.
        """
        if not group:
            return []
        counts = Counter([go_term for item in group
                          for go_term in self.mapping.left.get(item, [])])
        return self.test_counts(counts, len(group))
