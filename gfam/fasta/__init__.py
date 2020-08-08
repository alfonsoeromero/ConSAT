import re


def regexp_remapper(iterable, regexp=None, replacement=r'\g<id>'):
    """Regexp-based sequence ID remapper.

    This class takes a sequence iterator as an input and remaps
    each ID in the sequence using a call to `re.sub`.
    `iterable` must yield `SeqRecord` objects; `regexp`
    is the regular expression matching the part to be
    replaced, `replacement` is the replacement string that
    may contain backreferences to `regexp`.

    If `regexp` is ``None`` or an empty string, returns
    the original iterable.
    """
    if not regexp:
        for seq in iterable:
            yield seq
        return

    regexp = re.compile(regexp)
    for seq in iterable:
        seq.id = regexp.sub(replacement, seq.id)
        yield seq
