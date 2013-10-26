#!/usr/bin/env python

from gfam.interpro import *
from gfam.assignment import *
import sys

# AssignmentReader
def parse(l, seq_id, prots):
    # parses a list of assignments
    if len(l) == 0:
        return
    print "Architecture for %s" % seq_id
    print
    print "InterPro domains found: (", len(l), ")"
    seq = SequenceWithAssignments(l[0].id, l[0].length)
    for ass in l:
        print ass.short_repr(), " ", ass.source, " ", ass.comment
        seq.assign(ass)
    print
    # here remove overlapping domains
    l2 = seq.assignments

    #
    print "InterPro domains found after removin overlap: (", len(l2), ")"
    for ass in l2:
        print ass.short_repr(), "(", len(l2), ")", ass.source
    print
    t = TreeRepresentation(l2)
    print "Representation of the architecture (with positions):"
    print t.get_string_positions()
    print "Representation of the architecture (without positions): "
    print t.get_string()
    print "-----------------------------------------------"
    if prots == 1000:
        sys.exit(-1)

if __name__ == "__main__":
    prev_seq = ""
    seq = None 
    assignments = []
    prots = 0
    for ass in AssignmentReader(sys.argv[1]):
        if ass.id != prev_seq:
            parse(assignments, prev_seq, prots)
            assignments = []
            prots += 1

        assignments.append(ass)
        prev_seq = ass.id

    prots += 1
    parse(assignments, prev_seq, prots)
