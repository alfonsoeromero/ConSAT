#!/usr/bin/env python
"""Application that calculates the domain architecture of each gene and
outputs them in a simple text-based format.
"""

from collections import defaultdict

import operator
import optparse
import sys
import re

from gfam.assignment import AssignmentOverlapChecker, SequenceWithAssignments
from gfam.interpro import InterPro, InterProNames
from gfam.scripts import CommandLineApp
from gfam.utils import complementerset, redirected

__author__  = "Tamas Nepusz"
__email__   = "tamas@cs.rhul.ac.uk"
__copyright__ = "Copyright (c) 2010, Tamas Nepusz"
__license__ = "GPL"

class FindDomainArchitectureApp(CommandLineApp):
    """\
    Usage: %prog [options] interpro_file clustering_file

    Application that calculates the domain architecture of each gene and
    oputs them in a simple text-based format, given the filtered InterPro
    assignments and a clustering of the unknown regions in a separate file.
    """

    short_name = "find_domain_arch"

    def __init__(self, *args, **kwds):
        super(FindDomainArchitectureApp, self).__init__(*args, **kwds)
        self.seqcat = {}

    def create_parser(self):
        """Creates the command line parser for this application"""
        parser = super(FindDomainArchitectureApp, self).create_parser()
        parser.add_option("-s", "--min-size", dest="min_size",
                metavar="N",
                help="consider only clusters with at least N different sequences (not just fragments) as novel domains",
                config_key="analysis:find_domain_arch/min_novel_domain_size",   
                default=2, type=int)
        parser.add_option("-S", "--sequences",
                dest="sequences_file", metavar="FILE",
                help="FASTA file containing all the sequences of the representative gene model",
                config_key="file.input.sequences", default=None)
        parser.add_option("-i", "--interpro-parent-child-file",
                dest="interpro_parent_child_file",
                metavar="FILE",
                help="use the InterPro parent-child FILE to remap IDs",
                config_key="file.mapping.interpro_parent_child",
                default=None)
        parser.add_option("--details",
                dest="details", metavar="FILE",
                help="print more details about the domain architecture into FILE",
                config_key="generated/file.domain_architecture_details",
                default=None)
        parser.add_option("--stats",
                dest="stats", metavar="FILE",
                help="print genome-level statistics about the domain architectures into FILE",
                config_key="generated/file.domain_architecture_stats",
                default=None)
        parser.add_option("--new_domains_table",
                dest="new_domains_table", metavar="FILE",
                help="prints a table with the new domains, one per line, into FILE",
                config_key="generated/file.new_domains_table",
                default=None)
        parser.add_option("-n", "--names",
                dest="interpro_names_file",
                metavar="FILE",
                help="use the given FILE to assign InterPro IDs to names",
                config_key="file.mapping.interpro2name",
                default=None)
        parser.add_option("-r", "--seq-id-regexp", metavar="REGEXP",
                help="remap sequence IDs using REGEXP",
                config_key="sequence_id_regexp",
                dest="sequence_id_regexp")
        parser.add_option("--max-overlap", metavar="SIZE",
                help="sets the maximum overlap size allowed between "
                     "assignments of the same data source. Default: %default",
                config_key="max_overlap",
                dest="max_overlap", type=int, default=20)
        parser.add_option("--previous-table", metavar="FILE", dest="old_table",
                help="reads a previously built table from a file, trying to "
                "keep the same cluster names as much as possible",
                config_key="analysis:find_domain_arch/previous_domain_table",
                default=None)
        return parser

    def print_new_domains_table(self, table):
        """Prins the new domain table"""
        self.log.info("Printing the new domains table")
        
        table_file = open(self.options.new_domains_table, "w")
        with redirected(stdout=table_file):
            for cluster_name in sorted(table.keys()):
                print cluster_name + "\t" + "\t".join(table[cluster_name])
        table_file.close()

    def run_real(self):
        """Runs the applications"""
        if len(self.args) != 2:
            self.error("exactly two input files are expected")

        AssignmentOverlapChecker.max_overlap = self.options.max_overlap

        if self.options.interpro_parent_child_file:
            self.log.info("Loading InterPro parent-child assignments from %s..." % \
                    self.options.interpro_parent_child_file)
            self.interpro = InterPro.FromFile(self.options.interpro_parent_child_file)
        else:
            self.interpro = InterPro()

        self.interpro_names = InterProNames.FromFile(self.options.interpro_names_file)

        if self.options.details:
            self.details_file = open(self.options.details, "w")
        else:
            self.details_file = None

        if self.options.old_table:
            self.process_old_table(self.options.old_table)
            self.using_old_table = True
        else:
            self.using_old_table = False
            self.current_cluster_id = 1

        interpro_file, clustering_file = self.args
        self.process_interpro_file(interpro_file)
        table = self.process_clustering_file(clustering_file)
        self.sort_by_domain_architecture()

        if self.options.new_domains_table:
           self.print_new_domains_table(table)

        for seqs in self.domain_archs.itervalues():
            seqs.sort()

        self.domain_archs = self.domain_archs.items()
        self.domain_archs.sort(key=lambda x: len(x[1]), reverse=True)

        for domain_arch, members in self.domain_archs:
            if domain_arch:
                arch_str = ";".join(domain_arch)
            else:
                arch_str = "NO_ASSIGNMENT"
                arch_str_pos = "NO_ASSIGNMENT"
                arch_desc = "NO_DESCRIPTION"

            family_length = len(members)
            for member in members:
                seq = self.seqcat[member]
                if domain_arch:
                    # TODO: make a better representation, considering overlaps
                    # in this way: IPR1(IPR2) means that IPR2 is inserted into
                    # IPR1
                    arch_str_pos = ";".join(assignment.short_repr() \
                            for assignment in seq.assignments)
                    arch_desc = ";".join( \
                            self.interpro_names[assignment.domain]
                            for assignment in seq.assignments
                    )
                print "%s\t%d\t%s\t%d\t%s\t%s" % (member, seq.length, arch_str, \
                                              family_length, arch_str_pos, \
                                              arch_desc)

        self.details_file.close()

        if self.options.stats:
            stats_file = open(self.options.stats, "w")

            total_residues, covered_residues, covered_residues_nonnovel = 0.0, 0, 0
            nonnovel_sources = complementerset(["Novel"])

            for seq in self.seqcat.itervalues():
                total_residues += seq.length
                covered_residues += round(seq.coverage() * seq.length)
                covered_residues_nonnovel += round(seq.coverage(sources=nonnovel_sources) * seq.length)

            all_archs = set(arch for arch, _ in self.domain_archs)
            num_archs = len(all_archs)
            if "" in self.domain_archs:
                num_archs -= 1

            def exclude_novel_domains(domain_architecture):
                """Excludes novel domains from a domain architecture and returns
                the filtered domain architecture as a tuple."""
                return tuple(a for a in domain_architecture if not a.startswith("NOVEL"))

            archs_without_novel = set(exclude_novel_domains(arch)
                    for arch in all_archs)
            if () in archs_without_novel:
                archs_without_novel.remove(())
            num_archs_without_novel = len(archs_without_novel)

            num_seqs_with_nonempty_domain_arch = \
                    sum(len(value) for key, value in self.domain_archs if key)
            num_seqs_with_nonempty_domain_arch_ignore_novel = \
                    sum(len(value) for key, value in self.domain_archs
                        if exclude_novel_domains(key) in archs_without_novel)
            num_seqs_with_nonempty_nonnovel_domain_arch = \
                    sum(len(value) for key, value in self.domain_archs
                            if key and not any(a.startswith("NOVEL") for a in key))

            with redirected(stdout=stats_file):
                print "Domain architectures"
                print "===================="
                print ""
                print "Non-empty: %d" % num_archs
                print "Non-empty (when ignoring novel domains): %d" % num_archs_without_novel
                print ""
                print "Sequences"
                print "========="
                print ""
                print "Total: %d" % len(self.seqcat)
                print "With at least one domain: %d (%.4f%%)" %\
                        (num_seqs_with_nonempty_domain_arch,
                         100.0 * num_seqs_with_nonempty_domain_arch / len(self.seqcat))
                print "With at least one non-novel domain: %d (%.4f%%)" %\
                        (num_seqs_with_nonempty_domain_arch_ignore_novel,
                         100.0 * num_seqs_with_nonempty_domain_arch_ignore_novel / len(self.seqcat))
                print "With at least one domain and no novel domains: %d (%.4f%%)" %\
                        (num_seqs_with_nonempty_nonnovel_domain_arch,
                         100.0 * num_seqs_with_nonempty_nonnovel_domain_arch / len(self.seqcat))
                print ""
                print "Residues"
                print "========"
                print ""
                print "Total: %d" % total_residues
                print "Covered: %d (%.4f%%)" % (covered_residues, 100.0*covered_residues/total_residues)
                print "Covered by non-novel: %d (%.4f%%)" % (covered_residues_nonnovel, 100.0*covered_residues_nonnovel/total_residues)
            stats_file.close()

    def process_old_table(self, old_table_file):
        # we build 3 structures: a dict mapping
        # each fragment to its cluster id, and
        # a dict mapping each sequence id to the
        # list of proteins it belongs to.
        # For each cluster we store the set of sequences too.
        self.fragments_per_cluster = defaultdict(list)
        self.cluster_per_fragment = dict()
        self.fragments_per_seq = defaultdict(list)
        self.current_cluster_id = 1

        for line in open(old_table_file,'r'):
            fields = line.split()
            cluster_name, fragments = fields[0], fields[1:]
            cluster_num = int(re.findall('\d+', cluster_name)[0])
            self.current_cluster_id = max(self.current_cluster_id, cluster_num)

            self.fragments_per_cluster[cluster_name] = fragments
            for fragment in fragments:
                sequence = fragment.split(':')[0]
                self.cluster_per_fragment[fragment] = cluster_name
                self.fragments_per_seq[sequence].append(fragment)
                
    def process_interpro_file(self, interpro_file):
        from gfam.scripts.find_unassigned import FindUnassignedApp
        unassigned_app = FindUnassignedApp()
        unassigned_app.set_sequence_id_regexp(self.options.sequence_id_regexp)
        unassigned_app.process_sequences_file(self.options.sequences_file)
        unassigned_app.process_infile(interpro_file)
        self.seqcat = unassigned_app.seqcat
        for seq_id in set(unassigned_app.seq_ids_to_length.keys()) - set(self.seqcat.keys()):
            self.seqcat[seq_id] = SequenceWithAssignments(seq_id, \
                                  unassigned_app.seq_ids_to_length[seq_id])

    def find_domain_id(self, fragments):
        """Maps a set of fragments to the most likely
        cluster from the old_cluster_trable, if this is possible.
        The identifier (number) of this cluster is returned, if
        anyone is found, or -1 if no matching cluster is found.
        """
        # 1.- we vote each possible cluster
        num_fragments = len(fragments)
        num_unassigned = 0
        votes_per_cluster = defaultdict(int)
        for fragment in fragments:
            if fragment in self.cluster_per_fragment:
                # 1.1- if the fragment has a previously assigned cluster,
                # we vote on them
                votes_per_cluster[self.cluster_per_fragment[fragment]] += 1
            else:
                # 1.2.- if not, we check for fragments near this one in the
                # same sequence
                fields = fragment.split(':')
                sequence = fields[0]
                if sequence in fragments_per_seq:
                    # If there are other fragments in the same sequence...
                    matching_seq = defaultdict(float)
                    seq_from, seq_to = map(int, fields[1].split('-'))
                    for other_fragment in fragments_per_seq[sequence]:
                        oseq_from, oseq_to = map(int, other_fragment.split(':')[1].split('-'))
                        if oseq_from <= seq_to:
                            # there is overlap
                            matching_seq[other_fragment] = float(seq_to - oseq_from + 1) / float(oseq_to - seq_from + 1)
                        elif seq_from <= oseq_to:
                            matching_seq[other_fragment] = float(oseq_to - seq_from + 1) / float(seq_to - oseq_from + 1)
                    if len(matching_seq) > 0:
                        matched = max(matching_seq.iteritems(), key=lambda x: x[1])[0]
                        votes_per_cluster[self.cluster_per_fragment[matched]] += 1
                    else:
                        num_assigned += 1
                else:
                    # the sequence has not other fragments assigned to clusters
                    num_unassigned += 1

        # 2.- we count the votes and take a decision
        if len(votes_per_cluster) == 0:
            return -1
        elif len(votes_per_cluster) == 1:
            return votes_per_cluster.items()[0][0]
        else:
            return max(matching_seq.iteritems(), key=lambda x: x[1])[0]

    def process_clustering_file(self, fname):
        table = defaultdict()
        f = open(fname)
        for line in f:
            fragments = line.strip().split()
            if len(set([x.split(":",1)[0] for x in fragments])) < self.options.min_size:
                continue
            if not self.using_old_table:
                domain_name = "NOVEL%05d" % self.current_cluster_id
                self.current_cluster_id +=1
            else:
                # We find the name of the best matching cluster given this set of
                # sequences
                domain_id = self.find_domain_id(fragments)
                if domain_id != -1:
                    domain_name = "NOVEL%05d" % domain_id
                else:
                    domain_name = "NOVEL%05d" % self.current_cluster_id
                    self.current_cluster_id += 1
            sequences = []
            for id in ids:
                sequences.append(id)
                seq_id, _, limits = id.rpartition(":")
                start, end = map(int, limits.split("-"))
                self.seqcat[seq_id].assign_(start, end, domain_name)
            table[domain_name] = sequences
        f.close()
        return table

    def sort_by_domain_architecture(self):
        self.domain_archs = defaultdict(list)
        for seq_id, seq in self.seqcat.iteritems():
            assignments = sorted(seq.assignments, key=operator.attrgetter("start"))
            domains = []
            if self.details_file:
                print >>self.details_file, seq_id

            primary_source = set()

            new_assignments = []
            for assignment in assignments:
                new_assignment = assignment.resolve_interpro_ids(self.interpro)
                if assignment.comment == "1":
                    primary_source.add(assignment.source)
                domains.append(new_assignment.domain)
                new_assignments.append(new_assignment)
            self.domain_archs[tuple(domains)].append(seq_id)

            if not primary_source:
                primary_source = None
            else:
                primary_source = ", ".join(primary_source)

            if self.details_file:
                seq2 = SequenceWithAssignments(seq.name, seq.length)
                seq2.assignments = [assignment for assignment in assignments \
                                    if assignment.source != "Novel"]
                sources = sorted(set(assignment.source \
                        for assignment in assignments \
                        if assignment.source != "Novel"))

                print >>self.details_file, "    Primary assignment source:", primary_source
                print >>self.details_file, "    Number of data sources used:", len(sources)
                print >>self.details_file, "    Data sources: %s" % ", ".join(sources)
                print >>self.details_file, "    Coverage: %.3f" % seq.coverage()
                print >>self.details_file, "    Coverage w/o novel domains: %.3f" % seq2.coverage()
                for assignment in assignments:
                    attrs = assignment._asdict()
                    if assignment.comment is None and \
                       assignment.domain.startswith("NOVEL"):
                        attrs["comment"] = "novel"
                    row = "    %(start)4d-%(end)4d: %(domain)s "\
                          "(%(source)s, stage: %(comment)s)" % attrs
                    print >>self.details_file, row,
                    interpro_id = assignment.interpro_id
                    if not interpro_id and assignment.domain in self.interpro.mapping:
                        interpro_id = self.interpro.mapping[assignment.domain]
                    if interpro_id:
                        anc = self.interpro.tree.get_most_remote_ancestor(interpro_id)
                        if interpro_id == anc:
                            print >>self.details_file, "(InterPro ID: %s)" % anc
                        else:
                            print >>self.details_file, "(InterPro ID: %s --> %s)" % (interpro_id, anc)
                        if anc in self.interpro_names:
                            print >>self.details_file, " "*(row.index(":")+1), self.interpro_names[anc]
                    else:
                        print >>self.details_file, ""
                        if assignment.domain in self.interpro_names:
                            print >>self.details_file, " "*(row.index(":")+1), self.interpro_names[assignment.domain]
                print >>self.details_file, ""

            seq.assignments = new_assignments

if __name__ == "__main__":
    sys.exit(FindDomainArchitectureApp().run())

