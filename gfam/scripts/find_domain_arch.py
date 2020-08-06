#!/usr/bin/env python
"""Application that calculates the domain architecture of each gene and
outputs them in a simple text-based format.
"""
import operator
import re
import sys
from collections import defaultdict

from gfam.assignment import (AssignmentOverlapChecker, SequenceWithAssignments,
                             TreeRepresentation)
from gfam.interpro import InterProNames
from gfam.tasks.assignment_source_filter.interpro_file_factory import \
    InterproFileFactory
from gfam.tasks.find_domain_arch.base_find_domain_arch import \
    BaseFindDomainArch
from gfam.tasks.find_domain_arch.stats.architecture_stats import \
    ArchitectureStats
from gfam.tasks.find_domain_arch.stats.stats_file_printer import \
    StatsFilePrinter
from gfam.utils import complementerset, redirected

__author__ = "Tamas Nepusz"
__email__ = "tamas@cs.rhul.ac.uk"
__copyright__ = "Copyright (c) 2010, Tamas Nepusz"
__license__ = "GPL"


class FindDomainArchitectureApp(BaseFindDomainArch):
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
        self.interpro = None
        self.details_file = None
        self.using_old_table = False
        self.interpro_names = None
        self.current_cluster_id = -1
        self.domain_archs = None
        self.cluster_per_fragment = None
        self.fragments_per_cluster = None
        self.fragments_per_seq = None

    def create_parser(self):
        """Creates the command line parser for this application"""
        parser = super(FindDomainArchitectureApp, self).create_parser()
        self._add_common_options(parser)

        parser.add_option("--new_domains_table",
                          dest="new_domains_table", metavar="FILE",
                          help="prints a table with the new domains,"
                               " one per line, into FILE",
                          config_key="generated/file.new_domains_table",
                          default=None)
        parser.add_option("--previous-table", metavar="FILE", dest="old_table",
                          help="reads a previously built table from a file,"
                               " trying to "
                          "keep the same cluster names as much as possible",
                          config_key="analysis:find_domain_arch/"
                                     "previous_domain_table",
                          default=None)
        parser.add_option("-p", metavar="STRING", dest="prefix",
                          help="prefix for the new discovered domains ('NOVEL'"
                               " by default)",
                          config_key="analysis:find_domain_arch/prefix",
                          default='NOVEL')
        return parser

    def print_new_domains_table(self, table):
        """Prints the new domain table"""
        self.log.info("Printing the new domains table")
        table_file = open(self.options.new_domains_table, "w")
        with redirected(stdout=table_file):
            for cluster_name in sorted(table.keys()):
                print(cluster_name + "\t" + "\t".join(table[cluster_name]))
        table_file.close()

    def _set_app_parameters(self):
        AssignmentOverlapChecker.max_overlap = self.options.max_overlap
        AssignmentOverlapChecker.log = self.log

        if self.options.prefix:
            self.prefix = self.options.prefix
        else:
            self.prefix = "NOVEL"

        self.interpro = InterproFileFactory.get_from_file(
            self.options.interpro_file)
        self.interpro_names = InterProNames.from_file(
            self.options.interpro_names_file)

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

    def run_real(self):
        """Runs the applications"""
        if len(self.args) != 2:
            self.error("exactly two input files are expected")
        self._set_app_parameters()

        interpro_file, clustering_file = self.args
        self.process_interpro_file(interpro_file)
        table = self.process_clustering_file(clustering_file)
        self.sort_by_domain_architecture()

        if self.options.new_domains_table:
            self.print_new_domains_table(table)

        for seqs in self.domain_archs.values():
            seqs.sort()

        self.domain_archs = self.domain_archs.items()
        self.domain_archs.sort(key=lambda x: len(x[1]), reverse=True)

        for domain_arch, members in self.domain_archs:
            if domain_arch:
                arch_str = domain_arch
            else:
                arch_str = "NO_ASSIGNMENT"
                arch_str_pos = "NO_ASSIGNMENT"
                arch_desc = "NO_DESCRIPTION"

            family_length = len(members)
            for member in members:
                seq = self.seqcat[member]
                if domain_arch:
                    arch_str_pos = seq.architecture_pos
                    arch_desc = ";".join(
                        self.interpro_names[assignment.domain]
                        for assignment in seq.assignments)
                print("%s\t%d\t%d\t%s\t%d\t%s\t%s" % (member, seq.length,
                                                      seq.num_covered(),
                                                      arch_str,
                                                      family_length,
                                                      arch_str_pos,
                                                      arch_desc))
        self.details_file.close()

        if self.options.stats:
            architecture_stats = self._compute_architecture_stats()
            stats_printer = StatsFilePrinter()
            stats_printer.print_stats_file(self.options.stats,
                                           architecture_stats)

    def _compute_architecture_stats(self) -> ArchitectureStats:
        total_residues = 0.0
        covered_residues = 0
        covered_residues_nonnovel = 0
        nonnovel_sources = complementerset(["Novel"])

        for seq in self.seqcat.values():
            total_residues += seq.length
            covered_residues += round(seq.coverage() * seq.length)
            covered_residues_nonnovel += round(
                seq.coverage(sources=nonnovel_sources) * seq.length)

        all_archs = set(arch for arch, _ in self.domain_archs)
        num_archs = len(all_archs)
        if "" in self.domain_archs:
            num_archs -= 1

        def split_arch(arch):
            return [x for x in arch.replace("{", ";")
                                   .replace("}", ";")
                                   .split(";") if x]

        def exclude_novel_domains(domain_architecture):
            """Excludes novel domains from a domain architecture and returns
            the filtered domain architecture as a tuple."""
            return tuple(a for a in split_arch(domain_architecture)
                         if not a.startswith(self.prefix))

        archs_without_novel = set(exclude_novel_domains(arch)
                                  for arch in all_archs)
        if () in archs_without_novel:
            archs_without_novel.remove(())
        num_archs_without_novel = len(archs_without_novel)

        num_seqs_with_nonempty_domain_arch = sum(len(value) for ke, value
                                                 in self.domain_archs if ke
                                                 and ke != "NO_ASSIGNMENT")
        num_seqs_with_nonempty_domain_arch_ignore_novel = \
            sum(len(value) for key, value in self.domain_archs
                if exclude_novel_domains(key) in archs_without_novel
                and key != "NO_ASSIGNMENT")
        num_seqs_with_nonempty_nonnovel_domain_arch = \
            sum(len(value) for key, value in self.domain_archs
                if key and not any(a.startswith(self.prefix) for a in
                                   split_arch(key)) and
                key != "NO_ASSIGNMENT")
        return ArchitectureStats(
            num_archs,
            len(self.seqcat),
            total_residues,
            num_archs_without_novel,
            num_seqs_with_nonempty_domain_arch,
            num_seqs_with_nonempty_domain_arch_ignore_novel,
            num_seqs_with_nonempty_nonnovel_domain_arch,
            covered_residues,
            covered_residues_nonnovel)

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

        for line in open(old_table_file, 'r'):
            fields = line.split()
            cluster_name, fragments = fields[0], fields[1:]
            cluster_num = int(re.findall(r'\d+', cluster_name)[0])
            self.current_cluster_id = max(self.current_cluster_id, cluster_num)

            self.fragments_per_cluster[cluster_name] = fragments
            for fragment in fragments:
                sequence = fragment.split(':')[0]
                self.cluster_per_fragment[fragment] = cluster_name
                self.fragments_per_seq[sequence].append(fragment)

    def find_domain_id(self, fragments):
        """Maps a set of fragments to the most likely
        cluster from the old_cluster_trable, if this is possible.
        The identifier (number) of this cluster is returned, if
        anyone is found, or -1 if no matching cluster is found.
        """
        # 1.- we vote each possible cluster
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
                if sequence in self.fragments_per_seq:
                    # If there are other fragments in the same sequence...
                    matching_seq = defaultdict(float)
                    seq_from, seq_to = map(int, fields[1].split('-'))
                    for other_fragment in self.fragments_per_seq[sequence]:
                        oseq_from, oseq_to = map(
                            int, other_fragment.split(':')[1].split('-'))
                        if oseq_from <= seq_to:
                            # there is overlap
                            matching_seq[other_fragment] = \
                                float(seq_to - oseq_from + 1) / float(
                                    oseq_to - seq_from + 1)
                        elif seq_from <= oseq_to:
                            matching_seq[other_fragment] = float(
                                oseq_to - seq_from + 1) /\
                                float(seq_to - oseq_from + 1)
                    if matching_seq:
                        matched = max(matching_seq.items(),
                                      key=lambda x: x[1])[0]
                        clu = self.cluster_per_fragment[matched]
                        votes_per_cluster[clu] += 1

        # 2.- we count the votes and take a decision
        if not votes_per_cluster:
            return -1
        elif len(votes_per_cluster) == 1:
            return votes_per_cluster.items()[0][0]
        return max(matching_seq.items(), key=lambda x: x[1])[0]

    def process_clustering_file(self, fname: str):
        table = defaultdict()
        cluster_file = open(fname)
        for line in cluster_file:
            fragments = line.strip().split()
            if len(set([x.split(":", 1)[0] for x in fragments]))\
               < self.options.min_size:
                continue
            if not self.using_old_table:
                domain_name = self.prefix + "%05d" % self.current_cluster_id
                self.current_cluster_id += 1
            else:
                # We find the name of the best matching cluster given
                # this set of sequences
                domain_id = self.find_domain_id(fragments)
                if domain_id != -1:
                    domain_name = self.prefix + f"{domain_id:05}"
                else:
                    domain_name = self.prefix + f"{self.current_cluster_id:05}"
                    self.current_cluster_id += 1
            sequences = []
            for frag_id in fragments:  # something wrong here...
                sequences.append(frag_id)
                seq_id, _, limits = frag_id.rpartition(":")
                start, end = map(int, limits.split("-"))
                self.seqcat[seq_id].assign_(start, end, domain_name)
            table[domain_name] = sequences
        cluster_file.close()
        return table

    def sort_by_domain_architecture(self):
        self.domain_archs = defaultdict(list)

        for seq_id, seq in self.seqcat.items():
            assignments = sorted(seq.assignments,
                                 key=operator.attrgetter("start"))
            domains = []
            if self.details_file:
                print(seq_id, file=self.details_file)

            primary_source = set()

            new_assignments = []
            for assignment in assignments:
                new_assignment = assignment.resolve_interpro_ids(self.interpro)
                if assignment.comment == "1":
                    primary_source.add(assignment.source)
                domains.append(new_assignment.domain)
                new_assignments.append(new_assignment)
            tree_arch = TreeRepresentation(new_assignments, self.interpro)
            seq.architecture = tree_arch.get_string()
            seq.architecture_pos = tree_arch.get_string_positions()

            self.domain_archs[seq.architecture].append(seq_id)

            if not primary_source:
                primary_source = None
            else:
                primary_source = ", ".join(primary_source)

            if self.details_file:
                seq2 = SequenceWithAssignments(seq.name, seq.length)
                seq2.assignments = [assignment for assignment in assignments
                                    if assignment.source != "Novel"]
                sources = sorted(set(assignment.source
                                     for assignment in assignments
                                     if assignment.source != "Novel"))

                print("    Primary assignment source: {}".format(
                    primary_source), file=self.details_file)
                print("    Number of data sources used: {}".format(
                    len(sources)), file=self.details_file)
                print("    Data sources: %s" % ", ".join(sources),
                      file=self.details_file)
                print("    Coverage: %.3f" % seq.coverage(),
                      file=self.details_file)
                print("    Coverage w/o novel domains: %.3f"
                      % seq2.coverage(), file=self.details_file)
                for assignment in assignments:
                    attrs = assignment._asdict()
                    if assignment.comment is None and \
                       assignment.domain.startswith(self.prefix):
                        attrs["comment"] = "novel"
                    row = "    %(start)4d-%(end)4d: %(domain)s "\
                          "(%(source)s, stage: %(comment)s)" % attrs
                    print(row, file=self.details_file)
                    interpro_id = assignment.interpro_id
                    if not interpro_id\
                       and assignment.domain in self.interpro.mapping:
                        interpro_id = self.interpro.mapping[assignment.domain]
                    if interpro_id:
                        anc = self.interpro.tree.get_most_remote_ancestor(
                            interpro_id)
                        if interpro_id == anc:
                            print("(InterPro ID: %s)" % anc,
                                  file=self.details_file)
                        else:
                            print("(InterPro ID: %s --> %s)"
                                  % (interpro_id, anc), file=self.details_file)
                        if anc in self.interpro_names:
                            print("{}{}".format(" " * (row.index(":") + 1),
                                                self.interpro_names[anc]),
                                  file=self.details_file)
                    else:
                        print("", file=self.details_file)
                        if assignment.domain in self.interpro_names:
                            print("{}{}".format(
                                " " * (row.index(":") + 1),
                                self.interpro_names[assignment.domain]),
                                file=self.details_file)
                print("", file=self.details_file)

            seq.assignments = new_assignments


if __name__ == "__main__":
    sys.exit(FindDomainArchitectureApp().run())
