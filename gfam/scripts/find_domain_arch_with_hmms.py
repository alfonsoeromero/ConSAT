#!/usr/bin/env python
"""Application that calculates the domain architecture of each gene and
outputs them in a simple text-based format.
"""
import operator
import sys
from collections import defaultdict

from gfam.assignment_utils.assignment_overlap_checker import \
    AssignmentOverlapChecker
from gfam.assignment_utils.sequence_with_assignments import \
    SequenceWithAssignments
from gfam.assignment_utils.tree_representation import TreeRepresentation
from gfam.interpro import InterProNames
from gfam.tasks.assignment_source_filter.interpro_file_factory import \
    InterproFileFactory
from gfam.tasks.find_domain_arch.base_find_domain_arch import \
    BaseFindDomainArch
from gfam.tasks.find_domain_arch.stats.architecture_stats import \
    ArchitectureStats
from gfam.tasks.find_domain_arch.stats.stats_file_printer import \
    StatsFilePrinter
from gfam.utilities.complementer_set import ComplementerSet

__author__ = "Alfonso E. Romero, Tamas Nepusz"
__email__ = "aeromero@cs.rhul.ac.uk"
__copyright__ = "Copyright (c) 2012, Alfonso E. Romero, Tamas Nepusz"
__license__ = "GPL"

__all__ = ["FindDomainArchitectureWithHMMsApp"]


class FindDomainArchitectureWithHMMsApp(BaseFindDomainArch):
    """\
    Usage: %prog [options] interpro_file hmms_file

    Application that calculates the domain architecture of each gene and
    oputs them in a simple text-based format, given the filtered InterPro
    assignments and a file assignment.
    """

    short_name = "find_domain_architecture_with_hmms"

    def __init__(self, *args, **kwds):
        super(FindDomainArchitectureWithHMMsApp, self).__init__(*args, **kwds)
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
        self.hmm_domains = None

    def create_parser(self):
        """Creates the command line parser for this application"""
        parser = super(FindDomainArchitectureWithHMMsApp, self).create_parser()
        self._add_common_options(parser)
        return parser

    def run_real(self):
        """Runs the applications"""
        if len(self.args) != 2:
            self.error("exactly two input files are expected")

        AssignmentOverlapChecker.max_overlap = self.options.max_overlap

        self.interpro = InterproFileFactory.get_from_file(
            self.options.interpro_file)
        self.interpro_names = InterProNames.from_file(
            self.options.interpro_names_file)

        if self.options.details:
            self.details_file = open(self.options.details, "w")
        else:
            self.details_file = None

        interpro_file, hmmer_file = self.args
        self.process_interpro_file(interpro_file, True)
        self.process_hmmer_file(hmmer_file)
        self.sort_by_domain_architecture()

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
        covered_residues, covered_residues_nonnovel = 0, 0
        nonnovel_sources = ComplementerSet(["Novel"])

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
            return [x for x in arch.replace("{", ";").replace("}", ";")
                                                     .split(";") if x]

        def exclude_novel_domains(domain_architecture):
            """Excludes novel domains from a domain architecture and returns
            the filtered domain architecture as a tuple."""
            return tuple(a for a in split_arch(domain_architecture)
                         if a not in self.hmm_domains)

        archs_without_novel = set(exclude_novel_domains(arch)
                                  for arch in all_archs)
        if () in archs_without_novel:
            archs_without_novel.remove(())
        num_archs_without_novel = len(archs_without_novel)
        num_seqs_with_nonempty_domain_arch = sum(
            len(value) for key, value in self.domain_archs if key
            and key != "NO_ASSIGNMENT")
        num_seqs_with_nonempty_domain_arch_ignore_novel = sum(
            len(value) for key, value in self.domain_archs
            if exclude_novel_domains(key) in archs_without_novel
            and key != "NO_ASSIGNMENT")
        num_seqs_with_nonempty_nonnovel_domain_arch = sum(
            len(value) for ke, value in self.domain_archs
            if ke and not any(a in self.hmm_domains for a in ke)
            and ke != "NO_ASSIGNMENT")
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

    def process_hmmer_file(self, fname):
        self.hmm_domains = set()
        for line in open(fname):
            id_prot, model, _ = line.strip().split()
            # evalue is not used, at the moment
            seq_id, _, limits = id_prot.rpartition(":")
            start, end = map(int, limits.split("-"))
            self.seqcat[seq_id].assign_(start, end, model)
            self.hmm_domains.add(model)

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
                print("    Coverage w/o novel domains: %.3f" % seq2.coverage(),
                      file=self.details_file)
                for assignment in assignments:
                    attrs = assignment._asdict()
                    if assignment.comment is None and \
                       assignment.domain in self.hmm_domains:
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
                                  % (interpro_id, anc),
                                  file=self.details_file)
                        if anc in self.interpro_names:
                            print("{}{}".format(" "*(row.index(":")+1),
                                                self.interpro_names[anc]),
                                  file=self.details_file)
                    else:
                        print("", file=self.details_file)
                        if assignment.domain in self.interpro_names:
                            print("{}{}".format(
                                " "*(row.index(":")+1),
                                self.interpro_names[assignment.domain]),
                                file=self.details_file)
                print("", file=self.details_file)

            seq.assignments = new_assignments


if __name__ == "__main__":
    sys.exit(FindDomainArchitectureWithHMMsApp().run())
