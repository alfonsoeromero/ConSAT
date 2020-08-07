#!/usr/bin/env python
import os.path
import sys
from collections import defaultdict

from gfam.go import Tree as GOTree
from gfam.utilities.open_anything import open_anything

from result_file import ResultFileReader

__author__ = "Alfonso E. Romero"
__email__ = "aeromero@cs.rhul.ac.uk"
__copyright__ = "Copyright (c) 2013, Alfonso E. Romero"
__license__ = "GPL"


class Evaluation(object):
    """\
    Usage: %prog [options]

    Application that evaluates protein assignments to GO terms made by
    GFam.

    It is supposed to receive at least, a gene_ontology file, a GOA file
    with the ground truth (only experimental codes will be considered
    as truth), and an assignment of proteins to GO terms (with the
    same format as the overrepresentation file). The output measures
    will be written to an output directory.

    If two GO assignment files are included in the call, some comparison
    measures between the two will also be carried out
    """

    valid_ev_codes = ["EXP", "IDA", "IPI", "IMP", "IGI", "IEP", "TAS", "IC"]
    significance_levels = [1.0, 0.1, 0.05, 0.01, 0.001]
    cache_file = dict()

    def __init__(self, go_file, goa_file, output_directory, proteins):
        self.go_tree = GOTree.from_obo(go_file)
        self.terms_per_ontology = defaultdict(set)
        self.ontology_per_term = {}
        for _goterm, term in self.go_tree.terms.items():
            try:
                namespace = str(term.tags["namespace"][0])
                goterm_id = str(term.tags["id"][0])
                self.terms_per_ontology[namespace].add(goterm_id)
                self.ontology_per_term[goterm_id] = namespace
            except KeyError:
                # this is done to avoid reading the Typedef entries
                pass

        for ontology, terms in self.terms_per_ontology.items():
            print("{} {}".format(ontology, len(terms)))

        self.goa = self._read_goa_file(goa_file)
        self.output_dir = output_directory
        self._create_output_directory()
        self.proteins = set()
        with open(proteins, "r") as fin:
            for line in fin:
                self.proteins.add(line.strip())

    def _read_result_file(self, file_name, name):
        """Reads the result file into a dictionary structure where,
        for each protein identifier we get a list of tuples (goterm, pvalue)
        representing the assignment itself
        """
        if name not in self.cache_file:
            reader = ResultFileReader(file_name)
            self.cache_file[name] = reader.get_result_as_dict()
        return self.cache_file[name]

    def _get_all_proteins(self, assoc, name):
        """Returns all protein identifiers in the result file
        """
        identifier_list = self._read_result_file(assoc, name)
        return [prot for prot in identifier_list]

    def _get_orphan_proteins(self, assoc, name, significance=0.05):
        """Returns the set of protein ids without any
        annotation
        """
        li = self._read_result_file(assoc, name)
        return [prot for prot in li if
                len([p for p in li[prot] if p[1] <= significance]) == 0]

    def _evaluate_orphans(self, overrep_file, transfer_file):
        orphan_file = os.path.join(self.output_dir, "01_orphans")
        sys.stdout = open(orphan_file, "w")
        all_proteins = self.proteins
        num_all = float(len(all_proteins))
        prots_overrep = self._get_all_proteins(overrep_file, "overrep")
        num_prots_overrep = float(len(prots_overrep))
        prots_transfer = self._get_all_proteins(transfer_file, "transfer")
        num_prots_transfer = float(len(prots_transfer))
        print("Analysis of orphan proteins")
        print("=============================")
        print("Total number of proteins with some result: {}".format(
            len(all_proteins)))
        for significance in self.significance_levels:
            print("When selecting only assignments with p <= {}".format(
                significance))
            orphans_overrep = self._get_orphan_proteins(overrep_file,
                                                        "overrep",
                                                        significance)
            orphans_transfer = self._get_orphan_proteins(transfer_file,
                                                         "transfer",
                                                         significance)
            num_orphans_overrep = float(len(orphans_overrep))
            num_predicted_overrep = num_prots_overrep - num_orphans_overrep
            num_orphans_transfer = float(len(orphans_transfer))
            num_predicted_transfer = num_prots_transfer - num_orphans_transfer
            print("Number of orphans in overrep: %d" % (num_orphans_overrep))
            print("Number of proteins with some prediction in overrep: %d "
                  % (num_predicted_overrep))
            print("perc of orphans over total proteins | maximum_predicted:"
                  " %2.2f  | %2.2f " % (100.0*num_orphans_overrep / num_all,
                                        100.0*num_orphans_overrep /
                                        num_prots_overrep))
            print("Number of orphans in transfer: %d" % num_orphans_transfer)
            print("Number of proteins with some prediction in transfer: %d" %
                  (num_predicted_transfer))
            print("perc of orphans over total proteins | "
                  "maximum_predicted: %2.2f  | %2.2f " % (
                      100.0*num_orphans_transfer/num_all,
                      100.0*num_orphans_transfer/num_prots_transfer))
            print("Of which %d are common" % (len(set(orphans_overrep) &
                                                  set(orphans_transfer))))
        sys.stdout = sys.__stdout__

    def _evaluate_single(self, file_name, name):
        """We compare the obtained labels with those in the
        ground truth, for those GO terms obtained and for
        each one of the ontologies
        """
        single_file = os.path.join(self.output_dir,
                                   "02_evaluate_single_" + file_name)
        sys.stdout = open(single_file, "w")
        result = self._read_result_file(file_name, name)
        print("Evaluation of single file {}".format(file_name))
        print("(note we are ignoring the p-values)")
        print()

        # we ignore the thresholds (p-values)
        TP, FP, FN = 0, 0, 0
        Nempty, Nres = 0, 0
        avprec, avrec, avf1 = 0.0, 0.0, 0.0
        avprec_ont = defaultdict(int)
        avrec_ont, avf1_ont = defaultdict(int), defaultdict(int)
        TP_ont, FP_ont = defaultdict(int), defaultdict(int)
        FN_ont = defaultdict(int)
        for protein, terms_gs_ in self.goa.items():
            terms_gs = set(terms_gs_)
            # computation of global results
            if protein in result:
                terms_pred = set([res[0] for res in result[protein]])
                TPi = len(terms_pred & terms_gs)
                FPi = len(terms_pred - terms_gs)
                FNi = len(terms_gs - terms_pred)
                TP += TPi
                FP += FPi
                FN += FNi
                if TPi + FPi > 0:
                    prec = float(TPi) / float(TPi + FPi)
                else:
                    prec = 1.0
                rec = float(TPi) / float(TPi + FNi)
                if prec + rec > 0.0:
                    f1 = 2.0*prec*rec/(prec+rec)
                else:
                    f1 = 0.0
                avprec += prec
                avrec += rec
                avf1 += f1
                for onto, terms in self.terms_per_ontology.items():
                    # computation of results per ontology
                    terms_pred_onto = terms & terms_pred
                    terms_gs_onto = terms & terms_gs
                    TPi_onto = len(terms_pred_onto & terms_gs_onto)
                    FPi_onto = len(terms_pred_onto - terms_gs_onto)
                    FNi_onto = len(terms_gs_onto - terms_pred_onto)
                    if TPi_onto + FPi_onto > 0:
                        prec_onto = float(TPi_onto)/float(TPi_onto + FPi_onto)
                    else:
                        prec_onto = 1.0

                    if TPi_onto + FNi_onto:
                        rec_onto = float(TPi_onto)/float(TPi_onto + FNi_onto)
                    else:
                        rec_onto = 0.0
                    if prec_onto + rec_onto > 0.0:
                        f1_onto = 2.0*prec_onto*rec_onto/float(prec_onto +
                                                               rec_onto)
                    else:
                        f1_onto = 0.0
                    avprec_ont[onto] += prec_onto
                    avrec_ont[onto] += rec_onto
                    avf1_ont[onto] += f1_onto
                    TP_ont[onto] += TPi_onto
                    FP_ont[onto] += FPi_onto
                    FN_ont[onto] += FNi_onto
                Nres += 1
            else:
                Nempty += 1
                FN += len(terms_gs)
                for onto, terms in self.terms_per_ontology.items():
                    FN_ont[onto] += len(terms_gs & terms)
        print("Global results: ")
        print("Precision: {}".format(float(TP)/float(TP+FP)))
        print("Recall: {}".format(float(TP)/float(TP+FN)))
        print("F1: {}".format(2.0*TP/float(2.0*TP + FP + FN)))
        for ontology in self.terms_per_ontology:
            print("Results for ontology {}".format(ontology))
            prec_on = float(TP_ont[ontology])/float(
                TP_ont[ontology] + FP_ont[ontology])
            print("Precision: {}".format(prec_on))
            rec_on = float(TP_ont[ontology])/float(
                TP_ont[ontology] + FN_ont[ontology])
            print("Recall: {}".format(rec_on))
            f1_on = 2.0*TP_ont[ontology]/float(
                2.0*TP_ont[ontology] + FP_ont[ontology] + FN_ont[ontology])
            print("F1: {}".format(f1_on))
        print()
        print("Results per protein: ")
        print("Average precision: {}".format(avprec/float(Nempty + Nres)))
        print("Average precision only on predicted: {}".format(
            avprec/float(Nres)))
        print("Average recall: {}".format(avrec/float(Nempty + Nres)))
        print("Average recall only on predicted: {}".format(avrec/float(Nres)))
        print("Average F1: {}".format(avf1/float(Nempty + Nres)))
        print("Average F1 only on predicted: ", avf1/float(Nres))
        for ontology in self.terms_per_ontology:
            print("Results for ontology {}".format(ontology))
            print("Average precision: {}".format(
                avprec_ont[ontology]/float(Nempty + Nres)))
            print("Average precision only on predicted: {}".format(
                avprec_ont[ontology]/float(Nres)))
            print("Average recall: {}".format(
                avrec_ont[ontology]/float(Nempty + Nres)))
            print("Average recall only on predicted: {}".format(
                avrec_ont[ontology]/float(Nres)))
            print("Average F1: {}".format(
                avf1_ont[ontology]/float(Nempty + Nres)))
            print("Average F1 only on predicted: {}".format(
                avf1_ont[ontology]/float(Nres)))
        sys.stdout = sys.__stdout__

    def _join(self, name1, name2, final_name):
        out = dict()
        for protein, values in self.cache_file[name1].items():
            out[protein] = values

        for protein, values in self.cache_file[name2].items():
            if protein not in out:
                out[protein] = values
            else:
                out[protein].extend(values)
        self.cache_file[final_name] = out
        return out

    def evaluate_two_files(self, overrep_file, transfer_file):
        # first, evaluation of orphan proteins
        self._evaluate_orphans(overrep_file, transfer_file)
        self._evaluate_single(overrep_file, "overrep")
        self._evaluate_single(transfer_file, "transfer")
        self._join("overrep", "transfer", "both")
        self._evaluate_single("both", "both")

    def evaluate_single_file(self, assoc, name="association"):
        self._get_orphan_proteins(assoc, name)

    def _create_output_directory(self):
        """Creates the output directory checking before that it
        does not exist.
        """
        if os.path.exists(self.output_dir):
            print("ERROR: the specified output directory alread exists")
            sys.exit(-1)

        os.makedirs(self.output_dir)

    def _read_goa_file(self, goa_file):
        """Reads the GOA file, return a dict that,
        for each protein id it has a set of strings with the
        associated GO terms after up-propagation
        """
        d = {}
        self.go_terms_gs = set()
        for line in open_anything(goa_file):
            if not line.startswith("!") or line.startswith("#"):
                # split line, obtain protein_id, go_term, and ev_code
                fields = line.split("\t")
                prot_id, goterm, evcode = fields[1], fields[4], fields[6]
                if evcode in self.valid_ev_codes:
                    terms = map(str, self.go_tree.ancestors(goterm))
                    d[prot_id] = terms
                    self.go_terms_gs.update(terms)

        self.go_terms_gs_ont = {}
        for ontology, terms in self.terms_per_ontology.items():
            gs_terms_o = self.go_terms_gs & terms
            self.go_terms_gs_ont[ontology] = gs_terms_o
            print("# of Terms in {} in GOA file: {}".format(
                ontology, len(gs_terms_o)))

        return d


def usage():
    print("""ERROR. The list of arguments is the following:
             gene_ontology_file goa_file output_dir list_of_proteins_file
             association_file_1 [association_file_2]

             If two association files are given, the first is
             assumed to come from overrrepresentation analysis
             of InterPro labels, and the other one from transfer""")


if __name__ == "__main__":
    if len(sys.argv) < 6:
        usage()
        sys.exit(-1)

    go_file, goa_file, output_dir, list_prot, assoc1 = sys.argv[1:6]
    e = Evaluation(go_file, goa_file, output_dir, list_prot)

    if len(sys.argv) == 7:
        assoc2 = sys.argv[6]
        e.evaluate_two_files(assoc1, assoc2)
    elif len(sys.argv) == 5:
        e.evaluate_single_file(assoc1)
    else:
        usage()
        sys.exit(-1)
