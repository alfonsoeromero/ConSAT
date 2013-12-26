#!/usr/bin/env python
"""Automated script for running Gfam (or ConSAT) in a set of sequences without
having anything else than the sequences (and perhaps the InterPro output).

There are four use cases in this program:

    1.- Run GFam/ConSAT in UniProt: all the related UniProtKB files will be downloaded,
    including the InterPro data. The user will only have to provide a directory
    where the downloaded and the output files will be written.

    GFam will be run in the SwissProt proteins, and ConSAT will be run in the
    rest (TremBL). The models learnt by GFam in SwissProt will be used by ConSAT
    although the user can include a previously found set of models of putative
    new domains to extend this and maintain the names with respect to the old
    domains. The putative models found will be named GFAM00001, GFAM00002, etc.
    It is assumed that the provided models are also following this naming
    convention.

    2.- Runs only ConSAT on UniProt. All the resultad UniProtKB files will
    be downloaded including the InterPro data. The user will only have to
    provide a directory where the downloaded and the output files will be 
    written.

    ConSAT will be run in all proteins, and therefore the user will have to
    provide, as well a file with the HMMs.

    3.- Run GFam in a set of user-sequences plus InterPro output data. It is 
    assumed that the InterPro output corresponds to those sequences the user is
    inputing. A set of putative new domain models will be predicted for this set
    of sequences (NOVEL00001, NOVEL00002, etc.). This is equivalent to the
    ``classic`` GFam approach (with some minor tweaks) except of the fact that
    a set of HMM models for the new domains is provided.

    4.- Run ConSAT in a set of user-sequences plus InterPro output data. It is
    assumed, again, that the InterPro output corresponds to those sequences the
    user is inputing. Besides, the user __should__ provide a set of HMMs found
    in previous GFam executions to allow ConSAT tagging the unassigned pieces
    with these models.

All four cases will take a single positional argument, a route to a directory 
which will be created if it does not exists, and where three sub-directories
will be setup:
    data/ : here will be downloaded all the supporting files used to 
    compute our output. Also, the given files will be copied inside.

    work/ : work folder for ConSAT/GFam

    output/ : output folder for ConSAT/GFam

The GFam configuration file used will be created in the root directory of the
specified route.

In the first use case (UniProt), two configuration files (gfam.cfg and 
consat.cfg) will be created, and one additional directory level will be setup
in this case:
    data/ : all the input / given files

    work/gfam : work folder for GFam
    work/consat : work folder for ConSAT

    output/gfam : output folder for GFam
    output/consat : output folder for Consat

"""

from __future__ import with_statement

from gfam.scripts import CommandLineApp
from gfam.utils import redirected
from gfam.uniprot import XMLIprscanToTxt
import urllib2
import os
import gzip
import shutil
import logging
from subprocess import call
import sys
import math
import re

__author__  = "Alfonso E. Romero"
__email__   = "aromero@cs.rhul.ac.uk"
__copyright__ = "Copyright (c) 2013, Alfonso E. Romero"
__license__ = "GPL"

__all__ = ["AutomatedGFam"]

class AutomatedGFam(CommandLineApp):
    """\
        Usage: %prog [options] mode directory

            Performs an automated run of GFam/Consat in one of four `modes`:

                - uniprot
                - uniprot_consat
                - gfam
                - consat

            Each of the models will download all the required data in the
            specified `directory`

    1.- Run GFam/ConSAT on UniProt: all the related UniProtKB files will be downloaded,
    including the InterPro data. The user will only have to provide a directory
    where the downloaded and the output files will be written.

    GFam will be run in the SwissProt proteins, and ConSAT will be run in the
    rest (TremBL). The models learnt by GFam in SwissProt will be used by ConSAT
    although the user can include a previously found set of models of putative
    new domains to extend this and maintain the names with respect to the old
    domains. The putative models found will be named GFAM00001, GFAM00002, etc.
    It is assumed that the provided models are also following this naming
    convention.

    2.- Runs only ConSAT on UniProt. All the resultad UniProtKB files will
    be downloaded including the InterPro data. The user will only have to
    provide a directory where the downloaded and the output files will be 
    written.

    ConSAT will be run in all proteins, and therefore the user will have to
    provide, as well a file with the HMMs.

    3.- Run GFam in a set of user-sequences plus InterPro output data. It is 
    assumed that the InterPro output corresponds to those sequences the user is
    inputing. A set of putative new domain models will be predicted for this set
    of sequences (NOVEL00001, NOVEL00002, etc.). This is equivalent to the
    ``classic`` GFam approach (with some minor tweaks) except of the fact that
    a set of HMM models for the new domains is provided.

    4.- Run ConSAT in a set of user-sequences plus InterPro output data. It is
    assumed, again, that the InterPro output corresponds to those sequences the
    user is inputing. Besides, the user __should__ provide a set of HMMs found
    in previous GFam executions to allow ConSAT tagging the unassigned pieces
    with these models.

    All four cases will take a single positional argument, a route to a directory 
    which will be created if it does not exists, and where three sub-directories
    will be setup:

    data/ : here will be downloaded all the supporting files used to 
    compute our output. Also, the given files will be copied inside.

    work/ : work folder for ConSAT/GFam

    output/ : output folder for ConSAT/GFam

    The GFam configuration file used will be created in the root directory of the
    specfied route.

    In the first use case (UniProt), two configuration files (gfam.cfg and 
    consat.cfg) will be created, and one additional directory level will be setup
    in this case:

    data/ : all the input / given files

    work/gfam : work folder for GFam
    work/consat : work folder for ConSAT

    output/gfam : output folder for GFam
    output/consat : output folder for Consat

    """

    urls = {"file.mapping.gene_ontology": "http://purl.obolibrary.org/obo/go/go-basic.obo",
            "file.mapping.interpro2go" : "http://www.geneontology.org/external2go/interpro2go",
            "file.mapping.interpro_parent_child" : "ftp://ftp.ebi.ac.uk/pub/databases/interpro/ParentChildTreeFile.txt",
           }

    params = dict() 

    short_name = "automated"

    def __init__(self, *args, **kwds): 
        super(AutomatedGFam, self).__init__(*args, **kwds)
        #Use the logger from Modula
        self.log = logging.getLogger("automated") 

    def create_parser(self):
        """Creates the command line parser for the automated script"""
        parser = super(AutomatedGFam, self).create_parser()
        parser.add_option("-f", "--force", dest="force", action="store_true",
                help="force recalculation of results even when gfam thinks "\
                     "everything is up-to-date")
        parser.add_option("-m", "--models", dest="models", metavar="FILE",
                help="file with previously computed HMM models to maintain"
                "nomenclature (uniprot/gfam mode) or use them (consat mode)")
        parser.add_option("-s", "--user-sequences", dest="sequences", 
                metavar="FILE", help="file with user sequences (gfam/consat"\
                        " modes)")
        parser.add_option("-i", "--interpro", dest="interpro", metavar="FILE",
                help="interpro data for gfam/consat modes")
        parser.add_option("-r", "--regexp", dest="regexp", metavar="REGEX",
                help="regular expression to extract protein ids (gfam/consat)")
        parser.add_option("-n", "--no-download", dest="no_download", 
                action="store_true", help="does not download data, use if you"\
                        "want to continue a previous execution")
        parser.add_option("-g", "--goa", dest="goa", help="GOA file to transfer"\
                         "the function from")
        parser.add_option("-b", "--blast-route", dest="blast_route", help="route to"\
                        "blast executables if not in the path")
        parser.add_option("-p", "--pubmed-cache", dest="pubmed_cache", help="route to"\
                        " the cache of PubMed downloaded files")
        parser.add_option("-o", "--gene-ontology", dest="gene_ontology", help="route to"\
                        " the desired Gene Ontology file (instead of using the most recent one")
        parser.add_option("-x", "--lexicon", dest="lexicon", help="route to lexicon "\
                        "file to maintain ids compatibility")
        return parser

    def process_arguments(self):
        if len(self.args) != 2:
            self.error("Wrong argument number")

        self.mode, self.directory = self.args

        if not self.mode in ["uniprot", "uniprot_consat", "gfam", "consat"]:
            self.error("Wrong mode, must be 'uniprot', 'uniprot_consat', 'gfam' or 'consat'")

        if self.options.models:
            self.models = self.options.models
        elif self.mode == "consat" or self.mode == "uniprot_consat":
            self.error(self.mode + " mode should have a 'models' file (option -m)")

        if not self.mode.startswith("uniprot") and not self.options.sequences:
            self.error("No sequence file specified (-s option)")

        if not self.mode.startswith("uniprot") and not self.options.interpro:
            self.error("InterPro sequences were not specified (-i option)")

    def _download_file(self, url, filename):
        """
        Downloads a file from a given `url` to an specified `filename`
        route. It downloads first to a "temporary" file which is the same
        `filename` ended by "_tmp", which is finally renamed to the desired
        name.
        """
        if not os.path.isfile(filename) or self.options.force:
            print "Downloading from ", url
            f = urllib2.urlopen(url)
            block_sz = 8192
            with open(filename + "_tmp", "wb") as output:
                while True:
                    chunk = f.read(block_sz)
                    if not chunk:
                        break
                    output.write(chunk)
            os.rename(filename + "_tmp", filename)
            print "...Downloaded"

    def init_folders(self):
        self.log.info("Creating directories...")
        # 1.- create directory
        if not os.path.exists(self.directory):
            os.makedirs(self.directory)

        # 2.- create data, work and output dirs
        self.program_dir = dict()
        for dirname in ["data", "work", "output"]:
            self.program_dir[dirname] = os.path.join(self.directory, dirname)
            if not os.path.exists(self.program_dir[dirname]):
                os.makedirs(self.program_dir[dirname])

        if self.options.blast_route:
            self.params["folder.blast"] = self.options.blast_route

        if self.options.pubmed_cache:
            self.params["folder.pubmed_cache"] = self.options.pubmed_cache
        else:
            self.params["folder.pubmed_cache"]=os.path.join(self.program_dir["data"], "pubmed_cache")

        if self.mode == "uniprot":
            self.gfam_dir = dict()
            self.consat_dir = dict()
            for subdir, route in self.program_dir.items():
                self.gfam_dir[subdir] = os.path.join(route, "gfam")
                if not os.path.exists(self.gfam_dir[subdir]):
                    os.makedirs(self.gfam_dir[subdir])
                self.consat_dir[subdir] = os.path.join(route, "consat")
                if not os.path.exists(self.consat_dir[subdir]):
                    os.makedirs(self.consat_dir[subdir])

    def _gunzip(self, file_name):
        file_name2 = file_name.replace('.gz', '')
        if os.path.exists(file_name2):
            return
        file_name2_tmp = file_name2 + "_tmp"
        f = gzip.open(file_name, 'rb')
        with open(file_name2_tmp, "w") as out:
            for line in f:
                out.write(line)
        f.close()
        os.rename(file_name2 + "_tmp", file_name2)

    def download_data(self):
        """Downloads the required data files for running Gfam
        and Consta into the "data/" folder. This includes copying 
        the user sequences and InterPro data into that folder for
        the case of GFam and Consat executions.
        """
        self.log.info("Downloading data files...")
        # 1.- download general files
        data = self.program_dir["data"]
        for name, url in self.urls.items():
            filename = os.path.join(data, url.split("/")[-1])
            if name != "file.mapping.gene_ontology" or (name == "file.mapping.gene_ontology" and not self.options.gene_ontology):
                self._download_file(url, filename)
                self.params[name] = filename

        if self.options.gene_ontology:
            self.params["file.mapping.gene_ontology"] = self.options.gene_ontology

        print "YEHEEEEEE ", self.params["file.mapping.interpro_parent_child"]

        # 2.- execute the "download_names" script
        names = os.path.join(data, "names.dat.gz")
        self.params["file.mapping.interpro2name"] = names
        if not os.path.isfile(names) or self.options.force:
            import tempfile
            f = tempfile.NamedTemporaryFile(delete=False)
            f.close()
            from gfam.scripts.download_names import DownloadNamesApp
            stdout = sys.stdout
            with open(f.name, 'w') as sys.stdout:
                script = DownloadNamesApp()
                script.run_real()
            sys.stdout = stdout
            with open(f.name, 'rb') as orig_file:
                with gzip.open(names, 'wb') as zipped_file:
                    zipped_file.writelines(orig_file)

        if self.mode == "uniprot":
            # we download UniProt data: fasta files
            self.params_gfam = dict()
            self.params_consat = dict()
            #######################################
            self.log.info("Downloading sequences...")
            url_swissprot = "ftp://ftp.ebi.ac.uk/pub/databases/uniprot/knowledgebase/uniprot_sprot.fasta.gz"
            file_swissprot = os.path.join(self.program_dir["data"], "uniprot_sprot.fasta.gz")
            self._download_file(url_swissprot, file_swissprot)
            self._gunzip(file_swissprot)
            file_swissprot = file_swissprot.replace(".gz", "")

            url_trembl = "ftp://ftp.ebi.ac.uk/pub/databases/uniprot/knowledgebase/uniprot_trembl.fasta.gz"
            file_trembl = os.path.join(self.program_dir["data"], "uniprot_trembl.fasta.gz")
            self._download_file(url_trembl, file_trembl)
            self._gunzip(file_trembl)
            file_trembl = file_trembl.replace(".gz", "")

            file_swissprot_base = os.path.basename(file_swissprot)
            filename = os.path.join(self.gfam_dir["data"], file_swissprot_base)
            if not os.path.isfile(filename) or self.options.force:
                shutil.copy(file_swissprot, self.gfam_dir["data"])
            self.params_gfam["file.input.sequences"] = os.path.join(self.gfam_dir["data"], "uniprot_sprot.fasta")
            out = self._create_mask_fasta_file(self.params_gfam["file.input.sequences"])
            self.params_gfam["low_complexity_regions_file"] = out

            file_trembl_base = os.path.basename(file_trembl)
            filename = os.path.join(self.consat_dir["data"], file_trembl_base)
            if not os.path.isfile(filename) or self.options.force:
                shutil.copy(file_trembl, self.consat_dir["data"])
            self.params_consat["file.input.sequences"] = os.path.join(self.consat_dir["data"], "uniprot_trembl.fasta")
            out = self._create_mask_fasta_file(self.params_consat["file.input.sequences"])
            self.params_consat["low_complexity_regions_file"] = out

            # we write in a file the set of swissprot ids
            file_sprot_ids = os.path.join(self.gfam_dir["data"], "swissprot_ids")
            with open(file_sprot_ids, "w") as out:
                for line in open(self.params_gfam["file.input.sequences"], "r"):
                    if line.startswith(">"):
                        id_sprot = line.split("|")[1]
                        out.write(id_sprot + "\n")

            self.log.info("Downloading interpro...")

            # we download the InterPro data and prepare it
            url_interpro = "ftp://ftp.ebi.ac.uk/pub/databases/interpro/match_complete.xml.gz"
            file_interpro = os.path.join(self.program_dir["data"], "match_complete.xml.gz")
            self._download_file(url_interpro, file_interpro)
            self._gunzip(file_interpro)
            file_interpro = file_interpro.replace(".gz", "")

            interpro_gfam = os.path.join(self.gfam_dir["data"], "uniprot_sprot.interpro")
            self.params_gfam["file.input.iprscan"] = interpro_gfam
            interpro_consat = os.path.join(self.consat_dir["data"], "uniprot_trembl.interpro")
            self.params_consat["file.input.iprscan"] = interpro_consat

            self.log.info("Partitioning InterPro file...")

            # Conversion of the intepro file (XML) into text
            if not os.path.isfile(interpro_gfam) and not os.path.isfile(interpro_consat):
                conversor = XMLIprscanToTxt(file_interpro, self.params["file.mapping.interpro2go"],
                        file_sprot_ids, interpro_gfam, interpro_consat)
                conversor.run()
    
            # work and output folders
            self.params_gfam["folder.work"] = self.gfam_dir["work"]
            self.params_gfam["folder.output"] = self.gfam_dir["output"]
            self.params_consat["folder.work"] = self.consat_dir["work"]
            self.params_consat["folder.output"] = self.consat_dir["output"]

            # GOA
            url_goa = "ftp://ftp.ebi.ac.uk/pub/databases/GO/goa/UNIPROT/gene_association.goa_uniprot.gz"
            file_goa = os.path.join(self.program_dir["data"], "gene_association.goa_uniprot.gz")
            self._download_file(url_goa, file_goa)
            self._gunzip(file_goa)
            file_goa = file_goa.replace(".gz", "")
            self.params["file.function.goa_file"] = file_goa

            self.params["sequence_id_regexp"] = "(\w+\|)(?P<id>\w+)(\|\w+)+"
            self.params["num_cpu_cores"] = "12"

            if self.options.models:
                shutil.copy(self.options.models, data)
                name_models = os.path.basename(self.options.models)
                file_models = os.path.join(data, name_models)
                self.params_gfam["previous_domain_table"] = file_models

            computed_models = os.path.join(self.params_gfam["folder.output"], "hmms")
            computed_models = os.path.join(computed_models, "all.hmm")
            self.params_consat["file.input.hmms"] = computed_models
            self.params_consat["transfer_from_both"] = "True"

            old_lexicon = os.path.join(self.params_gfam["folder.output"], "lexicon")
            self.params_consat["file.lexicon"] = old_lexicon
            table = os.path.join(self.params_gfam["folder.output"], "domain_architectures.tab")
            self.params_consat["file.function.domain_arch_table"] = table

            # we download the RDF data
            url_rdf = "http://www.uniprot.org/citations/?query=citedin%3a(*)&compress=yes&format=rdf"
            rdf_file = os.path.join(self.program_dir["data"], "citations.rdf.gz")
            self._download_file(url_rdf, rdf_file)
            self._gunzip(rdf_file)
            rdf_file = rdf_file.replace(".gz", "")
            self.params["file.rdffile"] = rdf_file

            # we download the Idmapping data
            url_id_mapping = "ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/idmapping/idmapping_selected.tab.gz"
            id_mapping_file = os.path.join(self.program_dir["data"], "idmapping_selected.tag.gz")
            self._download_file(url_id_mapping, id_mapping_file)
            self._gunzip(id_mapping_file)
            id_mapping_file = id_mapping_file.replace(".gz", "")
            self.params["file.idmapping"] = id_mapping_file

            # we set the prefix to "GFAM"
            self.params["prefix"] = "GFAM"

        elif self.mode == "consat" or self.mode == "gfam": 
            # gfam or consat modes, we first copy two input files
            # into the data folder
            shutil.copy(self.options.sequences, data)
            name_seq = os.path.basename(self.options.sequences)
            self.params["file.input.sequences"] = os.path.join(data, name_seq)
            shutil.copy(self.options.interpro, data)
            name_ipr = os.path.basename(self.options.interpro)
            self.params["file.input.iprscan"] = os.path.join(data, name_ipr)
            out = self._create_mask_fasta_file(self.params["file.input.sequences"])
            self.params["low_complexity_regions_file"] = out
            if self.options.goa:
                shutil.copy(self.options.goa, data)
                name_goa = os.path.basename(self.otions.goa)
                self.params["file.function.goa_file"] = os.path.join(data, name_goa)
            if self.options.regexp:
                self.params["sequence_id_regexp"] = self.options.regexp
            if self.options.models:
                shutil.copy(self.options.models, data)
                name_models = os.path.basename(self.options.models)
                file_models = os.path.join(data, name_models)
                if self.mode == "gfam":
                    self.params["previous_domain_table"] = file_models
                else:
                    self.params["file.input.hmms"] = file_models 
            self.params["folder.work"] = self.program_dir["work"]
            self.params["folder.output"] = self.program_dir["output"]
        else: # mode == "uniprot_consat"
            # 1.- download the fasta files and join them
            self.log.info("Downloading sequences...")
            url_swissprot = "ftp://ftp.ebi.ac.uk/pub/databases/uniprot/knowledgebase/uniprot_sprot.fasta.gz"
            file_swissprot = os.path.join(self.program_dir["data"], "uniprot_sprot.fasta.gz")
            self._download_file(url_swissprot, file_swissprot)
            self._gunzip(file_swissprot)
            file_swissprot = file_swissprot.replace(".gz", "")

            url_trembl = "ftp://ftp.ebi.ac.uk/pub/databases/uniprot/knowledgebase/uniprot_trembl.fasta.gz"
            file_trembl = os.path.join(self.program_dir["data"], "uniprot_trembl.fasta.gz")
            self._download_file(url_trembl, file_trembl)
            self._gunzip(file_trembl)
            file_trembl = file_trembl.replace(".gz", "")

            sequences = os.path.join(data, "uniprot.fasta")
            with open(sequences, 'w') as outfile:
                for fname in (file_swissprot, file_trembl):
                    with open(fname) as infile:
                        for line in infile:
                            outfile.write(line)
            self.params["file.input.sequences"] = sequences

            # 2.- download an prepare the InterPro files
            self.log.info("Downloading interpro...")

            # we download the InterPro data and prepare it
            url_interpro = "ftp://ftp.ebi.ac.uk/pub/databases/interpro/match_complete.xml.gz"
            file_interpro = os.path.join(self.program_dir["data"], "match_complete.xml.gz")
            self._download_file(url_interpro, file_interpro)
            self._gunzip(file_interpro)
            file_interpro = file_interpro.replace(".gz", "")

            # we write in a file the set of swissprot ids
            file_sprot_ids = os.path.join(self.program_dir["data"], "swissprot_ids")
            with open(file_sprot_ids, "w") as out:
                for line in open(file_swissprot, "r"):
                    if line.startswith(">"):
                        id_sprot = line.split("|")[1]
                        out.write(id_sprot + "\n")

            interpro1 = os.path.join(self.program_dir["data"], "interpro1")
            interpro2 = os.path.join(self.program_dir["data"], "interpro2")
            interpro_out = os.path.join(self.program_dir["data"], "uniprot.interpro")

            # Conversion of the intepro file (XML) into text
            if not os.path.isfile(interpro1) and not os.path.isfile(interpro2):
                conversor = XMLIprscanToTxt(file_interpro, self.params["file.mapping.interpro2go"],
                        file_sprot_ids, interpro1, interpro2)
                conversor.run()
            with open(sequences, 'w') as outfile:
                for fname in (interpro1, interpro2):
                    with open(fname) as infile:
                        for line in infile:
                            outfile.write(line)
            self.params["file.input.iprscan"] = interpro_out

            out = self._create_mask_fasta_file(self.params["file.input.sequences"])
            self.params["low_complexity_regions_file"] = out
            self.params["sequence_id_regexp"] = "(\w+\|)(?P<id>\w+)(\|\w+)+"
            self.params["num_cpu_cores"] = "12"
    
            # GOA
            url_goa = "ftp://ftp.ebi.ac.uk/pub/databases/GO/goa/UNIPROT/gene_association.goa_uniprot.gz"
            file_goa = os.path.join(self.program_dir["data"], "gene_association.goa_uniprot.gz")
            self._download_file(url_goa, file_goa)
            self._gunzip(file_goa)
            file_goa = file_goa.replace(".gz", "")
            self.params["file.function.goa_file"] = file_goa

            if self.options.models:
                # should always be true...
                shutil.copy(self.options.models, data)
                name_models = os.path.basename(self.options.models)
                file_models = os.path.join(data, name_models)
                self.params["file.input.hmms"] = file_models
                
            if self.options.lexicon:
                # lexicon, to maintain ids compatibility
                shutil.copy(self.options.lexicon, data)
                self.params["file.lexicon"] = os.path.join(data, os.path.basename(self.options.lexicon))

            # we download the RDF data
            url_rdf = "http://www.uniprot.org/citations/?query=citedin%3a(*)&compress=yes&format=rdf"
            rdf_file = os.path.join(self.program_dir["data"], "citations.rdf.gz")
            self._download_file(url_rdf, rdf_file)
            self._gunzip(rdf_file)
            rdf_file = rdf_file.replace(".gz", "")
            self.params["file.rdffile"] = rdf_file

            # we download the Idmapping data
            url_id_mapping = "ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/idmapping/idmapping_selected.tab.gz"
            id_mapping_file = os.path.join(self.program_dir["data"], "idmapping_selected.tag.gz")
            self._download_file(url_id_mapping, id_mapping_file)
            self._gunzip(id_mapping_file)
            id_mapping_file = id_mapping_file.replace(".gz", "")
            self.params["file.idmapping"] = id_mapping_file

            # general folders
            self.params["folder.work"] = self.program_dir["work"]
            self.params["folder.output"] = self.program_dir["output"]

    def _create_mask_fasta_file(self, filename):
        """Creates a mask of a fasta file by running the `seg` utility
        on it. The results are stored in the same directory in a file with
        the same name but with the "mask" extension.
        """
        out = filename[0:filename.rfind('.')] + ".mask"
        if not os.path.isfile(out):
            call(["segmasker", "-in", filename, "-out", out])
        return out

    def create_conf_file(self):
        self.conf_file = []
        if self.mode == "gfam" or self.mode == "uniprot":
            from gfam.scripts.master import GFamMasterScript
            script = GFamMasterScript()
            conf_file = os.path.join(self.directory, "gfam.conf")
            if os.path.isfile(conf_file): 
                os.remove(conf_file)
            self.conf_file.append(conf_file)
            script.run(args=["init", "-s", "-c", conf_file])
            if self.mode == "uniprot":
                self._fill_config_file(conf_file, self.params_gfam)
            else:
                self._fill_config_file(conf_file)
        if self.mode == "consat" or self.mode == "uniprot" or self.mode == "uniprot_consat":
            from gfam.scripts.master_consat import ConSATMasterScript
            script = ConSATMasterScript()
            conf_file = os.path.join(self.directory, "consat.conf")
            if os.path.isfile(conf_file):
                os.remove(conf_file)
            self.conf_file.append(conf_file)
            script.run(args=["init", "-s", "-c", conf_file])
            self._fill_config_file(conf_file)
            if self.mode == "uniprot":
                self._fill_config_file(conf_file, self.params_consat)
            else:
                self._fill_config_file(conf_file)

    def _fill_config_file(self, filename, my_params=None):
        tmp_filename = filename + "_tmp"
        with open (tmp_filename, "w") as out:
            for line in open(filename, "r"):
                if "=" in line:
                    param = line.split("=")[0].strip()
                    if param in self.params:
                        value = self.params[param]
                        out.write(param + "=" + value + "\n")
                    elif my_params is not None and param in my_params:
                        value = my_params[param]
                        out.write(param + "=" + value + "\n")
                    else:
                        # pre-defined parameters
                        out.write(line)
                else:
                    out.write(line)
                    
        os.rename(tmp_filename, filename)

    def run_calculation(self):
        if self.mode == "gfam":
            self.log.info("Running GFam on users' sequences")
            from gfam.scripts.master import GFamMasterScript
            script = GFamMasterScript()
            script.run(args=["run", "-c", self.conf_file[0]])
        elif self.mode == "consat":
            self.log.info("Running ConSAT on users' sequences")
            from gfam.scripts.master_consat import ConSATMasterScript
            script = ConSATMasterScript()
            script.run(args=["run", "-c", self.conf_file[0]])
        elif self.mode == "uniprot_consat":
            self.log.info("Running ConSAT on UniProt")
            from gfam.scripts.master_consat import ConSATMasterScript
            script = ConSATMasterScript()
            script.run(args=["run", "-c", self.conf_file[0]])
        else:
            self.log.info("Running GFam on UniProtKB")
            from gfam.scripts.master import GFamMasterScript
            script = GFamMasterScript()
            script.run(args=["run", "-c", self.conf_file[0]])

            self.log.info("Running ConSAT on UniProtKB")
            from gfam.scripts.master_consat import ConSATMasterScript
            script2 = ConSATMasterScript()
            script2.run(args=["run", "-c", self.conf_file[1]])

    def prepare_package(self):
        pass

    def doc_freq_per_id(self, filename):
        doc = dict()
        pepe = False
        for line in open(filename, 'r'):
            _, id, df = line.split()
            doc[int(id)] = int(df)
        return doc

    def get_num_proteins(self, filename):
        N = 0
        for line in open(filename, 'r'):
            if line.startswith(">"):
                N+=1
        return N

    def get_unnormalized_vector(self, line, idf):
        v = dict()
        for token in line.split():
            id, weight = token.split(':')
            v[int(id)] = float(weight)/idf[int(id)]
        return v

    def join_text_representation(self):
        self.log.info("Joining UniProtKB (sp + tr) textual annotation in a single file")

        N_sp = self.get_num_proteins(self.params_gfam["file.input.sequences"])
        N_tr = self.get_num_proteins(self.params_consat["file.input.sequences"])
        N = N_sp + N_tr # Number of "documents" = Total number of proteins
        lexicon_sp = os.path.join(self.params_gfam["folder.output"], "lexicon")
        df_sp = self.doc_freq_per_id(lexicon_sp)
        lexicon_tr = os.path.join(self.params_consat["folder.output"], "lexicon")
        df_tr = self.doc_freq_per_id(lexicon_tr)
        df_total = dict(df_sp)
        for id, df in df_tr.items():
            if id in df_total:
                df_total[id] = df_total[id] + df
            else:
                df_total[id] = df
        idf_sp = dict([(id, math.log(N_sp/df)) for id, df in df_sp.items()])
        idf_tr = dict([(id, math.log(N_tr/df)) for id, df in df_tr.items()])
        idf = dict([(id, math.log(N/df)) for id, df in df_total.items()])

        text_file_sp = os.path.join(self.params_gfam["folder.output"], "weight_file_per_arch")
        text_file_tr = os.path.join(self.params_consat["folder.output"], "weight_file_per_arch")
        out_file = os.path.join(self.params_consat["folder.output"], "weight_file_per_arch_sp_tr")

        archs = dict()
        archs_sp = set()
        pattern = re.compile('\s+')
        for _line in open(text_file_sp):
            line = _line.strip()
            if pattern.findall(line):
                arch, vector = line.split(None, 1)
                archs_sp.add(arch)
            else:
                arch, vector = line, ""
            archs[arch] = self.get_unnormalized_vector(vector, idf_sp)

        visited = set()
        with open(out_file, "w") as out:
            for _line in open(text_file_tr):
                line = _line.strip()
                if pattern.findall(line):
                    arch, vector = line.split(None, 1)
                else:
                    arch, vector = line.strip(), ""
                v = self.get_unnormalized_vector(vector, idf)
                if arch in archs:
                    for term, weight in archs[arch].items():
                        if term in v:
                            v[term] = v[term] + weight
                        else:
                            v[term] = weight
                v_idf = dict([(id, idf[id]*w) for id, w in v.items()])
                out.write(arch + " ")
                out.write(" ".join(["{0}:{1}".format(id, w) for id, w in v_idf.items()]))
                out.write("\n")
                visited.add(arch)
            remaining = set(archs.keys()) - visited
            for arch in remaining:
                out.write(arch + " ")
                v = archs[arch]
                v_idf = dict([(id, idf[id]*w) for id, w in v.items()])
                out.write(" ".join(["{0}:{1}".format(id, w) for id, w in v_idf.items()]))
                out.write("\n")

    def run_real(self):
        """Runs the application"""
        self.process_arguments()
        self.init_folders()
        self.download_data()
        self.create_conf_file()
        self.run_calculation()
        if self.mode == "uniprot":
            self.prepare_package()
            self.join_text_representation()
        print "All things finished!"

