"""Utility classes for parsing UniProtKB files and related ones."""

from __future__ import print_function
import datetime
from collections import defaultdict
import re
import sys

__author__ = "Alfonso E. Romero"
__email__ = "aeromero@cs.rhul.ac.uk"
__copyright__ = "Copyright (c) 2013, Tamas Nepusz"
__license__ = "GPL"
__all__ = ["XMLIprscanToTxt"]


class XMLIprscanToTxt(object):
    """\
        Usage: %prog [options] INPUT_XML_FILE INTERPRO2GOFILE
            LIST_SWISSPROT_IDS OUTPUT_TXT_FILE_SWISSPROT
            OUTPUT_TXT_FILE_TREMBL

        Converts an InterProScan XML output to a TXT one, separated for
        SwissProt and TremBL. The XML file can be dowloaded from:
        ftp://ftp.ebi.ac.uk/pub/databases/interpro/match_complete.xml.gz

        Note that this program does not really "parse" the XML file, as
        its format is relatively easy. Instead of using a SAX parser, we
        have written a "cooked" XML parser based on regexes. For more
        details on which is the format of that XML, please check the DTD
        at ftp://ftp.ebi.ac.uk/pub/databases/interpro/match_complete.dtd
    """

    fields = ["protein_id", "protein_crc64", "protein_length", "match_evd",
              "match_id", "match_name", "lcn_start", "lcn_end",
              "lcn_score", "match_status", "date", "ipr_id", "ipr_name"]

    def __init__(self, xml_file, interpro2go, list_file, out_sp, out_tr):
        self.xml_file = xml_file
        self.interpro2go = interpro2go
        self.list_file = list_file
        self.out_sp = out_sp
        self.out_tr = out_tr
        now = datetime.datetime.now()
        self.date = now.strftime("%d-%b-%y").upper()
        self.match_data = {"date": self.date}
        self.header_is_read = False
        self.ipr2go = defaultdict(list)
        self.swissprot = set()
        self.attr = ""
        self.sp_file = ""
        self.tr_file = ""

    def _read_interpro2go(self):
        self.ipr2go = defaultdict(list)
        for line in open(self.interpro2go, "r"):
            if line.startswith("!"):
                continue
            ipr_term = line.split(":")[1].split()[0]
            go_term = line.split(">")[1].strip()
            self.ipr2go[ipr_term].append(go_term)

    def _read_protein_ids(self):
        self.swissprot = set()
        for line in open(self.list_file, "r"):
            self.swissprot.add(line.strip())

    def _clear_all(self):
        for field in self.fields:
            self.match_data[field] = ""
        self.match_data["date"] = self.date

    def _clear_data_for_tag(self, tag):
        for field in self.fields:
            if tag in field:
                self.match_data[field] = ""

    def _print_file_line(self):
        line = ["\t".join([self.match_data[val] for val in self.fields])]
        if "ipr_id" in self.match_data:
            line.append("\t")
            line.append(", ".join(self.ipr2go[self.match_data["ipr_id"]]))
        line.append("\n")
        # print to a file
        if "mobidb" in self.match_data["protein_id"].strip():
            return
        if self.match_data["protein_id"].strip() in self.swissprot:
            self.sp_file.write(''.join(line))
        else:
            self.tr_file.write(''.join(line))

    def _process_xml_tag(self, tag):
        """Process an entire XML tag from the InterPro
        match_complete.xml file

        `string` : string containing the tag to process
        """
        # two cases, opening or closing tag
        if tag[0] == '<' and tag[1] != '/' and tag.endswith(">"):
            # opening tag
            if self.header_is_read:
                # we process it only if the header has been read
                interior = tag.split("<")[1].split(">")[0].strip()
                tag_name = interior.split()[0]
                content = interior.split(" ", 1)[1].strip()
                for result in re.finditer(self.attr, content):
                    result_set = result.group()\
                                       .strip()\
                                       .replace("\"", "")\
                                       .split("=")
                    key, val = result_set
                    value_set = val
                    self.match_data[tag_name + "_" + key.strip()] = value_set

                if tag_name == "lcn":
                    self._print_file_line()
        elif tag.startswith("</") and tag.endswith(">"):
            # closing tag
            tag_name = tag.split("/")[1].split()[0].replace('>', '')
            if tag_name == "release":
                self.header_is_read = True
            elif tag_name == "match":
                for clear_tag in ["ipr", "lcn", "match"]:
                    self._clear_data_for_tag(clear_tag)
        elif self.header_is_read:
            print("ERROR in XML file, tag={}".format(tag))
            sys.exit(-1)

    def run(self):
        """ Does everything, call after the constructor
        """
        self._clear_all()
        self._read_interpro2go()
        self._read_protein_ids()

        self.sp_file = open(self.out_sp, "w")
        self.tr_file = open(self.out_tr, "w")

        self.attr = re.compile(r'\w+=\"[^"]+\"')
        pattern = re.compile(r'<[^>]*>')
        for line in open(self.xml_file, "r"):
            for result in re.finditer(pattern, line):
                self._process_xml_tag(result.group().strip())

        self.sp_file.close()
        self.tr_file.close()
