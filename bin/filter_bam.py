#!/usr/bin/env python

"""
- remove reads that are not cell-associated.
- remove reads that are not mapped to target genes.
- remove PCR duplicate reads
- add read group tag RG
"""

import argparse
from collections import defaultdict

import pysam
import utils


class FilterBam:
    def __init__(self, args):
        self.args = args
        self.match_barcode = set(utils.read_one_col(args.match_barcode_file))
        self.genes = set(args.genes.split(","))
        self.out_bam_file = f"{args.sample}.filtered.bam"
        self.dup_dict = utils.nested_defaultdict(dim=4)
        self.gene_read = defaultdict(int)
        self.cb_useRead = defaultdict(int)

    def filter(self, max_dup=1):
        """
        for each (barcode,UMI,reference_name,reference_start), keep at most max_duplicate_reads
        """
        with pysam.AlignmentFile(self.args.bam, "rb") as reader:
            header = reader.header.to_dict()
            # add RG to header
            header["RG"] = []
            for cb in self.match_barcode:
                header["RG"].append(
                    {
                        "ID": cb,
                        "SM": cb,
                    }
                )
            with pysam.AlignmentFile(self.out_bam_file, "wb", header=header) as writer:
                for record in reader:
                    gn = record.get_tag("GN")
                    cb = record.get_tag("CB")
                    umi = record.get_tag("UB")
                    if any(x == "-" for x in (gn, cb, umi)):
                        continue
                    if (cb not in self.match_barcode) or (gn not in self.genes):
                        continue
                    rn, rs = record.reference_name, record.reference_start
                    self.dup_dict[cb][umi][rn][rs] += 1
                    self.gene_read[gn] += 1
                    if self.dup_dict[cb][umi][rn][rs] > max_dup:
                        continue
                    self.cb_useRead[cb] += 1
                    record.set_tag(tag="RG", value=record.get_tag("CB"), value_type="Z")
                    writer.write(record)

    def write_stats(self):
        gene_fn = f"{self.args.sample}.scflvdj.gene.json"
        stats_fn = f"{self.args.sample}.scflvdj.filter_bam.stats.json"
        utils.write_json(self.gene_read, gene_fn)

        x = sum(self.cb_useRead.values()) / len(self.cb_useRead)
        stats = {"Mean Used Reads per Cell": x}
        d = self.dup_dict
        total = 0
        distinct = 0
        for cb in d:
            for umi in d[cb]:
                for rn in d[cb][umi]:
                    for rs in d[cb][umi][rn]:
                        distinct += 1
                        total += d[cb][umi][rn][rs]
        saturation = 1 - float(distinct) / total
        stats.update({"Saturation": saturation})

        utils.write_json(stats, stats_fn)

    def run(self):
        self.filter()
        self.write_stats()


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--bam", help="Input bam file", required=True)
    parser.add_argument("--match_barcode_file", help="File containing matched barcodes", required=True)
    parser.add_argument("--sample", help="Sample name", required=True)
    parser.add_argument("--genes", help="Target genes")
    args = parser.parse_args()

    FilterBam(args).run()
