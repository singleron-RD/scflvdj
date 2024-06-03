#!/usr/bin/env python

import argparse

import parse_protocol
import pyfastx
import utils

ASSAY = "scatac"


class Auto(parse_protocol.Auto):
    def __init__(
        self,
        fq1_list,
        sample,
    ):
        super().__init__(fq1_list, sample)

    def seq_protocol(self, seq):
        """
        Returns: protocol or None

        >>> runner = Auto([], "fake_sample")
        """
        for protocol in self.protocol_dict:
            if self.is_protocol(seq, protocol):
                return protocol


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--sample", required=True)
    parser.add_argument("--fq1", required=True)
    parser.add_argument("--fq2", required=True)
    parser.add_argument("--fq3", required=True)
    parser.add_argument("--assets_dir", required=True)
    parser.add_argument("--protocol", required=True)
    args = parser.parse_args()

    fq1_list = args.fq1.split(",")
    fq2_list = args.fq2.split(",")
    fq3_list = args.fq3.split(",")
    # protocol
    protocol_dict = parse_protocol.get_protocol_dict(args.assets_dir)
    if args.protocol == "auto":
        p = Auto(fq1_list, args.sample).run()
    else:
        p = args.protocol

    fn = f"{args.sample}.{ASSAY}.protocol.stats.json"
    utils.write_json({"Protocol": p}, fn)

    pmeta = protocol_dict[p]
    pattern_dict = pmeta["pattern_dict"]
    raw_list, mismatch_list = parse_protocol.get_raw_mismatch(pmeta["bc"], 1)

    # out_fq
    out_fq_fn = {x: f"{args.sample}_R{x}.fq.gz" for x in [1, 2]}
    outdict = {k: utils.openfile(v, "wt") for k, v in out_fq_fn.items()}

    raw_reads = 0
    valid_reads = 0
    corrected_reads = 0
    for fq1, fq2 in zip(fq1_list, fq2_list):
        fq1 = pyfastx.Fastx(fq1)
        fq2 = pyfastx.Fastx(fq2)

        for (name1, seq1, qual1), (name2, seq2, qual2) in zip(fq1, fq2):
            raw_reads += 1
            bc_list = [seq2[x] for x in pattern_dict["C"]]
            valid, corrected, corrected_seq = parse_protocol.check_seq_mismatch(bc_list, raw_list, mismatch_list)
            if valid:
                valid_reads += 1
                if corrected:
                    corrected_reads += 1
                read_name = f"{corrected_seq}:{raw_reads}"
                umi = parse_protocol.get_seq_str(seq1, pattern_dict["U"])
                qual_umi = parse_protocol.get_seq_str(qual1, pattern_dict["U"])
                qual_bc = parse_protocol.get_seq_str(qual1, pattern_dict["C"])
                bc = corrected_seq.replace("_", "")
                outdict[1].write(utils.fastq_str(read_name, bc + umi, qual_bc + qual_umi))
                outdict[2].write(utils.fastq_str(read_name, seq2, qual2))

    fn = f"{args.sample}.{ASSAY}.extract.stats.json"
    metrics = {"Raw Reads": raw_reads}
    metrics["Valid Reads"] = utils.get_frac(valid_reads / raw_reads)
    metrics["Corrected Barcodes"] = utils.get_frac(corrected_reads / valid_reads)
    utils.write_json(metrics, fn)
