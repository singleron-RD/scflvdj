#!/usr/bin/env python

import argparse
from collections import defaultdict

import pandas as pd
import pysam
from sccore import utils


def extract_fields(record, samples):
    """
    extract fields from a snpeff annotated VCF record; ANN columns are 0-based
    CHROM, POS, REF, ALT, ANN[*].GENE(col 3), ANN[*].HGVS_C(col 9), ANN[*].HGVS_P(col 10)
    https://pcingola.github.io/SnpEff/snpsift/extractfields/
    Returns:
        name: str
        genotypes: list of str,e.g. ['0/0', '0/1', '1/1', 'NA']

    """
    ann = record.info["ANN"][0].split("|")
    effect = ann[1]
    impact = ann[2]
    gene = ann[3]
    rank = ann[8]
    hgvs_c = ann[9]
    hgvs_p = ann[10]
    name = "-".join([gene, hgvs_c, hgvs_p])

    genotypes = []
    for sample in samples:
        genotype = record.samples[sample]["GT"]
        if genotype == (None,):
            genotype_str = "NA"
        else:
            g1, g2 = genotype
            genotype_str = "/".join([str(g1), str(g2)])
        genotypes.append(genotype_str)

    meta = {
        "Effect": effect,
        "Impact": impact,
        "Gene": gene,
        "Rank": rank,
        "HGVS_C": hgvs_c,
        "HGVS_P": hgvs_p,
    }

    return name, genotypes, meta


def parse_vcf(vcf_fn):
    """
    Returns
    """
    vcf = pysam.VariantFile(vcf_fn)
    samples = list(vcf.header.samples)
    name_gts = {}
    name_gt_count = defaultdict(dict)
    name_meta = {}
    for record in vcf:
        name, genotypes, meta = extract_fields(record, samples)
        name_meta[name] = meta
        name_gts[name] = genotypes
        for gt in ["0/0", "0/1", "1/1", "NA"]:
            name_gt_count[name][gt] = genotypes.count(gt)
    df = pd.DataFrame(name_gts, index=samples)
    return df, name_gt_count, name_meta


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--sample", required=True)
    parser.add_argument("--vcf", required=True)
    args = parser.parse_args()

    gt_csv_fn = f"{args.sample}.GT.csv"
    count_json_fn = f"{args.sample}.scflvdj.count.json"
    meta_fn = f"{args.sample}.scflvdj.meta.json"
    stats_fn = f"{args.sample}.scflvdj.vcf.stats.json"

    df, name_gt_count, name_meta = parse_vcf(args.vcf)
    df.to_csv(gt_csv_fn)
    utils.write_json(name_gt_count, count_json_fn)
    utils.write_json(name_meta, meta_fn)
    utils.write_json({"Number of variant after filtering": len(name_meta)}, stats_fn)


if __name__ == "__main__":
    main()
