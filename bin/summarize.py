#!/usr/bin/env python

import argparse
import copy
import math
from collections import defaultdict

import pandas as pd
import pysam
import utils
from __init__ import ASSAY

MAX_CELL = 2 * 10**5


def get_umi_count(bcs, umis, cell_barcodes, sample):
    """
    Args:
        bcs: raw barcodes
        umis: umi count
        cell_barcodes: cell barcodes
    """
    a = [(umi, bc) for umi, bc in zip(umis, bcs) if umi > 0]
    a.sort(reverse=True)
    cell_barcodes = set(cell_barcodes)
    plot_data = {}
    n = len(a)
    first_noncell = n - 1
    for i, (umi, bc) in enumerate(a):
        if bc not in cell_barcodes:
            first_noncell = i
            break
    print(f"first non-cell barcode rank: {first_noncell}")
    last_cell = 0
    for i in range(min(n - 1, MAX_CELL), -1, -1):
        bc = a[i][1]
        if bc in cell_barcodes:
            last_cell = i
            break
    pure = sample + ".cells.pure" + f"({first_noncell}/{first_noncell}, 100%)"
    bg_cells = n - first_noncell
    bg = sample + ".cells.background" + f"(0/{bg_cells}, 0%)"
    plot_data[pure] = {}
    plot_data[bg] = {}
    for i in range(first_noncell + 1):  # for plot
        plot_data[pure][i + 1] = int(a[i][0])

    n_mix = last_cell - first_noncell + 1
    if n_mix != 0:
        n_total = len(cell_barcodes)
        n_mix_cell = n_total - first_noncell
        mix_rate = round(n_mix_cell / n_mix * 100, 2)
        mix = sample + ".cells.mix" + f"({n_mix_cell}/{n_mix}, {mix_rate}%)"
        plot_data[mix] = {}
        for i in range(first_noncell, last_cell + 1):
            plot_data[mix][i + 1] = int(a[i][0])

    for i in range(last_cell + 1, min(MAX_CELL, n), 10):
        plot_data[bg][i + 1] = int(a[i][0])
    # do not record every umi count
    for i in range(MAX_CELL, n, 1000):
        plot_data[bg][i + 1] = int(a[i][0])
    return plot_data


class Auto:
    """
    threshold = top {percentile}% cell count / coef
    count is usually UMI count.
    >>> array = [50] * 100 + [30] * 100 + [10] * 100 + [4] * 100
    >>> Auto(array, coef=10).run()
    5
    >>> Auto(array, percentile=70, coef=3).run()
    10
    >>> Auto(array, percentile=50, coef=10, expected_cell_num=100).run()
    5
    >>> Auto([1, 2, 20, 30, 40], expected_cell_num=4, percentile=50, coef=10).run()
    2
    """

    def __init__(self, array, percentile=99, coef=3, expected_cell_num=None, **kwargs):
        self.array = [x for x in array if x > 0]
        self.percentile = percentile
        self.coef = int(coef)
        self.expected_cell_num = int(expected_cell_num)
        self.kwargs = kwargs

    @staticmethod
    def calculate_percentile(array, percentile):
        size = len(array)
        return sorted(array)[int(math.ceil((size * percentile) / 100)) - 1]

    def run(self):
        array = self.array
        if not array:
            return 1

        if not self.expected_cell_num:
            expected_cell_num = len(array)
        else:
            expected_cell_num = self.expected_cell_num
            if expected_cell_num > len(array):
                print("Warning: expected_cell_num > len(array)")
                expected_cell_num = len(array)

        sorted_counts = sorted(array, reverse=True)
        count_cell_percentile = self.calculate_percentile(sorted_counts[:expected_cell_num], self.percentile)
        threshold = int(count_cell_percentile / self.coef)

        return threshold


def get_vdj_metric(df, chains, pairs):
    """
    Add vdj metrics.
    """
    data_dict = {}
    fl_pro_pair_df = pd.DataFrame(df[df["productive"]].barcode.value_counts())
    fl_pro_pair_df = fl_pro_pair_df[fl_pro_pair_df["count"] >= 2]
    cell_nums = len(set(df["barcode"]))

    data_dict = {"Cells With Productive V-J Spanning Pair": utils.get_frac(fl_pro_pair_df.shape[0] / cell_nums)}

    for pair in pairs:
        chain1, chain2 = pair.split("_")[0], pair.split("_")[1]
        cbs1 = set(df[(df["full_length"]) & (df["productive"]) & (df["chain"] == chain1)].barcode)
        cbs2 = set(df[(df["full_length"]) & (df["productive"]) & (df["chain"] == chain2)].barcode)
        paired_cbs = len(cbs1.intersection(cbs2))

        data_dict.update(
            {f"Cells With Productive V-J Spanning ({chain1}, {chain2}) Pair": utils.get_frac(paired_cbs / cell_nums)}
        )

    for chain in chains:
        value = len(set(df[df["chain"] == chain].barcode))
        data_dict.update({f"Cells With {chain} Contig": utils.get_frac(value / cell_nums)})

        value = len(set(df[(df["chain"] == chain) & (df["cdr3"] is not None)].barcode))
        data_dict.update({f"Cells With CDR3-annotated {chain} Contig": utils.get_frac(value / cell_nums)})

        value = len(set(df[(df["full_length"]) & (df["chain"] == chain)].barcode))
        data_dict.update({f"Cells With V-J Spanning {chain} Contig": utils.get_frac(value / cell_nums)})

        value = len(set(df[(df["full_length"]) & (df["productive"]) & (df["chain"] == chain)].barcode))
        data_dict.update({f"Cells With Productive {chain} Contig": utils.get_frac(value / cell_nums)})

    return data_dict


def target_cell_calling(df_umi_sum, expected_target_cell_num=3000, coef=5, percentile=85, umi_col="umis"):
    """
    Args:
        df_UMI_sum: A dataframe with columns highest umi's contig and UMI.

    Returns:
        target_contigs_id: list
    """
    umi_threshold = Auto(
        list(df_umi_sum[umi_col]), expected_cell_num=expected_target_cell_num, coef=coef, percentile=percentile
    ).run()

    # avoid change the original dataframe
    df_temp = df_umi_sum.copy()
    target_contigs = set(df_temp.loc[df_temp[umi_col] >= umi_threshold].contig_id)

    return target_contigs


def parse_contig_file(sample, barcode_report, annot_fa):
    """
    Generate all_contig_annotation file.
    Generate all_contig_fasta file.
    """

    df = pd.read_csv(barcode_report)

    df["productive"] = df["full_length"]
    contig_set = set(df.contig_id)

    # generate all contig fasta file
    # add length of each contig.
    len_dict = dict()
    all_fa = open(f"{sample}_all_contig.fasta", "w")
    with pysam.FastxFile(annot_fa) as fa:
        for read in fa:
            len_dict[read.name] = read.comment.split(" ")[0]
            if read.name in contig_set:
                sequence = read.sequence
                all_fa.write(">" + read.name + "\n" + sequence + "\n")
    all_fa.close()
    df["length"] = df["contig_id"].apply(len_dict.get)

    return df


def cell_calling(df, seqtype, trust_report, expected_target_cell_num, coef):
    """
    Common filtering based on CDR3:
    Filter nonfunctional CDR3(shown 'out_of_frame' in cdr3 report), or CDR3 sequences containing "N" in the nucleotide sequence.
    Keep CDR3aa start with C.
    Keep CDR3aa length >= 5.
    Keep no stop codon in CDR3aa.
    Filter low abundance contigs based on a umi cut-off.
    """
    df.sort_values(by="umis", ascending=False, inplace=True)
    if seqtype == "BCR":
        df_chain_heavy = df[df["chain"] == "IGH"]
        df_chain_light = df[(df["chain"] == "IGK") | (df["chain"] == "IGL")]
    else:
        df_chain_heavy = df[df["chain"] == "TRA"]
        df_chain_light = df[df["chain"] == "TRB"]
    df_chain_heavy = df_chain_heavy.drop_duplicates(["barcode"])
    df_chain_light = df_chain_light.drop_duplicates(["barcode"])
    df_for_clono = pd.concat([df_chain_heavy, df_chain_light], ignore_index=True)

    # Common filtering
    trust_report = pd.read_csv(trust_report, sep="\t")
    correct_cdr3 = set(df_for_clono.cdr3).intersection(set(trust_report.CDR3aa))
    correct_cdr3 = [i for i in correct_cdr3 if i.startswith("C")]
    correct_cdr3 = [i for i in correct_cdr3 if len(i) >= 5]
    correct_cdr3 = [i for i in correct_cdr3 if "UAG" or "UAA" or "UGA" not in i]
    df_for_clono = df_for_clono[df_for_clono["cdr3"].isin(correct_cdr3)]

    # Filter low abundance contigs based on a umi cut-off
    if seqtype == "BCR":
        df_chain_heavy = df_for_clono[df_for_clono["chain"] == "IGH"]
        df_chain_light = df_for_clono[(df_for_clono["chain"] == "IGK") | (df_for_clono["chain"] == "IGL")]
    else:
        df_chain_heavy = df_for_clono[df_for_clono["chain"] == "TRA"]
        df_chain_light = df_for_clono[df_for_clono["chain"] == "TRB"]

    filtered_congtigs_id = set()
    for _df in [df_chain_heavy, df_chain_light]:
        target_contigs = target_cell_calling(_df, expected_target_cell_num=expected_target_cell_num, coef=coef)
        filtered_congtigs_id = filtered_congtigs_id | target_contigs

    df_for_clono = df_for_clono[df_for_clono.contig_id.isin(filtered_congtigs_id)]

    df_for_clono_pro = df_for_clono[df_for_clono["productive"]]
    cell_barcodes, filtered_contig = set(df_for_clono_pro["barcode"]), set(df_for_clono_pro["contig_id"])

    return df_for_clono, cell_barcodes, filtered_contig


def filter_fasta(sample, cell_barcodes):
    """Filter all contig fasta file by barcodes which are identified to be cell.

    :param cell_barcodes: all barcodes identified to be cell.
    """
    all_contig_fasta = f"{sample}_all_contig.fasta"
    filter_contig_fasta = f"{sample}_filtered_contig.fasta"

    filter_contig_fasta = open(filter_contig_fasta, "w")
    with pysam.FastxFile(all_contig_fasta) as fa:
        for read in fa:
            name = read.name
            cb = "_".join(name.split("_")[:-1])  # remove the suffix num in "bc1_bc2_bc3_num"
            sequence = read.sequence
            if cb in cell_barcodes:
                filter_contig_fasta.write(">" + name + "\n" + sequence + "\n")
    filter_contig_fasta.close()


def parse_clonotypes(sample, df, df_for_clono, cell_barcodes, filtered_contig):
    """Parse clonotypes from CDR3 and manually add clonotype id for each contig.

    :param df: original contig file.
    :param df_for_clono: contig info after filter.
    :param cell_barcodes: all barcodes identified to be cell.
    :return df_filter_contig: filtered contigs by cell barcodes.
    """
    df_for_clono_pro = df_for_clono[df_for_clono["productive"]].copy()
    df_for_clono_pro["chain_cdr3aa"] = df_for_clono_pro.loc[:, ["chain", "cdr3"]].apply(":".join, axis=1)
    df_for_clono_pro["chain_cdr3nt"] = df_for_clono_pro.loc[:, ["chain", "cdr3_nt"]].apply(":".join, axis=1)

    cbs = set(df_for_clono_pro["barcode"])
    clonotypes = open(f"{sample}_clonotypes.csv", "w")
    clonotypes.write("barcode\tcdr3s_aa\tcdr3s_nt\n")
    for cb in cbs:
        temp = df_for_clono_pro[df_for_clono_pro["barcode"] == cb]
        temp = temp.sort_values(by="chain", ascending=True)
        aa_chain = ";".join(list(temp["chain_cdr3aa"]))
        nt_chain = ";".join(list(temp["chain_cdr3nt"]))
        clonotypes.write(f"{cb}\t{aa_chain}\t{nt_chain}\n")
    clonotypes.close()

    df_clonotypes = pd.read_csv(f"{sample}_clonotypes.csv", sep="\t", index_col=None)
    contig_with_clonotype = copy.deepcopy(df_clonotypes)
    df_dict = df_clonotypes[["cdr3s_nt", "cdr3s_aa"]].set_index("cdr3s_nt").to_dict(orient="dict")["cdr3s_aa"]
    df_clonotypes = df_clonotypes.groupby("cdr3s_nt", as_index=False).agg({"barcode": "count"})
    df_clonotypes.rename(columns={"barcode": "frequency"}, inplace=True)
    sum_f = df_clonotypes["frequency"].sum()

    df_clonotypes["proportion"] = df_clonotypes["frequency"].apply(lambda x: x / sum_f)
    df_clonotypes.sort_values(by="frequency", ascending=False, inplace=True)
    df_clonotypes["clonotype_id"] = [f"clonotype{i}" for i in range(1, df_clonotypes.shape[0] + 1)]
    df_clonotypes["cdr3s_aa"] = df_clonotypes["cdr3s_nt"].apply(lambda x: df_dict[x])
    df_clonotypes = df_clonotypes.reindex(columns=["clonotype_id", "frequency", "proportion", "cdr3s_aa", "cdr3s_nt"])
    df_clonotypes.to_csv(f"{sample}_clonotypes.csv", sep=",", index=False)
    used_for_merge = df_clonotypes[["cdr3s_nt", "clonotype_id"]]

    df_merge = pd.merge(used_for_merge, contig_with_clonotype, on="cdr3s_nt", how="outer")
    df_merge = df_merge[["barcode", "clonotype_id"]]
    df_all_contig = pd.merge(df_merge, df, on="barcode", how="outer")
    df_all_contig.fillna("", inplace=True)
    df_all_contig = df_all_contig[
        [
            "barcode",
            "is_cell",
            "contig_id",
            "high_confidence",
            "length",
            "chain",
            "v_gene",
            "d_gene",
            "j_gene",
            "c_gene",
            "full_length",
            "productive",
            "cdr3",
            "cdr3_nt",
            "reads",
            "umis",
            "clonotype_id",
        ]
    ]
    df_filter_contig = df_all_contig[df_all_contig["barcode"].isin(cell_barcodes)]
    for _df in [df_all_contig, df_filter_contig]:
        _df.loc[~_df.contig_id.isin(filtered_contig), "clonotype_id"] = ""

    df_all_contig.to_csv(f"{sample}_all_contig.csv", sep=",", index=False)
    df_filter_contig.to_csv(f"{sample}_filtered_contig.csv", sep=",", index=False)


def gen_summary(sample, fq2, assembled_reads, seqtype, df_for_clono):
    """Generate metrics in html"""
    df_for_clono_pro = df_for_clono[df_for_clono["productive"]]
    cell_barcodes = set(df_for_clono_pro["barcode"])
    total_cells = len(cell_barcodes)

    read_count = 0
    read_count_all = 0
    umi_dict = defaultdict(set)
    umi_count = defaultdict()
    with pysam.FastxFile(fq2) as fq:
        for read in fq:
            read_count_all += 1
            cb = read.name.split(":")[0]
            umi = read.name.split(":")[1]
            umi_dict[cb].add(umi)
            if cb in cell_barcodes:
                read_count += 1

    for cb in umi_dict:
        umi_count[cb] = len(umi_dict[cb])
    df_umi = pd.DataFrame.from_dict(umi_count, orient="index", columns=["UMI"])
    df_umi["barcode"] = df_umi.index
    df_umi = df_umi.reset_index(drop=True)
    df_umi = df_umi.reindex(columns=["barcode", "UMI"])
    df_umi = df_umi.sort_values(by="UMI", ascending=False)
    df_umi["mark"] = df_umi["barcode"].apply(lambda x: "CB" if x in cell_barcodes else "UB")
    df_umi.to_csv(f"{sample}.count.txt", sep="\t", index=False)
    # self.add_data(chart=get_plot_elements.plot_barcode_rank(f'{self.outdir}/count.txt'))

    used_read = 0
    with pysam.FastxFile(assembled_reads) as fa:
        for read in fa:
            bc = read.name.split(":")[0]
            if bc in cell_barcodes:
                used_read += 1

    if seqtype == "TCR":
        chains = ["TRA", "TRB"]
    else:
        chains = ["IGH", "IGL", "IGK"]

    data_dict = {
        "Estimated Number of Cells": total_cells,
        "Mean Read Pairs per Cell": int(read_count / total_cells),
        "Mean Used Read Pairs per Cell": int(used_read / total_cells),
        "Fraction of Reads in Cells": utils.get_frac(used_read / read_count_all),
    }

    for c in chains:
        temp_df = df_for_clono_pro[df_for_clono_pro["chain"] == c]

        try:
            median_umi = int(temp_df["umis"].median())
        except ValueError:
            # ValueError: cannot convert float NaN to integer
            median_umi = 0

        data_dict.update({f"Median {c} UMIs per Cell": median_umi})
    fn = f"{args.sample}.{ASSAY}.summarize.stats.json"
    utils.write_json(data_dict, fn)

    # UMI count
    plot_data = get_umi_count(df_umi["barcode"], df_umi["UMI"], cell_barcodes, args.sample)
    fn = f"{args.sample}.{ASSAY}.umi_count.json"
    utils.write_json(plot_data, fn)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--sample", required=True)
    parser.add_argument("--seqtype", required=True)
    parser.add_argument("--coef", required=True)
    parser.add_argument("--expected_target_cell_num", required=True)
    parser.add_argument("--fq2", required=True)
    parser.add_argument("--assembled_reads", required=True)
    parser.add_argument("--filter_report_tsv", required=True)
    parser.add_argument("--annot_fa", required=True)
    parser.add_argument("--barcode_report", required=True)
    args = parser.parse_args()

    if args.seqtype == "BCR":
        chains = ["IGH", "IGL", "IGK"]
        paired_groups = ["IGK_IGH", "IGL_IGH"]
    else:
        chains = ["TRA", "TRB"]
        paired_groups = ["TRA_TRB"]

    # cell calling
    original_df = parse_contig_file(args.sample, args.barcode_report, args.annot_fa)
    df_for_clono, cell_barcodes, filtered_contig = cell_calling(
        original_df, args.seqtype, args.filter_report_tsv, args.expected_target_cell_num, args.coef
    )
    filter_fasta(args.sample, cell_barcodes)
    parse_clonotypes(args.sample, original_df, df_for_clono, cell_barcodes, filtered_contig)
    gen_summary(args.sample, args.fq2, args.assembled_reads, args.seqtype, df_for_clono)

    # vdj metrics
    df = pd.read_csv(f"{args.sample}_filtered_contig.csv")
    data_dict = get_vdj_metric(df, chains, paired_groups)
    fn = f"{args.sample}.{ASSAY}.annotation.stats.json"
    utils.write_json(data_dict, fn)
