#!/usr/bin/env python

import argparse
import csv

import utils


def parse_summary(f):
    parsed_data = {}
    reader = csv.reader(f)
    for row in reader:
        parsed_data[row[0]] = row[1]
    return parsed_data


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Starsolo summary")
    parser.add_argument("--summary", help="starsolo summary file")
    parser.add_argument("--sample", help="sample name")
    args = parser.parse_args()

    # summary
    parsed_data = parse_summary(open(args.summary))
    data = {}
    data["Raw Reads"] = parsed_data["Number of Reads"]
    data["Valid Reads"] = parsed_data["Reads With Valid Barcodes"]
    summary_file = args.sample + ".scflvdj.read.stats.json"
    utils.write_json(data, summary_file)
