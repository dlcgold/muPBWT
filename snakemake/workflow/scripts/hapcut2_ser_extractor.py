#!/usr/bin/env python3

import sys
import csv


def main(argv):
    fieldnames = [
        "mutation perc",
        "tool",
        "switch error rate",
        "mismatch rate",
        "flat rate",
        "phased count",
        "AN50",
        "N50",
        "num snps max block",
    ]
    res = {
        "mutation perc": argv[0],
        "tool": argv[1],
    }
    # print(sys.stdin)
    for line in sys.stdin:
        tokens = line.split(":")
        if len(tokens) != 2:
            continue
        tokens[1] = tokens[1].strip()
        if tokens[0] == "switch rate":
            res["switch error rate"] = float(tokens[1])
        elif tokens[0] == "mismatch rate":
            res["mismatch rate"] = float(tokens[1])
        elif tokens[0] == "flat rate":
            res["flat rate"] = float(tokens[1])
        elif tokens[0] == "phased count":
            res["phased count"] = int(tokens[1])
        elif tokens[0] == "AN50":
            res["AN50"] = float(tokens[1])
        elif tokens[0] == "N50":
            res["N50"] = float(tokens[1])
        elif tokens[0] == "num snps max blk":
            res["num snps max block"] = int(tokens[1])
        else:
            continue
    writer = csv.DictWriter(sys.stdout, fieldnames=fieldnames)
    writer.writeheader()
    writer.writerow(res)


if __name__ == "__main__":
    main(sys.argv[1:])
