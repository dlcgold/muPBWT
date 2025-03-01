#!/usr/bin/env python3

import sys
import random


def mutate_genotype(genotype):
    if genotype == "0/0":
        return random.choice(["1/1", "0/1"])
    elif genotype == "1/1":
        return random.choice(["0/0", "0/1"])
    elif genotype == "0/1":
        return random.choice(["0/0", "1/1"])
    return genotype


def main():

    percentage = int(sys.argv[2]) / 100

    with open(sys.argv[1], "r") as infile, open(sys.argv[3], "w") as outfile:
        for line in infile:
            if line.startswith("#"):
                outfile.write(line)
            else:
                fields = line.strip().split("\t")
                genotype = fields[9].split(":")[0]
                # print(fields)
                if random.random() < percentage:
                    fields[9] = mutate_genotype(genotype)
                ###print(fields)
                outfile.write("\t".join(fields) + "\n")


if __name__ == "__main__":
    main()
