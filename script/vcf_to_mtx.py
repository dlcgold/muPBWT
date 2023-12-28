#!/usr/bin/env python3

import sys


def main(argv):
    vcf = argv[0]
    haps = []
    with open(vcf, "r") as f:
        for line in f:
            if line[0] != "#":
                t = "".join(i for i in line.split()[9:] if i != "\n")

                haps.append(t.replace("|", ""))
    num_rows = len(haps)
    num_cols = len(haps[0])

    transpose_matrix = [[0] * num_rows for _ in range(num_cols)]

    # Riempie la matrice trasposta
    for i in range(num_rows):
        for j in range(num_cols):
            # print(int(haps[i][j]), end="")
            transpose_matrix[j][i] = int(haps[i][j])

    # print(transpose_matrix)
    for row in transpose_matrix:
        print("".join(map(str, row)))


if __name__ == "__main__":
    main(sys.argv[1:])
