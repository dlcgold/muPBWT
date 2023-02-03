import argparse
from pysam import VariantFile
import pysam 
def main(args):
    input = args.input[0]
    panel = None
    if args.panel is not None:
        panel = args.panel[0]
    queries = None
    if args.queries is not None:
        queries = args.queries[0]
    output = args.output[0]
    save = pysam.set_verbosity(0)
    panel_d = {}
    queries_d = {}
    if panel is not None:
        tmp = VariantFile(panel)
        tmp_s = list((tmp.header.samples))
        j = 0
        for i in range(0, len(tmp_s)*2, 2):
            panel_d[i] = f"{tmp_s[j]}_0"
            panel_d[i + 1] = f"{tmp_s[j]}_1"
            j += 1
    if queries is not None:
        tmp = VariantFile(queries)
        tmp_s = list((tmp.header.samples))
        j = 0
        for i in range(0, len(tmp_s)*2, 2):
            queries_d[i] = f"{tmp_s[j]}_0"
            queries_d[i + 1] = f"{tmp_s[j]}_1"
            j += 1
     
    
    with open(input, "r") as i:
        with open(output, "w") as o:
            lines = i.readlines()
            for line in lines:  
                tokens = line.split()
                if len(queries_d) != 0:
                    tokens[1] = queries_d[int(tokens[1])]
                if len(panel_d) != 0:
                    tokens[2] = panel_d[int(tokens[2])]
                o.write("\t".join(tokens))
                o.write("\n")
    pysam.set_verbosity(save)

if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("-i", "--input", nargs=1, help="SMEM file in Durbin's format")
    parser.add_argument("-p", "--panel", nargs=1, help="panel as VCF/BCF (optional)")
    parser.add_argument("-q", "--queries", nargs=1, help="queries as VCF/BCF (optional)")
    parser.add_argument("-o", "--output", nargs=1, help="output file")

    args = parser.parse_args()

    if args.input is None:
        print("No input file")
        exit(1)
    if args.output is None:
        print("No output file")
        exit(1)
    if args.panel is None and args.queries is None:
        print("Nothing to do: specify panel and/or queries")
        exit(1)
    main(args)