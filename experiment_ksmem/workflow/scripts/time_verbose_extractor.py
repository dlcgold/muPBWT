import sys
import csv

def main(argv):
    fieldnames = ['samples', 'sites', 'panel', 'n_query', 'tool', 'variant', 'user_time', 'sys_time', 'max_mem', 'wall_clock']
    res = {
        'tool': argv[0],
        'samples': int(argv[1]),
        'sites': int(argv[2]),
        'panel': argv[3],
        'n_query': int(argv[4]),
    }
    #print(sys.stdin)
    for line in sys.stdin:
        line = line[1:-1]
        tokens = line.split(sep=":")
        if tokens[0] == "Command being timed":
            # out.write(tokens[1].lstrip())
            # out.write("\n")
            if '/rlpbwt' in tokens[1]:
                if '-N' in tokens[1]:
                    res['variant'] = "naive"
                if '-B' in tokens[1]:
                    res['variant'] = "bitvectors"
                if '-P' in tokens[1] and '-r' not in tokens[1]:
                    if '-e' in tokens[1]:
                        res['variant'] = "panel extended"
                    else:
                        res['variant'] = "panel"
                if '-P' in tokens[1] and '-r' in tokens[1]:
                    if '-e' in tokens[1]:
                        res['variant'] = "panel extended raw"
                if '-S' in tokens[1] and '-r' not in tokens[1]:
                    if '-e' in tokens[1] and '-t' in tokens[1]:
                        res['variant'] = "slp_thr extended"
                    elif '-e' in tokens[1] and '-t' not in tokens[1]:
                        res['variant'] = "slp_no_thr extended"
                    elif '-e' not in tokens[1] and '-t' not in tokens[1]:
                        res['variant'] = "slp_no_thr "
                    else:
                        res['variant'] = "slp_thr"
                if '-S' in tokens[1] and '-r' in tokens[1]:
                    if '-e' in tokens[1] and '-t' in tokens[1]:
                        res['variant'] = "slp_thr_raw extended"
                    elif '-e' in tokens[1] and '-t' not in tokens[1]:
                        res['variant'] = "slp_no_thr_raw extended"
                    

            if '/pbwt' in tokens[1]:
                if '-matchNaive' in tokens[1]:
                    res['variant'] = "original"
                if '-matchIndexed' in tokens[1]:
                    res['variant'] = "indexed"
                if '-matchDynamic' in tokens[1]:
                    res['variant'] = "dynamic"

        if tokens[0] == "User time (seconds)":
            res['user_time'] = float(tokens[1].lstrip())
        if tokens[0] == "System time (seconds)":
            res['sys_time'] = float(tokens[1].lstrip())
        if tokens[0] == "Maximum resident set size (kbytes)":
            res['max_mem'] = int(tokens[1].lstrip())
        if tokens[0] == "Elapsed (wall clock) time (h":
            tot = 0.0
            for x in tokens[4:]:
                tot = tot*60 + float(x.lstrip())
                res['wall_clock'] = tot

    writer = csv.DictWriter(sys.stdout, fieldnames=fieldnames)
    writer.writeheader()
    writer.writerow(res)

if __name__ == "__main__":
    main(sys.argv[1:])
