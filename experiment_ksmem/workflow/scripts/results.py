from turtle import color
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.ticker import StrMethodFormatter, NullFormatter
import numpy as np
import os
import statistics
import sys

def get_size(file_path, unit='bytes', r = 2):
    file_size = os.path.getsize(file_path)
    exponents_map = {'bytes': 0, 'kb': 1, 'mb': 2, 'gb': 3}
    if unit not in exponents_map:
        raise ValueError("Must select from \
        ['bytes', 'kb', 'mb', 'gb']")
    else:
        size = file_size / 1024 ** exponents_map[unit]
        return round(size, r)


def main(argv):
    if os.path.exists('~/.matplotlib/stylelib/nord-light.mplstyle'):
        plt.style.use('~/.matplotlib/stylelib/nord-light.mplstyle')

    PANELS = ['1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15', '16', '17', '18', '19', '20', '21', '22']
    PANELS.reverse()
    kb_to_gb = 9.5367431640625e-7
    b_to_gb = 9.313225746154785e-10
    m_to_gb = 0.00097656
    
    input_dir = argv[0]
    output_dir = argv[1]
    results_dir = argv[2]
    
    if input_dir[-1] == '/':
        input_dir = input_dir[:-1]
    if output_dir[-1] == '/':
        output_dir = output_dir[:-1]
    if results_dir[-1] == '/':
        results_dir = results_dir[:-1]
        
    data_dir = f"{results_dir}/data"
    plot_dir = f"{results_dir}/plots"
    table_dir = f"{results_dir}/tables"

    if not os.path.exists(plot_dir):
        os.makedirs(plot_dir)
    if not os.path.exists(table_dir):
        os.makedirs(table_dir)

    bcf_size = []
    mupbwt_size = []
    syllable_size = []

    bgt_size = []
    
    for e in PANELS:
        bcf_size.append(get_size(f"{input_dir}/{e}/panel_{e}.bcf",'gb'))
        mupbwt_size.append(get_size(f"{output_dir}/mupbwt/{e}/{e}.ser",'gb'))
        syllable_size.append(get_size(f"{output_dir}/syllable/{e}/chr_{e}-save0.txt",'gb'))
        bgt_size.append(get_size(f"{output_dir}/bgt/{e}/{e}.bgt.bcf",'gb') +
                        get_size(f"{output_dir}/bgt/{e}/{e}.bgt.bcf.csi",'gb') +
                        get_size(f"{output_dir}/bgt/{e}/{e}.bgt.pbf",'gb') +
                        get_size(f"{output_dir}/bgt/{e}/{e}.bgt.spl",'gb'))

    df_make = pd.read_csv(f"{data_dir}/make-time.csv")
    df_exe = pd.read_csv(f"{data_dir}/exe-time.csv")

    df_make_pbwt = df_make[df_make["tool"] == "pbwt"]
    df_make_pbwt = df_make_pbwt[::-1].reset_index(drop=True)
    
    df_make_mupbwt = df_make[df_make["tool"] == "mupbwt"]
    df_make_mupbwt = df_make_mupbwt[::-1].reset_index(drop=True)

    df_make_syllable = df_make[df_make["tool"] == "syllable"]
    df_make_syllable = df_make_syllable[::-1].reset_index(drop=True)

    df_make_bgt = df_make[df_make["tool"] == "bgt"]
    df_make_bgt = df_make_bgt[::-1].reset_index(drop=True)

    make_mem_pbwt = []
    make_mem_mupbwt = []
    make_mem_syllable = []
    make_mem_bgt = []
    make_time_pbwt = []
    make_time_mupbwt = []
    make_time_syllable = []
    make_time_bgt = []
    
    for e in df_make_pbwt["max_mem"]:
        make_mem_pbwt.append(round(e * kb_to_gb, 2))
        
    for e in df_make_mupbwt["max_mem"]:
        make_mem_mupbwt.append(round(e * kb_to_gb, 2))
        
    for e in df_make_syllable["max_mem"]:
        make_mem_syllable.append(round(e * kb_to_gb, 2))

    for e in df_make_bgt["max_mem"]:
        make_mem_bgt.append(round(e * kb_to_gb, 2))
        
    for e in df_make_pbwt["wall_clock"]:
        make_time_pbwt.append(e)

    for e in df_make_mupbwt["wall_clock"]:
        make_time_mupbwt.append(e)

    for e in df_make_syllable["wall_clock"]:
        make_time_syllable.append(e)

    for e in df_make_bgt["wall_clock"]:
        make_time_bgt.append(e)


    df_exe_pbwt = df_exe[df_exe["tool"] == "pbwt"]
    df_exe_pbwt = df_exe_pbwt[df_exe_pbwt["variant"] == "indexed"]
    df_exe_pbwt = df_exe_pbwt[::-1].reset_index(drop=True)

    df_exe_mupbwt = df_exe[df_exe["tool"] == "mupbwt"]
    df_exe_mupbwt = df_exe_mupbwt[::-1].reset_index(drop=True)

    
    exe_mem_pbwt = []
    exe_mem_mupbwt = []
    exe_time_pbwt = []
    exe_time_mupbwt = []
    
    for e in df_exe_pbwt["max_mem"]:
        exe_mem_pbwt.append(round(e * kb_to_gb, 2))
        
    for e in df_exe_mupbwt["max_mem"]:
        exe_mem_mupbwt.append(round(e * kb_to_gb, 2))
        
    for e in df_exe_pbwt["wall_clock"]:
        exe_time_pbwt.append(e)

    for e in df_exe_mupbwt["wall_clock"]:
        exe_time_mupbwt.append(e)


    load_d  = {}
    chr_sites_d = {}
    chr_runt_d = {}
    chr_run_d = {}
    map_mem_d = {}
    thr_mem_d = {}
    sam_mem_d = {} 
    phi_mem_d = {} 
    for e in PANELS:
        with open(f"{output_dir}/mupbwt/{e}/stat.txt", "r") as f:
            l = f.readlines()
            load_d[e] = l[0].split()[2]
            chr_sites_d[e] = l[4].split()[2]
            chr_runt_d[e] = l[6].split()[2]
            chr_run_d[e] = l[7].split()[2]
            
            map_mem_d[e] = l[13].split()[2]
            thr_mem_d[e] = l[10].split()[1]
            sam_mem_d[e] = l[12].split()[1]
            phi_mem_d[e] = l[16].split()[6]
        

    load = []
    chr_sites = []
    chr_runt = []
    chr_run = []
    map_mem = []
    thr_mem = []
    sam_mem = [] 
    phi_mem = [] 
    for e in PANELS:
        load.append(float(load_d[e]))
        chr_sites.append(int(chr_sites_d[e]))
        chr_runt.append(int(chr_runt_d[e]))
        chr_run.append(int(chr_run_d[e]))
        thr_mem.append(float(thr_mem_d[e]))
        sam_mem.append(float(sam_mem_d[e]))
        map_mem.append(float(map_mem_d[e]) - float(thr_mem_d[e]) - float(sam_mem_d[e]))
        phi_mem.append(float(phi_mem_d[e]))

    with open(f"{table_dir}/mupbwt_build.tex", "w") as f:
        f.write("\\begin{tabular}{c|r|r|r|r|r|r}\n\tChr & Samples & Sites & Runs & BCF & $\mu$-PBWT & Constuction time & Construction max mem\\\\ \n\t\\hline\n")
        i = 0
        tot_s = 0
        tot_r = 0
        tot_b = 0
        tot_se = 0
        tot_m = 0
        for e in PANELS[::-1]:
            f.write(f"\t{e} & 2454 & {chr_sites_d[e]} & {chr_run_d[e]} & {round(bcf_size[len(bcf_size)-i-1], 2)} & {round(mupbwt_size[len(mupbwt_size)-i-1], 2)} & 00:{int(round(make_time_mupbwt[len(make_time_mupbwt)-i-1]/60,  0))} & {round(make_mem_mupbwt[len(make_mem_mupbwt)-i-1], 2)} \\\\\n")
        
            tot_s += int(chr_sites_d[e])
            tot_r += int(chr_run_d[e])
            tot_b += round(bcf_size[i], 2)
            tot_se += round(mupbwt_size[i], 2)
            tot_m += make_mem_mupbwt[len(make_mem_mupbwt)-i-1]
            i = i+1
        f.write(f"\tTotal & 2454 & {tot_s} & {tot_r/i} & {round(tot_b, 2)} & {round(tot_se, 2)} & - & {tot_m} \n \\end{{tabular}}")
        f.write("\n")


    with open(f"{table_dir}/mupbwt_compare.tex", "w") as f:
        f.write("\\begin{tabular}{c|r|r|r|r|r}\n\tChr & Samples & Sites & Runs & BCF & $\mu$-PBWT & Syllable-PBWT\\\\ \n\t\\hline\n")
        i = 0
        tot_s = 0
        tot_r = 0
        tot_b = 0
        tot_se = 0
        tot_m = 0
        for e in PANELS[::-1]:
            f.write(f"\t{e} & 2454 & {chr_sites_d[e]} & {chr_run_d[e]} & {round(bcf_size[len(bcf_size)-i-1], 2)} & {round(mupbwt_size[len(mupbwt_size)-i-1], 2)} & {round(syllable_size[len(syllable_size)-i-1], 2)}  \\\\\n")
        
            tot_s += int(chr_sites_d[e])
            tot_r += int(chr_run_d[e])
            tot_b += round(bcf_size[i], 2)
            tot_se += round(mupbwt_size[i], 2)
            tot_m += round(syllable_size[i], 2)
            i = i+1
        f.write(f"\tTotal & 2454 & {tot_s} & {int(round(tot_r/i, 0))} & {round(tot_b, 2)} & {round(tot_se, 2)} & {tot_m} \n \\end{{tabular}}")
        f.write("\n")

    with open(f"{table_dir}/mupbwt_compare_bgt.tex", "w") as f:
        f.write("\\begin{tabular}{c|r|r|r|r|r|r}\n\tChr & Samples & Sites & Runs & BCF & $\mu$-PBWT & Syllable-PBWT & BGT\\\\ \n\t\\hline\n")
        i = 0
        tot_s = 0
        tot_r = 0
        tot_b = 0
        tot_se = 0
        tot_m = 0
        tot_h = 0
        for e in PANELS[::-1]:
           f.write(f"\t{e} & 2454 & {chr_sites_d[e]} & {chr_run_d[e]} & {round(bcf_size[len(bcf_size)-i-1], 2)} & {round(mupbwt_size[len(mupbwt_size)-i-1], 2)} & {round(syllable_size[len(syllable_size)-i-1], 2)} & {round(bgt_size[len(bgt_size)-i-1], 2)}  \\\\\n")
            tot_s += int(chr_sites_d[e])
            tot_r += int(chr_run_d[e])
            tot_b += round(bcf_size[i], 2)
            tot_se += round(mupbwt_size[i], 2)
            tot_m += round(syllable_size[i], 2)
            tot_h += round(bgt_size[i], 2)
            i = i+1
        f.write(f"\tTotal & 2454 & {tot_s} & {int(round(tot_r/i, 0))} & {round(tot_b, 2)} & {round(tot_se, 2)} & {tot_m} & {tot_h} \n \\end{{tabular}}")
        f.write("\n")


        p = []
        for e in PANELS:
            p.append(int(e))
        pi = p[::-1]


        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10,5))


        L2 = ax1.plot(pi, make_mem_mupbwt, marker="d", color = '#d08770')[0]
        
        L3 = ax1.plot(pi, exe_mem_mupbwt, marker="v", color= '#8fbcbb')[0]

        L5 = ax1.plot(pi, exe_mem_pbwt, marker="o", color = "#81a1c1")[0]

        L7 = ax1.plot(pi, make_mem_syllable, marker=">",  color = '#4c566a')[0]

        ax1.xaxis.set_major_formatter(StrMethodFormatter('{x:.0f}'))
        ax1.xaxis.set_minor_formatter(NullFormatter())
        ax1.tick_params(axis='x', labelrotation = 65)
        ax1.set_title(f"a)1KGP memory usage", fontweight="bold")
        ax1.set_xticks(pi)
        ax1.set_xticklabels(p)
        
        ax1.set_xlabel("Chromosome", fontweight="bold")
        ax1.set_ylabel("Memory (gigabytes)", fontweight="bold")
        ax1.semilogy()



        L2 = ax2.plot(pi, make_time_mupbwt, marker="d", color = '#d08770')[0]
        
        L3 = ax2.plot(pi, exe_time_mupbwt, marker="v", color = '#8fbcbb')[0]

        L5 = ax2.plot(pi, exe_time_pbwt, marker="o", color = '#81a1c1')[0]

        L7 = ax2.plot(pi, make_time_syllable, marker=">",  color = '#4c566a')[0]


        ax2.xaxis.set_major_formatter(StrMethodFormatter('{x:.0f}'))
        ax2.xaxis.set_minor_formatter(NullFormatter())
        ax2.tick_params(axis='x', labelrotation = 65)
        ax2.set_title(f"b) 1KGP time", fontweight="bold")
        ax2.set_xticks(pi)
        ax2.set_xticklabels(p)
        ax2.set_xlabel("Chromosome", fontweight="bold")
        ax2.set_ylabel("Time (seconds)", fontweight="bold")
        ax2.semilogy()
    
        plt.tight_layout()

        line_labels = ["μ-PBWT building", "μ-PBWT quering", "PBWT", "Syllable-PBWT building"]
        fig.legend(handles= [ L2, L3, L5, L7],     # The line objects
                   labels=line_labels,   # The labels for each line
                   loc="lower center",   # Position of legend
                   #borderaxespad=0.1,    # Small spacing around legend box
                   bbox_to_anchor=(0.5, -0.10),
                   ncol=5
                   #title="Legend Title"  # Title for the legend
                   )
        plt.savefig(f"{plot_dir}/1kgb.pdf", dpi=500, bbox_inches='tight')



        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10,5))


        L2 = ax1.plot(pi, make_mem_mupbwt, marker="d", color = '#d08770')[0]
        
        L3 = ax1.plot(pi, exe_mem_mupbwt, marker="v", color= '#8fbcbb')[0]

        L7 = ax1.plot(pi, make_mem_syllable, marker=">",  color = '#4c566a')[0]

        ax1.xaxis.set_major_formatter(StrMethodFormatter('{x:.0f}'))
        ax1.xaxis.set_minor_formatter(NullFormatter())
        ax1.tick_params(axis='x', labelrotation = 65)
        ax1.set_title(f"a)1KGP memory usage", fontweight="bold")
        ax1.set_xticks(pi)
        ax1.set_xticklabels(p)
        
        ax1.set_xlabel("Chromosome", fontweight="bold")
        ax1.set_ylabel("Memory (gigabytes)", fontweight="bold")
        ax1.semilogy()



        L2 = ax2.plot(pi, make_time_mupbwt, marker="d", color = '#d08770')[0]
        
        L3 = ax2.plot(pi, exe_time_mupbwt, marker="v", color = '#8fbcbb')[0]

        L5 = ax2.plot(pi, exe_time_pbwt, marker="o", color = '#81a1c1')[0]

        L7 = ax2.plot(pi, make_time_syllable, marker=">",  color = '#4c566a')[0]


        ax2.xaxis.set_major_formatter(StrMethodFormatter('{x:.0f}'))
        ax2.xaxis.set_minor_formatter(NullFormatter())
        ax2.tick_params(axis='x', labelrotation = 65)
        ax2.set_title(f"b) 1KGP time", fontweight="bold")
        ax2.set_xticks(pi)
        ax2.set_xticklabels(p)
        ax2.set_xlabel("Chromosome", fontweight="bold")
        ax2.set_ylabel("Time (seconds)", fontweight="bold")
        #ax2.semilogy()
    
        plt.tight_layout()

        line_labels = ["μ-PBWT building", "μ-PBWT quering", "PBWT", "Syllable-PBWT building"]
        fig.legend(handles= [ L2, L3, L5, L7],     # The line objects
                   labels=line_labels,   # The labels for each line
                   loc="lower center",   # Position of legend
                   #borderaxespad=0.1,    # Small spacing around legend box
                   bbox_to_anchor=(0.5, -0.10),
                   ncol=5
                   #title="Legend Title"  # Title for the legend
                   )
        plt.savefig(f"{plot_dir}/1kgb_nodu.pdf", dpi=500, bbox_inches='tight')


        fig, (ax2) = plt.subplots(1, 1, figsize=(5,5))

        L1 = ax2.plot(pi, mupbwt_size , marker="d", color = '#b48ead')[0]
        
        L2 = ax2.plot(pi, syllable_size, marker="s", color = '#d08770')[0]

        L3 = ax2.plot(pi, bcf_size , marker="p", color = '#8fbcbb')[0]

        ax2.xaxis.set_major_formatter(StrMethodFormatter('{x:.0f}'))
        ax2.xaxis.set_minor_formatter(NullFormatter())
        ax2.tick_params(axis='x', labelrotation = 65)
        ax2.set_title(f"a) 1KGP file sizes", fontweight="bold")
        ax2.set_xticks(pi)
        ax2.set_xticklabels(p)
        ax2.set_xlabel("Chromosome", fontweight="bold")
        ax2.set_ylabel("Memory (gigabytes)", fontweight="bold")
        ax2.semilogy()


        plt.tight_layout()

        line_labels = ["μ-PBWT serialization", "Syllable-PBWT serialization", "BCF file"]
        fig.legend(handles= [L1, L2, L3],     # The line objects
                   labels=line_labels,   # The labels for each line
                   loc="lower center",   # Position of legend
                   #borderaxespad=0.1,    # Small spacing around legend box
                   bbox_to_anchor=(0.5, -0.05),
                   ncol=3
                   #title="Legend Title"  # Title for the legend
                   )
        plt.savefig(f"{plot_dir}/files_size.pdf", dpi=500, bbox_inches='tight')


        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10,5))


        L2 = ax1.plot(pi, make_mem_mupbwt, marker="d", color = '#d08770')[0]
        
        L3 = ax1.plot(pi, exe_mem_mupbwt, marker="v", color= '#8fbcbb')[0]

        L5 = ax1.plot(pi, exe_mem_pbwt, marker="o", color = "#81a1c1")[0]

        L7 = ax1.plot(pi, make_mem_syllable, marker=">",  color = '#4c566a')[0]

        L8 = ax1.plot(pi, make_mem_bgt, marker="<",  color = '#a3be8c')[0]

        ax1.xaxis.set_major_formatter(StrMethodFormatter('{x:.0f}'))
        ax1.xaxis.set_minor_formatter(NullFormatter())
        ax1.tick_params(axis='x', labelrotation = 65)
        ax1.set_title(f"a)1KGP memory usage", fontweight="bold")
        ax1.set_xticks(pi)
        ax1.set_xticklabels(p)
        
        ax1.set_xlabel("Chromosome", fontweight="bold")
        ax1.set_ylabel("Memory (gigabytes)", fontweight="bold")
        ax1.semilogy()
        ax1.yaxis.set_major_formatter(ticker.FuncFormatter(myLogFormat))


        L2 = ax2.plot(pi, make_time_mupbwt, marker="d", color = '#d08770')[0]
        
        L3 = ax2.plot(pi, exe_time_mupbwt, marker="v", color = '#8fbcbb')[0]

        L5 = ax2.plot(pi, exe_time_pbwt, marker="o", color = '#81a1c1')[0]

        L7 = ax2.plot(pi, make_time_syllable, marker=">",  color = '#4c566a')[0]

        L8 = ax2.plot(pi, make_time_bgt, marker="<",  color = '#a3be8c')[0]


        ax2.xaxis.set_major_formatter(StrMethodFormatter('{x:.0f}'))
        ax2.xaxis.set_minor_formatter(NullFormatter())
        ax2.tick_params(axis='x', labelrotation = 65)
        ax2.set_title(f"b) 1KGP time", fontweight="bold")
        ax2.set_xticks(pi)
        ax2.set_xticklabels(p)
        ax2.set_xlabel("Chromosome", fontweight="bold")
        ax2.set_ylabel("Time (seconds)", fontweight="bold")
        #ax2.semilogy()
    
        plt.tight_layout()

        line_labels = ["μ-PBWT building", "μ-PBWT quering", "PBWT", "Syllable-PBWT building",  "BGT building"]
        fig.legend(handles= [ L2, L3, L5, L7, L8],     # The line objects
                   labels=line_labels,   # The labels for each line
                   loc="lower center",   # Position of legend
                   #borderaxespad=0.1,    # Small spacing around legend box
                   bbox_to_anchor=(0.5, -0.10),
                   ncol=6
                   #title="Legend Title"  # Title for the legend
                   )
        plt.savefig(f"{plot_dir}/1kgb_bgt.pdf", dpi=500, bbox_inches='tight')



        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10,5))


        L2 = ax1.plot(pi, make_mem_mupbwt, marker="d", color = '#d08770')[0]
        
        L3 = ax1.plot(pi, exe_mem_mupbwt, marker="v", color= '#8fbcbb')[0]

        L7 = ax1.plot(pi, make_mem_syllable, marker=">",  color = '#4c566a')[0]

        L8 = ax1.plot(pi, make_mem_bgt, marker="<",  color = '#a3be8c')[0]

        ax1.xaxis.set_major_formatter(StrMethodFormatter('{x:.0f}'))
        ax1.xaxis.set_minor_formatter(NullFormatter())
        ax1.tick_params(axis='x', labelrotation = 65)
        ax1.set_title(f"a)1KGP memory usage", fontweight="bold")
        ax1.set_xticks(pi)
        ax1.set_xticklabels(p)
        
        ax1.set_xlabel("Chromosome", fontweight="bold")
        ax1.set_ylabel("Memory (gigabytes)", fontweight="bold")
        #ax1.semilogy()



        L2 = ax2.plot(pi, make_time_mupbwt, marker="d", color = '#d08770')[0]
        
        L3 = ax2.plot(pi, exe_time_mupbwt, marker="v", color = '#8fbcbb')[0]

        L5 = ax2.plot(pi, exe_time_pbwt, marker="o", color = '#81a1c1')[0]

        L7 = ax2.plot(pi, make_time_syllable, marker=">",  color = '#4c566a')[0]

        L8 = ax2.plot(pi, make_time_bgt, marker="<",  color = '#a3be8c')[0]

        ax2.xaxis.set_major_formatter(StrMethodFormatter('{x:.0f}'))
        ax2.xaxis.set_minor_formatter(NullFormatter())
        ax2.tick_params(axis='x', labelrotation = 65)
        ax2.set_title(f"b) 1KGP time", fontweight="bold")
        ax2.set_xticks(pi)
        ax2.set_xticklabels(p)
        ax2.set_xlabel("Chromosome", fontweight="bold")
        ax2.set_ylabel("Time (seconds)", fontweight="bold")
        #ax2.semilogy()
    
        plt.tight_layout()

        line_labels = ["μ-PBWT building", "μ-PBWT quering", "PBWT", "Syllable-PBWT building", "BGT building"]
        fig.legend(handles= [ L2, L3, L5, L7, L8],     # The line objects
                   labels=line_labels,   # The labels for each line
                   loc="lower center",   # Position of legend
                   #borderaxespad=0.1,    # Small spacing around legend box
                   bbox_to_anchor=(0.5, -0.10),
                   ncol=6
                   #title="Legend Title"  # Title for the legend
                   )
        plt.savefig(f"{plot_dir}/1kgb_nodu_bgt.pdf", dpi=500, bbox_inches='tight')


        fig, (ax2) = plt.subplots(1, 1, figsize=(5,5))

        L1 = ax2.plot(pi, mupbwt_size , marker="d", color = '#b48ead')[0]
        
        L2 = ax2.plot(pi, syllable_size, marker="s", color = '#d08770')[0]

        L3 = ax2.plot(pi, bcf_size , marker="p", color = '#8fbcbb')[0]

        L4 = ax2.plot(pi, bgt_size , marker="<", color = '#a3be8c')[0]

        ax2.xaxis.set_major_formatter(StrMethodFormatter('{x:.0f}'))
        ax2.xaxis.set_minor_formatter(NullFormatter())
        ax2.tick_params(axis='x', labelrotation = 65)
        ax2.set_title(f"1KGP file sizes", fontweight="bold")
        ax2.set_xticks(pi)
        ax2.set_xticklabels(p)
        ax2.set_xlabel("Chromosome", fontweight="bold")
        ax2.set_ylabel("Memory (gigabytes)", fontweight="bold")
        ax2.semilogy()
        ax2.yaxis.set_major_formatter(ticker.FuncFormatter(myLogFormat))

        plt.tight_layout()

        line_labels = ["μ-PBWT serialization", "Syllable-PBWT serialization", "BCF file", "BGT files"]
        fig.legend(handles= [L1, L2, L3, L4],     # The line objects
                   labels=line_labels,   # The labels for each line
                   loc="lower center",   # Position of legend
                   #borderaxespad=0.1,    # Small spacing around legend box
                   bbox_to_anchor=(0.5, -0.05),
                   ncol=4
                   #title="Legend Title"  # Title for the legend
                   )
        plt.savefig(f"{plot_dir}/files_size_bgt.pdf", dpi=500, bbox_inches='tight')

        
        pm = np.array(phi_mem)
        mm = np.array(map_mem)
        sm = np.array(sam_mem)
        tm = np.array(thr_mem)
        
        fig, (ax2) = plt.subplots(1, 1, figsize=(5,5))

        L4 = ax2.bar(pi, pm+mm+sm+tm  , color = '#ebcb8b',  label="Thresholds")[0]
        L3 = ax2.bar(pi,  pm+mm+sm , color = '#8fbcbb', label="PA/DA samples")[0]
        L2 = ax2.bar(pi, pm+mm, color = '#d08770', label= "Mapping structure")[0]
        L1 = ax2.bar(pi, pm ,  color = '#b48ead',  label="Φ structure")[0]

        ax2.xaxis.set_major_formatter(StrMethodFormatter('{x:.0f}'))
        ax2.xaxis.set_minor_formatter(NullFormatter())
        ax2.tick_params(axis='x', labelrotation = 65)
        ax2.set_title(f"a) 1KGP components sizes", fontweight="bold")
        ax2.set_xticks(pi)
        ax2.set_xticklabels(p)
        ax2.set_xlabel("Chromosome", fontweight="bold")
        ax2.set_ylabel("Memory (megabytes)", fontweight="bold")
        ax2.semilogy()


        line_labels = ["Φ structure", "Mapping structure", "PA/DA samples", "Thresholds"]
        fig.legend(handles= [L1, L2, L3, L4],     # The line objects
                   labels=line_labels,   # The labels for each line
                   loc="lower center",   # Position of legend
                   #borderaxespad=0.1,    # Small spacing around legend box
                   bbox_to_anchor=(0.5, -0.07),
                   ncol=4
                   #title="Legend Title"  # Title for the legend
                   )
        plt.savefig(f"{plot_dir}/components.pdf", dpi=500, bbox_inches='tight')
        
if __name__ == "__main__":
    main(sys.argv[1:])
