rule extractVerboseQuery:
    input:
        os.path.join(bench_folder, "time", "{tool}", "mutate_{q}", "prediction.time"),
    output:
        os.path.join(bench_folder, "time", "{tool}", "mutate_{q}", "prediction.csv"),
    shell:
        """
        python workflow/scripts/time_verbose_extractor.py {wildcards.tool} {wildcards.q} < {input} > {output}
        """

rule Hapcut2Csv:
    input:
       os.path.join(bench_folder, "ser_hapcut2", "{tool}", "mutate_{q}", "prediction.ser.txt"),
    output:
       os.path.join(bench_folder, "ser_hapcut2", "{tool}", "mutate_{q}", "prediction.ser.csv"),
    shell:
        """
        python workflow/scripts/hapcut2_ser_extractor.py {wildcards.tool} {wildcards.q} < {input} > {output}
        """

rule mergeTime:
    input:
        expand(
            os.path.join(bench_folder, "time", "{tool}", "mutate_{q}", "prediction.csv"),
            tool = ["mupbwt", "mupbwt_u", "beagle"],
            q = mutations_perc
        ),
    output:
        os.path.join(results_folder, "time.csv"),
    conda: "../envs/csvkit.yml"
    shell:
        """
        csvstack {input} > {output}
        """

rule mergeHapcut2:
    input:
        expand(
            os.path.join(bench_folder, "ser_hapcut2", "{tool}", "mutate_{q}", "prediction.ser.csv"),
            tool = ["mupbwt", "mupbwt_u", "beagle"],
            q = mutations_perc
        ),
    output:
        os.path.join(results_folder, "hapcut2.csv"),
    conda: "../envs/csvkit.yml"
    shell:
        """
        csvstack {input} > {output}
        """
