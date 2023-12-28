rule makeMupbwt:
    input:
        panel = os.path.join(input_folder, "{chr}", "panel_{chr}.bcf"),
    output:
        ser = os.path.join(output_folder, "mupbwt", "{chr}_{k}", "{chr}.ser"),
    log:
        log = os.path.join(output_folder, "mupbwt", "{chr}_{k}", "make.log"),
        time = os.path.join(output_folder, "mupbwt", "{chr}_{k}", "make.time"),
    conda: "../envs/mupbwt.yml"
    shell:
        """
        /usr/bin/time --verbose -o {log.time} mupbwt -i {input.panel} -s {output.ser} -k {wildcard.k} &> {log.log}
        """
        
rule runMupbwt:
    input:
        mupbwt = os.path.join(output_folder, "mupbwt", "{chr}_{k}", "{chr}.ser"),
        query = os.path.join(input_folder, "{chr}_{k}", "query_{chr}.bcf"),
    output:
        res = os.path.join(output_folder, "mupbwt", "{chr}_{k}", "output.txt"),
        stat = os.path.join(output_folder, "mupbwt", "{chr}_{k}", "stat.txt"),
    log:
        log = os.path.join(output_folder, "mupbwt", "{chr}_{k}", "exe.log"),
        time = os.path.join(output_folder, "mupbwt", "{chr}_{k}", "exe.time"),
    conda: "../envs/mupbwt.yml"
    shell:
        """
        OMP_NUM_THREADS=1 /usr/bin/time --verbose -o {log.time}  mupbwt -l {input.mupbwt} -q {input.query} -o {output.res} -k {wildcard.k} &> {log.log}
        mupbwt -l {input.mupbwt} -d -k {wildcard.k} > {output.stat}
        """
