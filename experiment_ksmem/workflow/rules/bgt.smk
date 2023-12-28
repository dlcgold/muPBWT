rule makeBgt:
    input:
        panel = os.path.join(input_folder, "{chr}", "panel_{chr}.bcf"),
    params:
        bgt = os.path.join(output_folder, "bgt", "{chr}", "{chr}.bgt"),
    output:
         os.path.join(output_folder, "bgt", "{chr}", "{chr}.bgt.bcf"),
         os.path.join(output_folder, "bgt", "{chr}", "{chr}.bgt.bcf.csi"),
         os.path.join(output_folder, "bgt", "{chr}", "{chr}.bgt.pbf"),
         os.path.join(output_folder, "bgt", "{chr}", "{chr}.bgt.spl"),
    log:
        log = os.path.join(output_folder, "bgt", "{chr}", "make.log"),
        time = os.path.join(output_folder, "bgt", "{chr}", "make.time"),
    conda: "../envs/bgt.yml"
    shell:
        """
        /usr/bin/time --verbose -o {log.time} bgt import {params.bgt} {input.panel} &> {log.log}
        """
        
