rule makePanelPbwt:
    input:
        panel =  os.path.join(input_folder, "{chr}", "panel_{chr}.bcf"),
    output:
        panel = os.path.join(output_folder, "pbwt", "{chr}", "panel.pbwt"),
    params:
         panel = os.path.join(output_folder, "pbwt", "{chr}", "panel.sites"),
    shadow: "shallow"
    log:
        panel = os.path.join(output_folder, "pbwt", "{chr}",  "panel.time"),
	panel_l = os.path.join(output_folder, "pbwt", "{chr}",  "panel.log"),
    conda: "../envs/pbwt.yml"
    shell:
        """
        /usr/bin/time --verbose -o {log.panel} pbwt -readVcfGT {input.panel} -write {output.panel} -writeSites {params.panel} &> {log.panel_l}
        """
rule makeQueryPbwt:
    input:
        query =  os.path.join(input_folder, "{chr}", "query_{chr}.bcf"),
    output:
        query = os.path.join(output_folder, "pbwt", "{chr}", "query.pbwt"),
    params:
         query = os.path.join(output_folder, "pbwt", "{chr}", "query.sites"),
    shadow: "shallow"
    log:
        query = os.path.join(output_folder, "pbwt", "{chr}",  "query.time"),
        query_l = os.path.join(output_folder, "pbwt", "{chr}",  "query.log"),
    conda: "../envs/pbwt.yml"
    shell:
        """
        /usr/bin/time --verbose -o {log.query} pbwt -readVcfGT {input.query} -write {output.query} -writeSites {params.query} &> {log.query_l}
        """

##### thanks to https://stackoverflow.com/questions/51977436/restrict-number-of-jobs-by-a-rule-in-snakemake for load trick #####
rule runPbwtIndexed:
    input:
        panel = os.path.join(output_folder, "pbwt", "{chr}", "panel.pbwt"),
        query = os.path.join(output_folder, "pbwt", "{chr}", "query.pbwt"),
    output:
        out = os.path.join(output_folder, "pbwt", "{chr}", "pbwtIndexed.txt"),
    log:
        log = os.path.join(output_folder, "pbwt", "{chr}", "pbwtIndexed.log"),
        time = os.path.join(output_folder, "pbwt", "{chr}", "pbwtIndexed.time"),
    conda: "../envs/pbwt.yml"
    resources:
        load=100
    shell:
        """
        /usr/bin/time --verbose -o {log.time} pbwt -read {input.panel} -matchIndexed {input.query} > {output.out} 2> {log.log}
        """

rule runPbwtDynamic:
    input:
        panel = os.path.join(output_folder, "pbwt", "{chr}", "panel.pbwt"),
        query = os.path.join(output_folder, "pbwt", "{chr}", "query.pbwt"),
    output:
        out = os.path.join(output_folder, "pbwt", "{chr}", "pbwtDynamic.txt"),
    log:
        log = os.path.join(output_folder, "pbwt", "{chr}", "pbwtDynamic.log"),
        time = os.path.join(output_folder, "pbwt", "{chr}", "pbwtDynamic.time"),
    conda: "../envs/pbwt.yml"
    shell:
        """
        /usr/bin/time --verbose -o {log.time} pbwt -read {input.panel} -matchDynamic {input.query} > {output.out} 2> {log.log}
        """
