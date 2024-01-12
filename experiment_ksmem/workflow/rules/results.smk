rule extractVerbose:
    input:
        #os.path.join(output_folder, "pbwt", "{chr}",  "panel.time"),
        #os.path.join(output_folder, "pbwt", "{chr}",  "query.time"),
        #os.path.join(output_folder, "pbwt", "{chr}", "pbwtIndexed.time"),
        #os.path.join(output_folder, "pbwt", "{chr}", "pbwtDynamic.time"),
        os.path.join(output_folder, "mupbwt", "{chr}_{k}", "make.time"),
        os.path.join(output_folder, "mupbwt", "{chr}_{k}", "exe.time"),
        #os.path.join(output_folder, "syllable", "{chr}", "make.time"),
        time = os.path.join(output_folder, "mupbwt", "{chr}_{k}", "{file}.time"),
        #nqueries = n_queries
    output:
        os.path.join(output_folder, "mupbwt", "{chr}_{k}", "{file}.time.csv")
    shell:
        """
        python workflow/scripts/time_verbose_extractor.py {wildcards.k} 4908 0 {wildcards.chr} {n_queries}  < {input.time} > {output}
        """

rule mergeExeTime:
    input:
        #expand(
        #    os.path.join(output_folder, "pbwt", "{chr}", "{approach}.time.csv"),
        #    approach = APPROACHES_PBWT,
        #    chr = PANELS,
        #),
        expand(
            os.path.join(output_folder, "mupbwt", "{chr}_{k}", "exe.time.csv"),
            chr = PANELS,
            k = KS,
        )
    output:
        os.path.join(results_folder, "data", "exe-time.csv"),
    conda: "../envs/csvkit.yml"
    shell:
        """
        csvstack {input} > {output}
        """

rule mergeMakeTime:
    input:
        #expand(
        #    os.path.join(output_folder, "pbwt", "{chr}", "panel.time.csv"),
        #    chr = PANELS,
        #),
        expand(
            os.path.join(output_folder, "mupbwt", "{chr}_{k}", "make.time.csv"),
            chr = PANELS,
            k = KS,
        ),
       # expand(
       #     os.path.join(output_folder, "syllable", "{chr}", "make.time.csv"),
       #     chr = PANELS,
       # ),
    output:
        os.path.join(results_folder, "data", "make-time.csv"),
    conda: "../envs/csvkit.yml"
    shell:
        """
        csvstack {input} > {output}
        """
        
rule makePlot:
    input:
         os.path.join(results_folder, "data", "make-time.csv"),
         os.path.join(results_folder, "data", "exe-time.csv"),
    params:
        input_dir = input_folder,
        output_dir = output_folder,
        results_dir = results_folder,
    output:
        os.path.join(results_folder, "plots", "1kgb.pdf"),
        os.path.join(results_folder, "plots", "1kgb_nodu.pdf"),
        os.path.join(results_folder, "plots", "components.pdf"),
        os.path.join(results_folder, "plots", "files_size.pdf"),
        os.path.join(results_folder, "tables", "mupbwt_build.tex"),
        os.path.join(results_folder, "tables", "mupbwt_compare.tex"),
    conda: "../envs/results.yml"
    shell:
        """
        python workflow/scripts/results.py {params.input_dir} {params.output_dir} {params.results_dir}
        """
