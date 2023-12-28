##### TODO REMOVE vcf.gz #####
rule download1KGP:
    params:
        url = "https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr{chr}.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz",
        full_s = os.path.join(input_folder, "{chr}", "{chr}.txt"),
        nqueries = str(int(int(n_queries) / 2))
    output:
        gz = os.path.join(input_folder,  "{chr}", "{chr}.vcf.gz"),
        panel_s = os.path.join(input_folder, "{chr}", "panel_{chr}.txt"),
        query_s = os.path.join(input_folder, "{chr}", "query_{chr}.txt")
    conda: "../envs/bcftools.yml"
    shell:
        """
        curl -L {params.url} -s > {output.gz}
        bcftools query -l {output.gz} > {params.full_s}
        (head -{params.nqueries} > {output.query_s}; cat > {output.panel_s}) < {params.full_s}
        rm {params.full_s}
        """

        
rule makePanel1KGP:
    input:
        gz = os.path.join(input_folder,  "{chr}", "{chr}.vcf.gz"),
        panel_s = os.path.join(input_folder, "{chr}", "panel_{chr}.txt"),
    output:
        panel = os.path.join(input_folder, "{chr}", "panel_{chr}.bcf"), 
    conda: "../envs/bcftools.yml"
    shell:
        """
        bcftools view -m2 -M2  -v snps -S {input.panel_s} {input.gz} -Ob  > {output.panel}
        rm {input.panel_s}
        """

rule makeQuery1KGP:
    input:
        gz = os.path.join(input_folder,  "{chr}", "{chr}.vcf.gz"),
        query_s = os.path.join(input_folder, "{chr}", "query_{chr}.txt"),
    output:
        query = os.path.join(input_folder, "{chr}", "query_{chr}.bcf"), 
    conda: "../envs/bcftools.yml"
    shell:
        """
        bcftools view -m2 -M2  -v snps -S {input.query_s} {input.gz} -Ob  > {output.query}
        rm {input.query_s}
        """
