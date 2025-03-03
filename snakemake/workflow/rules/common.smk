rule extractTarget:
    input:
        os.path.join(reference),
    output:
        refgz = os.path.join(reference + ".gz"),
        refcsi = os.path.join(reference + ".gz.csi"),
        vcf = os.path.join(input_folder, "target.vcf"),
        vcfgz = os.path.join(input_folder, "target.vcf.gz"),
        csi = os.path.join(input_folder, "target.vcf.gz.csi"),
    params:
        t = os.path.join(input_folder, "target"),
        s = os.path.join(input_folder, "samples"),
    conda:
         "../envs/gen.yml"
    shell:
        """
        bcftools view -Oz {input} > {output.refgz}
        bcftools index {output.refgz}
        bcftools query -l {input} > {params.s}
        shuf -n 1 {params.s} > {params.t} 
        #tail -n1 {params.s} > {params.t}
        bcftools view -S {params.t} {input} > {output.vcf}
        bcftools view -Oz {output.vcf} >  {output.vcfgz}
        bcftools index {output.vcfgz}
        #rm  {params.s}  {params.t}
        """

rule createGenotype:
    input:
        vcf = os.path.join(input_folder, "target.vcf"),
    output:
        vcf = os.path.join(input_folder, "target_unphased.vcf"),
        vcfgz = os.path.join(input_folder, "target_unphased.vcf.gz"),
        csi = os.path.join(input_folder, "target_unphased.vcf.gz.csi"),
    conda:
        "../envs/gen.yml"
    shell:
        """
        sed -E '/^##/! s/\\b1\\|0\\b/0\\/1/g; \
                /^##/! s/\\b0\\|1\\b/0\\/1/g; \
                /^##/! s/\\b1\\|1\\b/1\\/1/g; \
                /^##/! s/\\b0\\|0\\b/0\\/0/g' {input} > {output.vcf}
        bcftools view -Oz {output.vcf} >  {output.vcfgz}
        bcftools index {output.vcfgz}
        """

rule mutateGenotype:
    input:
        os.path.join(input_folder, "target_unphased.vcf"),
    output:
        vcf = os.path.join(input_folder, "mutate_{q}", "target_unphased.vcf"),
        vcfgz = os.path.join(input_folder, "mutate_{q}", "target_unphased.vcf.gz"),
        csi = os.path.join(input_folder, "mutate_{q}", "target_unphased.vcf.gz.csi"),
    shell:
        """
        python workflow/scripts/mutate.py  {input} {wildcards.q} {output.vcf}
        bcftools view {output.vcf} -Oz > {output.vcfgz}
        bcftools index {output.vcfgz}
        """
