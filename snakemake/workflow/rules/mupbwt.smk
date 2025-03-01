rule runMupbwt:
    input:
        ref = os.path.join(reference),
        tar = os.path.join(input_folder, "mutate_{q}", "target_unphased.vcf"),
    output:
        bcf = os.path.join(output_folder, "mupbwt", "mutate_{q}", "prediction.bcf"),
        vcf = os.path.join(output_folder, "mupbwt", "mutate_{q}", "prediction.vcf"),
        vcfgz = os.path.join(output_folder, "mupbwt", "mutate_{q}", "prediction.vcf.gz"),
        csi = os.path.join(output_folder, "mupbwt", "mutate_{q}", "prediction.vcf.gz.csi"),
    log:
        os.path.join(bench_folder, "time", "mupbwt", "mutate_{q}", "prediction.time"),
    conda:
        "../envs/gen.yml"
    shell:
        """
        /usr/bin/time -vo {log} {mupbwt} -i {input.ref} -r {input.ref} -q {input.tar} -o {output.bcf} -u &> /dev/null
        bcftools view {output.bcf} > {output.vcf}
        bcftools view {output.vcf} -Oz > {output.vcfgz}
        bcftools index {output.vcfgz}
        """


rule runMupbwt_u:
    input:
        ref = os.path.join(reference),
        tar = os.path.join(input_folder, "mutate_{q}", "target_unphased.vcf"),
    output:
        bcf = os.path.join(output_folder, "mupbwt_u", "mutate_{q}", "prediction.bcf"),
        vcf = os.path.join(output_folder, "mupbwt_u", "mutate_{q}", "prediction.vcf"),
        vcfgz = os.path.join(output_folder, "mupbwt_u", "mutate_{q}", "prediction.vcf.gz"),
        csi = os.path.join(output_folder, "mupbwt_u", "mutate_{q}", "prediction.vcf.gz.csi"),
    log:
        os.path.join(bench_folder, "time", "mupbwt_u", "mutate_{q}", "prediction.time"),
    conda:
        "../envs/gen.yml"
    shell:
        """
        /usr/bin/time -vo {log} {mupbwt} -i {input.ref} -r {input.ref} -q {input.tar} -o {output.bcf} &> /dev/null
        bcftools view {output.bcf} > {output.vcf}
        bcftools view {output.vcf} -Oz > {output.vcfgz}
        bcftools index {output.vcfgz} 
        """

