rule serMupbwt:
    input:
        vcf = os.path.join(output_folder, "mupbwt", "mutate_{q}", "prediction.vcf"),
        ref = os.path.join(input_folder, "target.vcf"),
    output:
        os.path.join(bench_folder, "ser_hapcut2", "mupbwt", "mutate_{q}", "prediction.ser.txt"),
    shell:
        """
        python workflow/scripts/calculate_haplotype_statistics.py -v1 {input.vcf} -v2 {input.ref} > {output}
        """

rule serMupbwt_u:
    input:
        vcf = os.path.join(output_folder, "mupbwt_u", "mutate_{q}", "prediction.vcf"),
        ref = os.path.join(input_folder, "target.vcf"),
    output:
        os.path.join(bench_folder, "ser_hapcut2", "mupbwt_u", "mutate_{q}", "prediction.ser.txt"),
    shell:
        """
        python workflow/scripts/calculate_haplotype_statistics.py -v1 {input.vcf} -v2 {input.ref} > {output}
        """


rule serBeagle:
    input:
        vcf = os.path.join(output_folder, "beagle", "mutate_{q}", "prediction.vcf"),
        ref = os.path.join(input_folder, "target.vcf"),
    output:
        os.path.join(bench_folder, "ser_hapcut2", "beagle", "mutate_{q}", "prediction.ser.txt"),
    shell:
        """
        python workflow/scripts/calculate_haplotype_statistics.py -v1 {input.vcf} -v2 {input.ref} > {output}
        """

