rule runBeagle:
    input:
        ref = os.path.join(reference + ".gz"),
        map_b = os.path.join(beagle_map),
        tar = os.path.join(input_folder, "mutate_{q}", "target_unphased.vcf.gz"),
    output:
        vcf = os.path.join(output_folder, "beagle", "mutate_{q}", "prediction.vcf"),
        vcfgz = os.path.join(output_folder, "beagle", "mutate_{q}", "prediction.vcf.gz"),
        csi = os.path.join(output_folder, "beagle", "mutate_{q}", "prediction.vcf.gz.csi"),
    params:
        p = os.path.join(output_folder, "beagle", "mutate_{q}", "prediction"),
    log:
        os.path.join(bench_folder, "time", "beagle", "mutate_{q}", "prediction.time"),
    conda:
        "../envs/gen.yml"
    resources:
        load=100
    shell:
        """
        /usr/bin/time -vo {log} java -Djava.io.tmpdir=tmp/ -Xmx160g -jar  {beagle} chrom=chr{chrom} gt={input.tar} ref={input.ref} out={params.p}  map={input.map_b} nthreads={threads}
        bcftools index {output.vcfgz}
        bcftools view {output.vcfgz} >  {output.vcf}
        """
