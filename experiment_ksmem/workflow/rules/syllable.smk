##### it will probbaly beak if they update server.cpp but there is not a tag or release #####
rule downloadSyllablle:
    output:
        d=directory(syllable_folder),
        exe=os.path.join(syllable_folder, "server")
    conda: "../envs/compilation.yml"
    shell:
        """
        git clone https://github.com/ZhiGroup/Syllable-PBWT {output.d}
        cd {output.d}
        awk \'NR==123{{print \"exit(0);\"}}1\' server.cpp > tmp && mv tmp server.cpp
        make 
        """

rule makeSyllableInput:
    input:
        panel = os.path.join(input_folder, "{chr}", "panel_{chr}.bcf"),
    output:
        vcf =  os.path.join(output_folder, "syllable", "{chr}", "panel_{chr}.vcf"),
    conda: "../envs/bcftools.yml"
    shell:
        """
        bcftools view -m2 -M2 -v snps {input.panel} > {output.vcf}
        """

rule makeSyllable:
    input:
        exe = os.path.join(syllable_folder, "server"),
        vcf =  os.path.join(output_folder, "syllable", "{chr}", "panel_{chr}.vcf"),
    params:
        vcf = os.path.join("panel_{chr}.vcf"),
        fifo = os.path.join("chr_{chr}"),
        fol = os.path.join(output_folder, "syllable", "{chr}"),
    output:
        os.path.join(output_folder, "syllable", "{chr}", "chr_{chr}-save0.txt"),
        os.path.join(output_folder, "syllable", "{chr}", "make.time")
    shell:
        """
        cd {params.fol}
        /usr/bin/time --verbose -o make.time ../../../../{input.exe} -f {params.fifo} -i {params.vcf} -s &> exe.log
        rm {params.vcf}
        """
