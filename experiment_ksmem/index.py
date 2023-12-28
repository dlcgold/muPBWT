import os
import subprocess
import asyncio
import time

def subprocess_cmd(command):
    process = subprocess.Popen(command,stdout=subprocess.PIPE, shell=True)
    proc_stdout = process.communicate()[0].strip()
    print(proc_stdout)
    
def background(f):
    def wrapped(*args, **kwargs):
        return asyncio.get_event_loop().run_in_executor(None, f, *args, **kwargs)

    return wrapped

@background
def op(name):
    #print(f"bcftools {name}")
    subprocess_cmd(f"mkdir resources/input_index/{name}")
    subprocess_cmd(f"curl -L https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr{name}.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz > resources/input_index/{name}/{name}.vcf.gz")
    subprocess_cmd(f"bcftools view -m2 -M2  -v snps resources/input_index/{name}/{name}.vcf.gz -Ob > resources/input_index/{name}/{name}.bcf")
    #print(f"mupbwt {name}")
    subprocess_cmd(f"mupbwt -i resources/input_index/{name}/{name}.bcf -s 1kgp_index/{name}.ser &> /dev/null")
    subprocess_cmd(f"tar -czvf 1kgp_index_gz/{name}.ser.tar.gz  1kgp_index/{name}.ser")
    #subprocess_cmd(f"rm -r resources/input_index/{name}")
curr_dir = "."
panels = ['1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15', '16', '17', '18', '19', '20', '21', '22']
# for name in dirs:
#     op(name)

loop = asyncio.get_event_loop()                                              # Have a new event loop

looper = asyncio.gather(*[op(i) for i in panels])         # Run the loop
                               
results = loop.run_until_complete(looper) 
# for i in panels:
#     op(i)
