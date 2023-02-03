[![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?style=flat)](http://bioconda.github.io/recipes/mupbwt/README.html)
[![Conda](https://img.shields.io/conda/v/bioconda/mupbwt?color=green)](https://anaconda.org/bioconda/mupbwt)
[![Conda](https://img.shields.io/conda/dn/bioconda/mupbwt?color=green&label=conda%20%7C%20downloads)](https://anaconda.org/bioconda/mupbwt)
[![GitHub stars](https://img.shields.io/github/stars/dlcgold/muPBWT.svg)](https://github.com/dlcgold/muPBWT/stargazers)
# μ-PBWT
A PBWT-based light index  for UK Biobank scale genotype data.

## Conda install
μ-PBWT is available for Gnu/Linux on [conda](https://docs.conda.io/en/latest/) ([bioconda](https://bioconda.github.io/) channel):
```shell
conda install -c bioconda mupbwt
```

## Build from source
Prepare the cmake for building the current project in ‘.’ into the ‘build’ folder
```shell
cmake -S . -B build 
```
Build μ-PBWT:
```shell
cmake --build build
```
## Install from source (optional)
Install μ-PBWT (default in `/usr/local/bin/`, `sudo` required):
```shell
cmake --install build
```
Use `--prefix <path>` for custom path.

## Usage
File format supported:
- BCF/VCF
- [MaCS](https://github.com/gchen98/macs)
```shell
cd build
```
```shell

Usage: ./mupbwt [options]

Options:
  -i, --input_file <path>	 vcf/bcf file for panel
  -s, --save <path>	  path to save index
  -l, --load <path>	 path to load index
  -o, --output <path>	 path to query output
  -q, --query <path>	 path to query file (vcf/bcf)
  -m, --macs	use macs as file format for both input and query file
  -v, --verbose	 extra prints
  -d, --details	 print memory usage details
  -h, --help	 show this help message and exit
```

Build the index:
```shell
./mupbwt -i <input file> -s <index file>
```
Query the index:
```shell
./mupbwt -l <index file> -q <query file> -o <output file> 
```
Query without save the index:
```shell
./mupbwt -i <input file> -q <query file> -o <output file>
```
Query and  save the index:
```shell
./mupbwt -i <input file> -s <index file> -q <query file> -o <output file>
```
Using examples in `sample_data`:
```shell
./mupbwt -i sample_data/panel.bcf -s sample_data/index.ser
./mupbwt -l sample_data/index.ser -q sample_data/query.bcf -o sample_data/sample_data_results 
./mupbwt -i sample_data/panel.bcf -q sample_data/query.bcf -o sample_data/sample_data_results
./mupbwt -i sample_data/panel.bcf -s sample_data/index.ser -q sample_data/query.bcf -o sample_data/sample_data_results
```

Load the index and print details to stdout:
```shell
./mupbwt -l <index file> -d
```
An output example is:
```shell
> ./mupbwt -l sample_data/index.ser -d
built/loaded in: 0.015628 s

----
Total haplotypes: 900
Total sites: 499
----
Total runs: 27512
Average runs: 55
----
run: 0.0386925 megabytes
thr: 0.0387306 megabytes
uv: 0.0380135 megabytes
samples: 0.0833178 megabytes
rlpbwt (mapping): 0.201148 megabytes
phi panels: 0.414757 megabytes
phi support: 0.126385 megabytes
phi data structure (panels + support): 0.541142 megabytes
rlpbwt: 0.74229 megabytes
----
estimated dense size: 36.4132 megabytes
----
```

### Input
Only bialleic case is supported. In case of vcf/bcf [bcftools](https://github.com/samtools/bcftools) can be used to filter the input:
```shell
bcftools view -m2 -M2  -v snps <input vcf/vcf> > <filtered vcf/bcf file>
```
### Output
Output file follow the standard proposed in [Durbin's PBWT](https://github.com/richarddurbin/pbwt). 
Each row contain a SMEM:
```
MATCH   <query index>   <row index> <staring column>    <ending column> <SMEM length>
```
For example:
```
MATCH	99	150	414	430	17
```
Row index and query index are incrementally so the name of the sample and the precise haplotype can be calculated using the output of [bcftools](https://github.com/samtools/bcftools). 

The command:
```shell
bcftools query -l <input vcf/bcf> > samples.txt
```
store in `samples.txt` all the samples name, in order. So, for example, row indices 0 and 1 corresponds to the two haplotypes of the first sample, row indices 2 and 3 to the second one etc...

Eventually you can use `script/mem_sample.py`:
```shell
> python mem_sample.py -h
usage: mem_sample.py [-h] [-i INPUT] [-p PANEL] [-q QUERIES] [-o OUTPUT]

options:
  -h, --help            show this help message and exit
  -i INPUT, --input INPUT
                        SMEM file in Durbin's format
  -p PANEL, --panel PANEL
                        panel as VCF/BCF (optional)
  -q QUERIES, --queries QUERIES
                        queries as VCF/BCF (optional)
  -o OUTPUT, --output OUTPUT
                        output file
```
Esample:
```
python mem_sample.py -i sample_data/sample_data_results -p sample_data/panel.bcf -q sample_data/query.bcf -o sample_data/sample_data_results_new
```
Only one between `PANEL` and `QUERIES` can be specified.
New SMEM file will contain in each row:
```
MATCH   <querySample_haplotype>   <panelSample_haplotype> <staring column>    <ending column> <SMEM length>
```
For example, assuming both `PANEL` and `QUERIES`:
```
MATCH	1318026_1	4919834_0	414	430	17
```
## Results
Results on high-coverage whole genome sequencing data from UK Biobank (chromosome 20):

| **Region**              | **#Samples** | **#Sites**   | **Size BCF (GB)** | **μ-PBWT (GB)** | **Construction time (hh:mm)** | **Construction memory peak (GB)** |
|-------------------------|--------------|--------------|-------------------|-----------------|-------------------------------|-----------------------------------|
| chr20:60061-4060065     | 150119       | 865267       | 1.9               | 0.88            | 06:25                         | 2.27                              |
| chr20:4060066-8060066   | 150119       | 880899       | 2                 | 0.85            | 06:28                         | 2.22                              |
| chr20:8060067-12515479  | 150119       | 961591       | 2.1               | 0.77            | 07:04                         | 2.05                              |
| chr20:12515480-16768988 | 150119       | 917468       | 2                 | 0.73            | 06:47                         | 1.97                              |
| chr20:16768989-21050967 | 150119       | 931010       | 2                 | 0.71            | 06:53                         | 1.92                              |
| chr20:21050968-31549151 | 150119       | 1919134      | 4.2               | 1.20            | 13:54                         | 3.06                              |
| chr20:31549152-38282825 | 150119       | 1436549      | 2.8               | 0.99            | 10:25                         | 2.63                              |
| chr20:38282826-43181963 | 150119       | 1056144      | 2.2               | 0.76            | 07:42                         | 2.06                              |
| chr20:43181964-47619489 | 150119       | 955970       | 2                 | 0.79            | 06:56                         | 2.09                              |
| chr20:47619490-51789198 | 150119       | 923178       | 2                 | 0.80            | 06:44                         | 2.12                              |
| chr20:51789199-55789212 | 150119       | 911452       | 2                 | 0.81            | 06:45                         | 2.13                              |
| chr20:55789213-59874964 | 150119       | 925442       | 2                 | 0.84            | 06:49                         | 2.20                              |
| chr20:59874965-64334101 | 150119       | 1096089      | 2.4               | 0.93            | 08:00                         | 2.42                              |
| **Total**               | **150119**   | **13780193** | **29.6**          | **11.08**       | **-**                         | **29.15**                         |

Results on 1000 Genome Project phase 3 data including the average number of runs per site:

| **Chr**   | **#Samples** | **#Sites**   | **#Runs/site** | **Size BCF (GB)** | **μ-PBWT (GB)** | **Construction time (hh:mm)** | **Construction memory peak (GB)** |
|-----------|--------------|--------------|----------------|-------------------|-----------------|-------------------------------|-----------------------------------|
| 1         | 2454         | 6196151      | 11             | 0.78              | 1.44            | 00:19                         | 4.59                              |
| 2         | 2454         | 6786300      | 10             | 0.84              | 1.47            | 00:21                         | 4.76                              |
| 3         | 2454         | 5584397      | 10             | 0.71              | 1.20            | 00:18                         | 4.24                              |
| 4         | 2454         | 5480936      | 10             | 0.71              | 1.19            | 00:17                         | 4.28                              |
| 5         | 2454         | 5037955      | 9              | 0.63              | 1.08            | 00:16                         | 4.22                              |
| 6         | 2454         | 4800101      | 10             | 0.64              | 1.06            | 00:15                         | 4.28                              |
| 7         | 2454         | 4517734      | 10             | 0.58              | 1.03            | 00:14                         | 4.34                              |
| 8         | 2454         | 4417368      | 10             | 0.56              | 0.97            | 00:14                         | 4.30                              |
| 9         | 2454         | 3414848      | 11             | 0.43              | 0.81            | 00:11                         | 2.54                              |
| 10        | 2454         | 3823786      | 10             | 0.50              | 0.87            | 00:12                         | 2.77                              |
| 11        | 2454         | 3877543      | 10             | 0.49              | 0.84            | 00:12                         | 2.71                              |
| 12        | 2454         | 3698099      | 10             | 0.47              | 0.82            | 00:12                         | 2.63                              |
| 13        | 2454         | 2727881      | 10             | 0.35              | 0.60            | 00:9                          | 2.14                              |
| 14        | 2454         | 2539149      | 11             | 0.32              | 0.58            | 00:8                          | 2.18                              |
| 15        | 2454         | 2320474      | 12             | 0.29              | 0.57            | 00:7                          | 2.30                              |
| 16        | 2454         | 2596072      | 12             | 0.32              | 0.63            | 00:8                          | 2.28                              |
| 17        | 2454         | 2227080      | 12             | 0.28              | 0.55            | 00:7                          | 2.32                              |
| 18        | 2454         | 2171378      | 11             | 0.28              | 0.51            | 00:7                          | 2.23                              |
| 19        | 2454         | 1751878      | 13             | 0.23              | 0.45            | 00:6                          | 1.43                              |
| 20        | 2454         | 1739315      | 11             | 0.22              | 0.41            | 00:5                          | 1.30                              |
| 21        | 2454         | 1054447      | 14             | 0.14              | 0.30            | 00:3                          | 1.26                              |
| 22        | 2454         | 1055454      | 14             | 0.14              | 0.29            | 00:3                          | 1.24                              |
| **Total** | **2454**     | **77818346** | **11**         | **9.91**          | **17.67**       | **-**                         | **64.34**                         |

Note that total building times are not printed due to the fact that all the computations have been done in parallel.

The pipeline for 1000 Genome Project phase 3 data is available at [dlcgold/muPBWT-1KGP-workflow](https://github.com/dlcgold/muPBWT-1KGP-workflow).

WIP
