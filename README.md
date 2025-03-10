# μ-PBWT with k-SMEMs
Phase genotype with μ-PBWT.

## Build from source
Prepare the cmake for building the current project in ‘.’ into the ‘build’ folder
```shell
cmake -S . -B build 
```
Build μ-PBWT:
```shell
cmake --build build
```


## Usage
File format supported:
- BCF/VCF

```shell
cd build
```
```shell

Usage: ./mupbwt [options]

Options:
  -i, --input_file <path>	vcf/bcf file for panel
  -r, --ref <path>	vcf/bcf file for full reference (for imputation)
  -o, --output <path>	path to phased bcf
  -q, --query <path>	path target vcf/bcf
  -u, --unphased	leave ambigous sites unphased
  -d, --details	print memory usage details
  -v, --verbose	extra prints
  -h, --help	show this help message and exit
```
If input and ref panels are the same we perform phasing, unless (the ref panel
presents more variation sites) we try impute the target.

Phase (not impute) `sample_data/target_unphased.vcf` using `sample_data/reference.vcf` as reference:
```shell
./mupbwt -i sample_data/reference.vcf -r sample_data/reference.vcf -q sample_data/target_unphased.vcf -o target_phased.bcf -u
```



