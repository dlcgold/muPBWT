# μ-PBWT with k-SMEMs
Compute k-SMEMs with μ-PBWT.
At [branch](https://github.com/dlcgold/muPBWT/tree/k-smem-live) implementation with k
values not pre-computed (less memory but slower).

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
  -k, --ksmem <k value> compute k-smems 
  -v, --verbose	 extra prints
  -d, --details	 print memory usage details
  -h, --help	 show this help message and exit
```

Query the index for kSMEMs with pre-computed values:
```shell
./mupbwt -l <index file> -q <query file> -o <output file> -k <k value>
```

Using examples in `sample_data`:
```shell
./mupbwt -i sample_data/panel.bcf -s sample_data/index.ser -q sample_data/query.bcf -o sample_data/sample_data_results -k 30
```



## Reference
Bibtex:
```
@article {smems,
	author = {Bonizzoni, Paola and Boucher, Christina and Cozzi, Davide and  Gagie, Travis and Pirola, Yuri }},
	title = {Solving the Minimal Positional Substring Cover problem in sublinear space},
	year = {2024},
	doi = {10.4230/LIPIcs.CPM.2024.2},
	booktitle={35th Annual Symposium on Combinatorial Pattern Matching (CPM 2024)},
  organization={Schloss Dagstuhl-Leibniz-Zentrum f{\"u}r Informatik}
}
```
