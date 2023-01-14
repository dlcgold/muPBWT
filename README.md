# rePBWT
A PBWT-based light index  for biobank scale genotype data.

## Build
Build:
```
cmake -S . -B build 
cmake --build build
```
## Usage
```
cd build
./repbwt [options]
```
```
Usage: RLPBWT [options]

Options:
  -i, --input_file <path>	 macs file for panel
  -m, --memorize <path>	  path to save serialization 
  -l, --load <path>	 path to load serialization
  -o, --output <path>	 path to query output
  -q, --query <path>	 path to query file
  -e, --extend	 extend matches 
  -b, --bcf	use vcf/bcf as file format for both input and query file
  -v, --verbose	 extra prints
  -V, --fverbose	 extra prints for functions (cautions)
  -h, --help	 show this help message and exit
```

example using RLPBWT (slp mode) with lce queries and matches extended:

```
./repbwt -i <input matrix> -m <serialization file> -e -b
./repbwt -l <serialization file> -o <output file> -q <query file> -e -b```

