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
```
```
Usage: ./rePBWT [options]

Options:
  -i, --input_file <path>	 macs file for panel
  -s, --save <path>	  path to save serialization 
  -l, --load <path>	 path to load serialization
  -o, --output <path>	 path to query output
  -q, --query <path>	 path to query file
  -m, --macs	use macs as file format for both input and query file
  -v, --verbose	 extra prints
  -V, --fverbose	 extra prints for functions (cautions)
  -h, --help	 show this help message and exit
```

example using rePBWT (slp mode) with lce queries and matches extended:

```
./repbwt -i <input matrix> -s <serialization file>
./repbwt -l <serialization file> -o <output file> -q <query file>```

