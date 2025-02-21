#include "htslib/vcf.h"
#include "include/rlpbwt_int.h"
#include <getopt.h>
#include <iostream>
#include <sdsl/sd_vector.hpp>
#include <set>
#include <unordered_set>
#include <vector>
void printHelp() {
  std::cout << "Usage: mupbwt [options]\n" << std::endl;
  std::cout << "Options:" << std::endl;
  std::cout << "  -i, --input_file <path>\t macs file for panel" << std::endl;
  std::cout << "  -s, --save <path>\t  path to save serialization "
            << std::endl;
  std::cout << "  -l, --load <path>\t path to load serialization" << std::endl;
  std::cout << "  -o, --output <path>\t path to query output" << std::endl;
  std::cout << "  -q, --query <path>\t path to macs query file" << std::endl;
  //    std::cout << "  -e, --extend\t extend matches"
  //              << std::endl;
  std::cout << "  -m, --macs\t use macs format for panel and queries"
            << std::endl;
  std::cout << "  -d, --details\t print memory usage details" << std::endl;
  std::cout << "  -v, --verbose\t extra prints" << std::endl;
  std::cout << "  -h, --help\t show this help message and exit" << std::endl;
}
std::vector<unsigned int> getIntersection(const std::set<unsigned int> &set1,
                                          const std::set<unsigned int> &set2) {
  std::vector<unsigned int> intersection;

  std::set_intersection(set1.begin(), set1.end(), set2.begin(), set2.end(),
                        std::back_inserter(intersection));

  return intersection;
}

std::tuple<std::vector<sdsl::sd_vector<>>, std::vector<std::string>,
           std::vector<unsigned int>>
build_ref(std::string filename) {

  std::vector<sdsl::bit_vector> panel_tmp;
  std::vector<sdsl::sd_vector<>> panel;
  std::vector<std::string> samples;
  std::vector<unsigned int> sites;

  htsFile *fp = hts_open(filename.c_str(), "rb");
  std::cout << "Reading reference VCF file...\n";
  if (fp == NULL) {
    throw FileNotFoundException{};
  }

  bcf_hdr_t *hdr = bcf_hdr_read(fp);
  bcf1_t *rec = bcf_init();
  auto height = bcf_hdr_nsamples(hdr) * 2;
  auto width = 0;
  for (int i = 0; i < bcf_hdr_nsamples(hdr); i++) {
    samples.push_back(std::string(hdr->samples[i]));
  }
  fp = hts_open(filename.c_str(), "rb");
  hdr = bcf_hdr_read(fp);
  rec = bcf_init();

  // this->cols = std::vector<rl_column>(this->width + 1);

  unsigned int count = 0;
  std::string last_col;
  std::string last_column;
  std::string new_column;

  // iterate each vcf record
  while (bcf_read(fp, hdr, rec) >= 0) {
    // if (auto search = ps.find(rec->pos); search == ps.end())
    //   continue;
    width++;
    sites.push_back(rec->pos);
    // std::cout << rec->pos << "\n";
    new_column = "";
    sdsl::bit_vector tmp(height);
    int bi = 0;
    bcf_unpack(rec, BCF_UN_ALL);
    // read SAMPLE
    int32_t *gt_arr = NULL, ngt_arr = 0;
    int i, j, ngt, nsmpl = bcf_hdr_nsamples(hdr);
    ngt = bcf_get_genotypes(hdr, rec, &gt_arr, &ngt_arr);
    int max_ploidy = ngt / nsmpl;
    for (i = 0; i < nsmpl; i++) {
      int32_t *ptr = gt_arr + i * max_ploidy;
      for (j = 0; j < max_ploidy; j++) {
        // if true, the sample has smaller ploidy
        if (ptr[j] == bcf_int32_vector_end)
          break;

        // missing allele
        if (bcf_gt_is_missing(ptr[j]))
          exit(-1);

        // the VCF 0-based allele index
        int allele_index = bcf_gt_allele(ptr[j]);
        tmp[bi] = allele_index;
        bi++;
      }
    }
    panel_tmp.push_back(tmp);
    free(gt_arr);
  }
  panel = tb(panel_tmp);

  bcf_hdr_destroy(hdr);
  hts_close(fp);
  bcf_destroy(rec);

  return {panel, samples, sites};
}
std::vector<std::string> extractAltRef(char *als) {
  std::vector<std::string> res;
  std::string str = "";
  int i = 0;
  int j = 0;
  while (j != 2) {
    if (als[i] == '\0') {
      res.push_back(str);
      str.clear();
      j++;
    } else {
      str += als[i];
    }
    i++;
  }
  return res;
}
std::vector<std::vector<std::string>> get_ref_info(std::string reference) {
  std::vector<std::vector<std::string>> data;
  htsFile *fp = hts_open(reference.c_str(), "rb");
  bcf_hdr_t *hdr = bcf_hdr_read(fp);
  bcf1_t *rec = bcf_init();
  // bcf_unpack(rec, BCF_UN_STR);
  while (bcf_read(fp, hdr, rec) >= 0) {
    bcf_unpack(rec, BCF_UN_ALL);

    std::vector<std::string> tmp;
    tmp.push_back(std::string(bcf_hdr_id2name(hdr, rec->rid)));
    tmp.push_back(std::to_string(rec->pos + 1));
    tmp.push_back(std::string(rec->d.id ? rec->d.id : "."));
    // tmp.push_back(std::string(rec->d.allele[0]));
    // char alt_buffer[1024] = "";
    // for (int i = 1; i < rec->n_allele; i++) {
    //   std::strcat(alt_buffer, rec->d.allele[i]);
    //   if (i < rec->n_allele - 1)
    //     std::strcat(alt_buffer, ",");
    // }
    // tmp.push_back(std::string(alt_buffer));

    std::vector<std::string> altRef = extractAltRef(rec->d.als);
    tmp.push_back(altRef[0]);
    tmp.push_back(altRef[1]);
    tmp.push_back(std::to_string(rec->qual));
    // filter
    tmp.push_back(".");
    // info
    tmp.push_back(".");
    tmp.push_back("GT");
    // std::cout << tmp[0] << " " << tmp[1] << " " << tmp[2] << tmp[3] << " "
    //           << tmp[4] << " " << tmp[5] << "\n";
    data.push_back(tmp);
  }
  // std::cout << ps.size() << "\n";
  bcf_hdr_destroy(hdr);
  hts_close(fp);
  bcf_destroy(rec);
  return data;
}

int main(int argc, char **argv) {
  // TODO option to print stats
  if (argc == 1) {
    printHelp();
    exit(EXIT_SUCCESS);
  }
  bool verbose = false;
  bool details = false;
  bool query = false;
  bool macs = false;
  std::string matrix_input = "";
  // std::string memorize_file = "";
  // std::string load_file = "";
  std::string output = "";
  std::string reference = "";
  std::string query_input = "";
  bool unsafe = false;
  int c;

  // while (true) {
  //   static struct option long_options[] = {
  //       {"input", required_argument, nullptr, 'i'},
  //       {"save", required_argument, nullptr, 's'},
  //       {"load", required_argument, nullptr, 'l'},
  //       {"output", required_argument, nullptr, 'o'},
  //       {"query", required_argument, nullptr, 'q'},
  //       //{"extend",   no_argument,       nullptr, 'e'},
  //       {"macs", no_argument, nullptr, 'm'},
  //       {"details", no_argument, nullptr, 'd'},
  //       {"verbose", no_argument, nullptr, 'v'},
  //       {"help", no_argument, nullptr, 'h'},
  //       {nullptr, 0, nullptr, 0}};

  //   int option_index = 0;
  //   c = getopt_long(argc, argv, "i:s:l:o:q:emvdh", long_options,
  //   &option_index);

  //   if (c == -1) {
  //     break;
  //   }

  //   switch (c) {
  //   case 'i':
  //     matrix_input = optarg;
  //     break;
  //   case 's':
  //     memorize_file = optarg;
  //     break;
  //   case 'l':
  //     load_file = optarg;
  //     break;
  //   case 'o':
  //     output = optarg;
  //     break;
  //   case 'q':
  //     query_input = optarg;
  //     break;
  //   case 'm':
  //     macs = true;
  //     break;
  //   case 'd':
  //     details = true;
  //     break;
  //   case 'v':
  //     verbose = true;
  //     break;
  //   case 'h':
  //     printHelp();
  //     exit(EXIT_SUCCESS);
  //   default:
  //     printHelp();
  //     exit(EXIT_FAILURE);
  //   }
  // }

  while (true) {
    static struct option long_options[] = {
        {"input", required_argument, nullptr, 'i'},
        {"ref", required_argument, nullptr, 'r'},
        {"output", required_argument, nullptr, 'o'},
        {"query", required_argument, nullptr, 'q'},
        //{"extend",   no_argument,       nullptr, 'e'},
        {"details", no_argument, nullptr, 'd'},
        {"verbose", no_argument, nullptr, 'v'},
        {"unsafe", no_argument, nullptr, 'u'},
        {"help", no_argument, nullptr, 'h'},
        {nullptr, 0, nullptr, 0}};

    int option_index = 0;
    c = getopt_long(argc, argv, "i:o:r:q:evudh", long_options, &option_index);

    if (c == -1) {
      break;
    }

    switch (c) {
    case 'i':
      matrix_input = optarg;
      break;
    case 'o':
      output = optarg;
      break;
    case 'q':
      query_input = optarg;
      break;
    case 'r':
      reference = optarg;
      break;
    case 'd':
      details = true;
      break;
    case 'v':
      verbose = true;
      break;
    case 'u':
      unsafe = true;
      break;
    case 'h':
      printHelp();
      exit(EXIT_SUCCESS);
    default:
      printHelp();
      exit(EXIT_FAILURE);
    }
  }

  if (matrix_input.empty()) {
    std::cerr << "Error: input required\n";
    exit(EXIT_FAILURE);
  }
  if (reference.empty()) {
    std::cerr << "Error: reference required\n";
    exit(EXIT_FAILURE);
  }
  if (output.empty()) {
    std::cerr << "Error: output file required\n";
    exit(EXIT_FAILURE);
  }
  if (query_input.empty()) {
    std::cerr << "Error: query file required\n";
    exit(EXIT_FAILURE);
  }

  auto data_ref = build_ref(reference);

  htsFile *fp1 = hts_open(query_input.c_str(), "rb");
  std::cout << "Reading VCF file...\n";
  if (fp1 == NULL) {
    throw FileNotFoundException{};
  }

  bcf_hdr_t *hdr1 = bcf_hdr_read(fp1);
  bcf1_t *rec1 = bcf_init();
  std::unordered_set<unsigned int> ps;
  std::vector<std::string> samples;
  std::vector<std::string> samples_r;
  std::vector<std::vector<std::string>> data;
  for (int i = 0; i < bcf_hdr_nsamples(hdr1); i++) {
    samples.push_back(std::string(hdr1->samples[i]));
  }
  bcf_hdr_destroy(hdr1);
  hts_close(fp1);
  bcf_destroy(rec1);

  htsFile *fp2 = hts_open(matrix_input.c_str(), "rb");
  bcf_hdr_t *hdr2 = bcf_hdr_read(fp2);
  bcf1_t *rec2 = bcf_init();
  for (int i = 0; i < bcf_hdr_nsamples(hdr2); i++) {
    samples_r.push_back(std::string(hdr2->samples[i]));
  }
  std::unordered_set<unsigned int> ps2;
  std::vector<std::string> tmp;
  std::vector<unsigned int> s;
  while (bcf_read(fp2, hdr2, rec2) >= 0) {
    bcf_unpack(rec2, BCF_UN_ALL);
    ps2.insert(rec2->pos);
    s.push_back(rec2->pos);
  }
  // std::cout << ps2.size() << "\n";
  bcf_hdr_destroy(hdr2);
  hts_close(fp2);
  bcf_destroy(rec2);

  htsFile *fp = hts_open(query_input.c_str(), "rb");
  bcf_hdr_t *hdr = bcf_hdr_read(fp);
  bcf1_t *rec = bcf_init();
  // bcf_unpack(rec, BCF_UN_STR);
  while (bcf_read(fp, hdr, rec) >= 0) {
    bcf_unpack(rec, BCF_UN_ALL);
    if (ps2.find(rec->pos) == ps2.end())
      continue;

    ps.insert(rec->pos);
    std::vector<std::string> tmp;
    tmp.push_back(std::string(bcf_hdr_id2name(hdr, rec->rid)));
    tmp.push_back(std::to_string(rec->pos + 1));
    tmp.push_back(std::string(rec->d.id ? rec->d.id : "."));

    // tmp.push_back(std::string(rec->d.allele[0]));
    // char alt_buffer[1024] = "";
    // for (int i = 1; i < rec->n_allele; i++) {
    //   std::strcat(alt_buffer, rec->d.allele[i]);
    //   if (i < rec->n_allele - 1)
    //     std::strcat(alt_buffer, ",");
    // }
    // tmp.push_back(std::string(alt_buffer));

    std::vector<std::string> altRef = extractAltRef(rec->d.als);
    tmp.push_back(altRef[0]);
    tmp.push_back(altRef[1]);
    tmp.push_back(std::to_string(rec->qual));
    // filter
    tmp.push_back(".");
    // info
    tmp.push_back(".");
    tmp.push_back("GT");
    // std::cout << tmp[0] << " " << tmp[1] << " " << tmp[2] << tmp[3] << " "
    //           << tmp[4] << " " << tmp[5] << "\n";
    data.push_back(tmp);
  }
  // std::cout << ps.size() << "\n";
  bcf_hdr_destroy(hdr);
  hts_close(fp);
  bcf_destroy(rec);

  for (auto d : data) {
    if (auto search = ps.find(stoi(d[1]) - 1); search == ps.end())
      std::cout << d[1] << "\n";
  }

  // std::cout << "Requested " << ps.size() << " sites\n";
  clock_t START = clock();
  rlpbwt_int rlpbwt(matrix_input.c_str(), ps, verbose, 2);
  auto time_build = (float)(clock() - START) / CLOCKS_PER_SEC;
  std::cout << "Reference index built in: " << time_build << " s\n";
  START = clock();
  rlpbwt.query_match(matrix_input.c_str(), query_input.c_str(), output.c_str(),
                     ps, ps2, samples, samples_r, s, data,
                     get_ref_info(reference), data_ref, unsafe, verbose);
  time_build = (float)(clock() - START) / CLOCKS_PER_SEC;
  std::cout << "Unphased panel phased in: " << time_build << " s\n";
  // if (!query_input.empty()) {
  //   query = true;
  // }
  // if (matrix_input.empty() && load_file.empty()) {
  //   std::cerr << "Error: input or load file required\n";
  //   exit(EXIT_FAILURE);
  // }
  // if (memorize_file.empty() && (query_input.empty() || output.empty())) {
  //   if (!details) {
  //     std::cerr << "Error: nothing to do\n";
  //     printHelp();
  //     exit(EXIT_FAILURE);
  //   }
  // }
  // if (query && output.empty()) {
  //   std::cerr << "Error: output file required if query requested\n";
  //   printHelp();
  //   exit(EXIT_FAILURE);
  // }

  // if (query) {
  //   if (output.empty()) {
  //     std::cerr << "Error: output file required\n";
  //     exit(EXIT_FAILURE);
  //   }
  //   if (query_input.empty()) {
  //     std::cerr << "Error: query file required\n";
  //     exit(EXIT_FAILURE);
  //   }
  // }

  // clock_t START = clock();
  // if (load_file.empty()) {

  //   rlpbwt_int rlpbwt(matrix_input.c_str(), verbose, macs, 2);
  //   if (!memorize_file.empty()) {
  //     std::ofstream outstream;
  //     outstream.open(memorize_file.c_str());
  //     rlpbwt.serialize(outstream);
  //     outstream.close();
  //   }

  //   auto time_build = (float)(clock() - START) / CLOCKS_PER_SEC;
  //   std::cout << "built/loaded in: " << time_build << " s\n";
  //   if (details) {
  //     auto runs = rlpbwt.get_run_number();
  //     std::cout << "\n----\nTotal haplotypes: " << rlpbwt.height << "\n";
  //     std::cout << "Total sites: " << rlpbwt.width << "\n";
  //     std::cout << "----\nTotal runs: " << runs << "\n----\n";
  //     std::cout << "Average runs: " << std::ceil(runs / rlpbwt.width) <<
  //     "\n"; auto s = rlpbwt.size_in_mega_bytes(true); std::cout << "rlpbwt:
  //     "
  //     << s << " megabytes\n----\n"; std::cout << "estimated dense size: "
  //               << dense_size_megabyte(rlpbwt.height, rlpbwt.width)
  //               << " megabytes\n----\n";
  //   }
  //   if (query) {
  //     std::cout << "start querying \n";
  //     START = clock();
  //     rlpbwt.query_match(query_input.c_str(), output.c_str(), verbose,
  //     macs);

  //     auto time_query = (float)(clock() - START) / CLOCKS_PER_SEC;
  //     std::cout << "queried in " << time_query << " s\n";
  //   }
  // } else {
  //   rlpbwt_int rlpbwt;
  //   std::ifstream load;
  //   load.open(load_file.c_str());
  //   rlpbwt.load(load);
  //   load.close();
  //   std::cerr << "k=" << rlpbwt.k_smem << "\n";
  //   auto time_build = (float)(clock() - START) / CLOCKS_PER_SEC;

  //   std::cout << "built/loaded in: " << time_build << " s\n";
  //   if (details) {
  //     auto runs = rlpbwt.get_run_number();
  //     std::cout << "\n----\nTotal haplotypes: " << rlpbwt.height << "\n";
  //     std::cout << "Total sites: " << rlpbwt.width << "\n";
  //     std::cout << "----\nTotal runs: " << runs << "\n";
  //     std::cout << "Average runs: " << std::ceil(runs / rlpbwt.width)
  //               << "\n----\n";
  //     auto s = rlpbwt.size_in_mega_bytes(true);
  //     std::cout << "rlpbwt: " << s << " megabytes\n----\n";
  //     std::cout << "estimated dense size: "
  //               << dense_size_megabyte(rlpbwt.height, rlpbwt.width)
  //               << " megabytes\n----\n";
  //   }
  //   if (query) {
  //     std::cout << "start querying \n";
  //     START = clock();
  //     rlpbwt.query_match(query_input.c_str(), output.c_str(), verbose,
  //     macs);

  //     auto time_query = (float)(clock() - START) / CLOCKS_PER_SEC;
  //     std::cout << "queried in " << time_query << " s\n";
  //   }
  // }

  return 0;
}
