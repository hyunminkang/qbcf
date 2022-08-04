#ifndef __BCF_UTILS_H
#define __BCF_UTILS_H

#include <cstdlib>
#include <cstring>
#include <string>
#include <vector>
#include <set>
#include "htslib/vcf.h"
#include "qgenlib/hts_utils.h"

class bcf_record {
 public:
  // pointers to htslib objects
  bcf_hdr_t* hdr;
  bcf1_t* rec;

  // site-level information
  int32_t an;
  std::vector<double> acs;
  std::string varID;    

  // sample-level information
  std::vector<std::string> sm_ids;
  std::vector<int32_t> sm_icols;
  int32_t* gts;  // GT
  int32_t n_gts; 
  int32_t* dps;  // DP
  int32_t n_dps;  
  int32_t* ads;  // AD
  int32_t n_ads;
  //std::vector<uint8_t> ploidies;

  bcf_record(bcf_hdr_t* _hdr);
  ~bcf_record();
  int32_t set_samples(std::set<std::string>& sample_ids);
  inline int32_t get_nsamples() { return(int32_t)sm_icols.size(); }
  inline const std::string& get_sample_id(int32_t idx) { return sm_ids[sm_icols[idx]]; }
  bool set_variant(bcf1_t* v);

  // method to parse the fields from VCF
  bool parse_genotypes();
  bool parse_total_depths(const char* name = "DP");
  bool parse_allele_depths(const char* name = "AD");

  std::string& get_var_ID();  
};

#endif // __BCF_UTILS_H
