#include "bcf_utils.h"
#include "qgenlib/hts_utils.h"
#include "qgenlib/qgen_error.h"
#include "qgenlib/phred_helper.h"

bcf_record::bcf_record(bcf_hdr_t* _hdr) : hdr(_hdr) {
  int32_t ns = bcf_hdr_nsamples(hdr);
  for(int32_t i=0; i < ns; ++i) {
    sm_icols.push_back(i);
    sm_ids.push_back(bcf_hdr_sample_id(hdr,i));
  }
  rec = NULL;
  n_gts = n_dps = n_ads = 0;
  gts = dps = ads = NULL;
}

bcf_record::~bcf_record() {
  if ( gts != NULL ) free(gts);
  if ( dps != NULL ) free(dps);
  if ( ads != NULL ) free(ads);
}

int32_t bcf_record::set_samples(std::set<std::string>& sample_ids) {
  sm_icols.clear();
  sm_ids.clear();

  std::set<std::string>::iterator it = sample_ids.begin();
  while( it != sample_ids.end() ) {
    int32_t idx = bcf_hdr_sample_index(hdr, it->c_str());
    if ( idx < 0 )
      error("ERROR: Cannot find sample ID %s from the BCF file", it->c_str());
    sm_icols.push_back(idx);
    sm_ids.push_back(bcf_hdr_sample_id(hdr, idx));
    ++it;
  }
  return get_nsamples();
}

bool bcf_record::set_variant(bcf1_t* v) {
  rec = v;
  n_gts = 0;
  n_ads = 0;
  n_dps = 0;
  return true;
}

bool bcf_record::parse_genotypes() {
  if ( n_gts > 0 )
    return false;

  bcf_unpack(rec, BCF_UN_ALL);
  if ( bcf_get_genotypes(hdr, rec, &gts, &n_gts) < 0 ) {
    notice("Failed to parse genotypes at %s:%d:%s:%s", bcf_hdr_id2name(hdr,rec->rid), rec->pos+1, rec->d.allele[0], rec->d.allele[1]);
    return false;
  }

  int32_t nsamples = bcf_hdr_nsamples(hdr);
  acs.resize(rec->n_allele);
  std::fill(acs.begin(), acs.end(), 0);
  an = 0;

  int32_t tmp_gt;
  for(int32_t i=0; i < get_nsamples(); ++i) {
    for(int32_t j=0; j < 2; ++j) { // need to be changed to account for ploidies
      tmp_gt = gts[sm_icols[i]*2+j];
      //if ( ( j < ploidies[i]) && ( bcf_gt_allele(tmp_gt) >= 0 ) ) {
      if ( bcf_gt_allele(tmp_gt) >= 0 ) {
        ++an;
        ++acs[bcf_gt_allele(tmp_gt)];
      }
    }
  }
  return true;
}

bool bcf_record::parse_total_depths(const char* name) {
  if ( n_dps > 0 )
    return false;

  bcf_unpack(rec, BCF_UN_ALL);
  if ( bcf_get_format_int32(hdr, rec, name, &dps, &n_dps) < 0 ) {
    notice("Failed to parse total depths at %s:%d:%s:%s", bcf_hdr_id2name(hdr,rec->rid), rec->pos+1, rec->d.allele[0], rec->d.allele[1]);
    return false;
  }

  int32_t nsamples = bcf_hdr_nsamples(hdr);
  if ( nsamples != n_dps )
    error("Numer of DP fields %d do not match with the sample size %d at %s:%d:%s:%s", n_dps, nsamples, bcf_hdr_id2name(hdr,rec->rid), rec->pos+1, rec->d.allele[0], rec->d.allele[1]);

  // handle missingness 
  for(int32_t i=0; i < n_dps; ++i) {
    if ( dps[i] < 0 ) {
      dps[i] = 0;
    }
  }
  return true;
}

bool bcf_record::parse_allele_depths(const char* name) {
  if ( n_dps > 0 )
    return false;

  bcf_unpack(rec, BCF_UN_ALL);
  if ( bcf_get_format_int32(hdr, rec, name, &ads, &n_ads) < 0 ) {
    notice("Failed to parse allele depths at %s:%d:%s:%s", bcf_hdr_id2name(hdr,rec->rid), rec->pos+1, rec->d.allele[0], rec->d.allele[1]);
    return false;
  }

  int32_t nsamples = bcf_hdr_nsamples(hdr);
  int32_t nalleles = rec->n_allele;
  if ( nsamples * nalleles != n_ads )
    error("Numer of AD fields %d do not match with the [sample size %d] * [number of alleles %d] at %s:%d:%s:%s", n_ads, nsamples, nalleles, bcf_hdr_id2name(hdr,rec->rid), rec->pos+1, rec->d.allele[0], rec->d.allele[1]);

  // handle missingness 
  for(int32_t i=0; i < n_ads; ++i) {
    if ( ads[i] < 0 )
      ads[i] = 0;
  }
  return true;
}

std::string& bcf_record::get_var_ID() {
  bcf_unpack(rec, BCF_UN_STR);
  varID.assign(bcf_get_chrom(hdr,rec));
  varID += ":";
  char buf[256];
  sprintf(buf,"%zd", rec->pos);
  varID += buf;
  for(int32_t i=0; i < rec->n_allele; ++i) {
    varID += ":";
    varID += rec->d.allele[i];
  }
  return varID;
}
