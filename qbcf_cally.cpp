#include "qbcf.h"
#include "qgenlib/tsv_reader.h"
#include "qgenlib/qgen_error.h"
#include "htslib/synced_bcf_reader.h"
#include "nuclear_pedigree.h"
#include "sex_ploidy_map.h"
#include "bcf_utils.h"
#include <ctime>
#include <set>
#include <map>
#include <cmath>
#include <sys/stat.h>
#include <sys/types.h>

#define MAX_LGAMMA_ARG 65536 // maximum value of lgamma(.) values to precalculate
#define AD_LIMIT 65000       // maximum number of allele depth

double softR(int32_t a, int32_t b, int32_t c, int32_t d, double offset) {
  double fa = a + offset;
  double fb = b + offset;
  double fc = c + offset;
  double fd = d + offset;
  double r = (fa*fd - fb*fc)/sqrt((fa+fb)*(fa+fc)*(fd+fb)*(fd+fc));
  return(r);
}

// A class that caches genotype calls given the pair of allele depths
class adCaller {
public:
  // adGeno is an internal class that represents the cached genotypes, it contains
  // geno_filt   : filtered genotype,   0 - MISSING, 1 - REF, 2 - ALT genotype (haploid)
  // geno_unfilt : unfiltered genotype, 0 - MISSING, 1 - REF, 2 - ALT genotype (haploid)
  // ad0 : allele depth of REF allele
  // ad1 : allele depth of ALT allele
  // lRef : likelihood of REF genotype
  // lAlt : likelihood of ALT genotype
  // lOth : likelihood of other(outlier) genotype
  class adGeno {
  public:
    int32_t geno_filt, geno_unfilt, ad0, ad1; // 0 : MISSING, 1 : REF, 2 : ALT
    double lRef, lAlt, lOth;
    adGeno(int32_t _geno_filt, int32_t _geno_unfilt, int32_t _ad0, int32_t _ad1, double _lRef, double _lAlt, double _lOth) :
      geno_filt(_geno_filt), geno_unfilt(_geno_unfilt), ad0(_ad0), ad1(_ad1), lRef(_lRef), lAlt(_lAlt), lOth(_lOth) {}
  };
  std::map<uint32_t, adGeno*> genoCache; // dictionary of AD pair to adGeno object : [16bit * 65536 + 16bit] -> yGeno* 

public:
  int32_t mode; // mode 0 : diploid, 1 : X, 2 : Y, 3: MT, 4...? - currently only support Y
  std::vector<adGeno*> genos; // list of genotypes for each individual
  int32_t alpha; // pseudocount parameter in Beta distribution for modeling outlier
  int32_t beta;  // pseudocount parameter in Beta distribution for modeling outlier
  double le0, le1, logdiff;  // le0 = log(1-error), le1 : log(error), logdiff : threshold for confident genotype call
  double lgammas[MAX_LGAMMA_ARG]; // stored values of lgamma
  double rDP, fMNZ, fFNZ, fMIS, fED, fOD; // parmaeters for variant filtering

  // perform chrY genotype calling based on a single pair of allele depths
  adGeno* ad2geno(int32_t ad0, int32_t ad1) {
    uint32_t adkey = (uint32_t)ad0 * AD_LIMIT + (uint32_t)ad1; // make a pair of AD into 
    std::map<uint32_t, adGeno*>::iterator it = genoCache.find(adkey);
    if ( it == genoCache.end() ) {  // update cache only when stored value is not found
      double lRef = le0 * ad0 + le1 * ad1;
      double lAlt = le1 * ad0 + le0 * ad1;
      double lOth = lgammas[ad0+1] + lgammas[ad1+1] + lgammas[ad0+alpha] + lgammas[ad1+beta] - lgammas[ad0+ad1+1] - lgammas[alpha] - lgammas[beta] - lgammas[ad0+alpha+ad1+beta];
      int32_t geno_filt = 0; // MISSING = 0
      if ( lRef > lAlt && lRef > lOth - logdiff ) geno_filt = 1;       // REF allele
      else if ( lAlt > lRef && lAlt > lOth - logdiff ) geno_filt = 2;  // ALT allele
      int32_t geno_unfilt = ad0 > ad1 ? 1 : (ad0 < ad1 ? 2 : 0);       // unfiltered genotypes follows the most frequent allele
      adGeno* p = new adGeno(geno_filt, geno_unfilt, ad0, ad1, lRef, lAlt, lOth);
      genoCache[adkey] = p;
      return p;
    }
    else {
      return it->second;
    }
  }

  // constructor - currently supports chrY only
  adCaller(int32_t _mode = 2, int32_t _alpha=5, int32_t _beta=5, double _error = 0.01, double _logdiff = 10) : mode(_mode), alpha(_alpha), beta(_beta), logdiff(_logdiff) {
    if ( mode != 2 ) { // only accepts chrY mode for now
      error("Cannot recognize the mode %d", _mode);
    }

    le1 = log(_error);
    le1 = log(1-_error);
    lgammas[0] = lgammas[1] = 0;
    for(int32_t i=2; i < MAX_LGAMMA_ARG; ++i) {
      lgammas[i] = lgammas[i-1] + log(i-1.0);
    }
    rDP = fMNZ = fFNZ = fMIS = fED = fOD = 0;
  }

  // destructor 
  ~adCaller() {
    std::map<uint32_t, adGeno*>::iterator it;
    for(it = genoCache.begin(); it != genoCache.end(); ++it) {
      delete it->second;
    }
  }

  // perform genotype calling
  int32_t callY(int32_t nsamps, int32_t* ads, int32_t* dps, int8_t* ploidies) {
    int32_t ad0, ad1, dp;
    // allocate genos as NULL vectors - which will be filled with genoCache
    if ( (int32_t)genos.size() != nsamps )
      genos.resize(nsamps, NULL);

    uint64_t sumDPs[2] = {0, 0}; // male, female depths
    uint64_t sumEDs[2] = {0, 0}; // sum of depths of 'incorrect' alleles matching with genotype, male and female separately
    uint64_t sumCDs[2] = {0, 0}; // sum of depths of 'correct'   alleles matching with genotype, male and female separately
    uint64_t sumODs[2] = {0, 0}; // sum of depths of 'other'     alleles not matching with either allele, male and female separately.
    int32_t  cntNZs[2] = {0, 0}; // male, female nonzero counts
    int32_t cntGenos[2][3] = { {0, 0, 0}, {0, 0, 0} };  // genotype counts, separated by sex
    int8_t pd;  // ploidy - temporary variable
    int32_t  nsexes[2] = {0, 0}; // number of individuals by sex
    
    for(int32_t i=0; i < nsamps; ++i) {
      ad0 = ads[2*i]   >= AD_LIMIT ? AD_LIMIT - 1 : ads[2*i];   // REF allele depth, capped by AD_LIMIT
      ad1 = ads[2*i+1] >= AD_LIMIT ? AD_LIMIT - 1 : ads[2*i+1]; // ALT allele depth, capped by AD_LIMIT
      dp  = dps[i] >= AD_LIMIT ? AD_LIMIT-1 : dps[i];           // total depth, capped by AD_LIMIT

      genos[i] = ad2geno(ad0, ad1); // call genotype (or use cached genotype) based on allele depth

      // TODO: perform ploidy-aware filtering
      pd = ploidies[i]; // pd is 0 for female, 1 for male
      ++nsexes[pd];     
      sumDPs[pd] += dps[i];       // update sumDP
      sumODs[pd] += ( dp > ad0 + ad1 ? dp - ad0 - ad1 : 0 ); // update sumED
      cntNZs[pd] += (dps[i] > 0); // number of nonzero 
      ++cntGenos[pd][genos[i]->geno_filt];

      // update sumCD and sumED only when 
      if ( genos[i]->geno_filt == 1 ) {  
        sumCDs[pd] += ad0;   
        sumEDs[pd] += ad1;
      } else if ( genos[i]->geno_filt == 2 ) {
        sumCDs[pd] += ad1;        
        sumEDs[pd] += ad0;
      }
    }

    // update filtering metrics
    rDP = ((sumDPs[0]+1e-10)/(nsexes[0]+2e-10)) / ((sumDPs[1]+1e-10)/(nsexes[1]+2e-10));
    fFNZ = (cntNZs[0]+1e-10)/(nsexes[0]+2e-10);
    fMNZ = (cntNZs[1]+1e-10)/(nsexes[1]+2e-10);    
    fMIS = (cntGenos[1][0]+1e-10) / (nsexes[1] + 2e-10);
    fED  = (sumEDs[1] + 1e-10)/(sumEDs[1] + sumCDs[1] + 2e-10);
    fOD  = (sumODs[1] + 1e-10)/(sumDPs[1] + 2e-10);

    //notice("%.5lf\t%.5lf\t%.5lf", rDP, fNZ, fMIS);
    return nsamps;
  }
};

// haploid genotype concordance between a pair of individuals
// 0 : MISSING, 1 : 
class HaploidGenoConcordance {
 public:
  std::map<std::string,int32_t> counts;  
  int32_t nDups;

  HaploidGenoConcordance(int32_t numDups = 0) : nDups(numDups) {}
  void clear() { counts.clear(); }

  void addGenotype(std::vector<int32_t>& gDups) {
    std::string s(gDups.size(),'0');
    for(int32_t i=0; i < (int32_t)gDups.size(); ++i) {
      int32_t g = gDups[i];
      if ( ( g < 0 ) || ( g > 2 ) )
        error("[E:%s:%d %s] Cannot recognize genotype %d - i=%d.. from HaploidGenoConcordance::addGenotype(%d..%d)",__FILE__,__LINE__,__FUNCTION__,g,i,gDups[0],gDups[i]);
      s[i] = '0' + g;
    }
    ++counts[s];    
  }
  void addGenotype(int32_t g1, int32_t g2) {
    std::vector<int32_t> gs(1,g1);
    gs.push_back(g2);
    addGenotype(gs);    
  }
  
  int32_t fillDupCount(int32_t index1, int32_t index2, std::vector<int32_t>& c9) {
    int32_t sum = 0;
    c9.resize(9,0);
    std::fill(c9.begin(), c9.end(), 0);
    std::map<std::string,int32_t>::iterator it;
    for(it = counts.begin(); it != counts.end(); ++it) {
      const char* p = it->first.c_str();
      c9[(p[index1] - '0') * 3 + (p[index2] - '0')] += it->second;
      sum += it->second;
    }
    return sum;
  }
};

// class to perform Mendelian concordance analysis on Y chromosome
class mendelY {
public:
  NuclearPedigree* ped;
  adCaller* pCaller;
  int32_t minDP;
  std::map<NuclearFamily*, HaploidGenoConcordance> famConc;       // family-level concordance metric
  std::map<NuclearFamilyPerson*, HaploidGenoConcordance> dupConc; // duplicate-level concordance metric
  HaploidGenoConcordance varFamConc;
  HaploidGenoConcordance varDupConc;  

  // constructor to intialize the object
  mendelY(NuclearPedigree* _ped, adCaller* _pCaller, int32_t _minDP = 0) : ped(_ped), pCaller(_pCaller), minDP(_minDP), varFamConc(2), varDupConc(2) {
    std::map<std::string, NuclearFamily*>::iterator itF;
    for(itF = ped->famIDmap.begin(); itF != ped->famIDmap.end(); ++itF) {
      NuclearFamily* pFam = itF->second; 
      if ( ( pFam->pDad ) && ( pFam->pDad->samples.size() > 1 ) ) // for males, make duplicate concordance table
        dupConc.insert( std::make_pair(pFam->pDad, HaploidGenoConcordance((int32_t)pFam->pDad->samples.size())) );
      //if ( ( pFam->pMom ) && ( pFam->pMom->samples.size() > 1 ) ) // do not make duplicate concordance table for females
      //  dupConc.insert( std::make_pair(pFam->pMom, DupConcordance((int32_t)pFam->pMom->samples.size())) );
    
      if ( !pFam->pKids.empty() ) {
        //famConc.insert( std::make_pair(pFam, FamilyConcordance((int32_t)pFam->pKids.size())) );
        famConc.insert( std::make_pair(pFam, HaploidGenoConcordance((int32_t)pFam->pKids.size()+1)) );
        for(int32_t i=0; i < (int32_t)pFam->pKids.size(); ++i) { // do something only for males
          if ( pFam->pKids[i]->sex == 1 ) { // focus on male offsprings
            if ( pFam->pKids[i]->samples.size() > 1 )
              dupConc.insert( std::make_pair(pFam->pKids[i], HaploidGenoConcordance((int32_t)pFam->pKids[i]->samples.size())) );
          }
        }
      }     
    }
  }

  // get a person's genotypes and depths (as sum of AD0 and AD1) from adGenos
  int32_t getPersonGenoDepth(std::vector<adCaller::adGeno*>& adGenos, NuclearFamilyPerson* pPerson, std::vector<int>& genos, std::vector<int>& depths) {
    genos.clear();
    depths.clear();
    if ( pPerson == NULL ) return 0;
    else {
      int32_t nSamples = 0;
      for(int32_t i=0; i < (int32_t)pPerson->samples.size(); ++i) {
        int32_t idx = pPerson->samples[i]->index;
        if ( idx >= 0 ) { // person found
          adCaller::adGeno* ag = adGenos[idx];
          genos.push_back(ag->geno_filt);
          depths.push_back(ag->ad0 + ag->ad1);
          ++nSamples;
        }
      }
      return nSamples;
    }
  }  

  // perform Mendelian concordance check for dups and families in chrY
  // chrY genotypes : males 0/1/2, females : NA
  // DupConc (males only) : 3 * 3 table - (0,0) and (1,1) are only concordant
  // FamConc (male kids only) : 3 * 3 table (with dad) - (0,0) and (1,1) are only concordant
  int32_t mendelCheckY() {
    std::map<std::string, NuclearFamily*>::iterator itF;

    std::vector<int32_t> dadGTs;
    std::vector<int32_t> nSons;
    std::vector< std::vector<int32_t> > sonGTs;
    std::vector<int32_t> dadDPs;
    std::vector< std::vector<int32_t> > sonDPs;
    std::vector<int32_t> famGTs;
    std::vector<int32_t> duoGTs(2);    

    int32_t i, j, k, m;

    varDupConc.clear(); // varDupConc stores duplicate concordance across all duplicate pairs for the variant
    varFamConc.clear(); // varFamConc stores duplicate concordance across all father-offspring pairs for the variant
    
    // iterate each family
    for(itF = ped->famIDmap.begin(); itF != ped->famIDmap.end(); ++itF) {
      NuclearFamily* pFam = itF->second;

      // calculate genotype concordance
      famGTs.clear(); // famGTs represents the genotypes of current (male) family members

      int32_t nDad = getPersonGenoDepth(pCaller->genos, pFam->pDad, dadGTs, dadDPs); // get the Dad's genotype
      for(i=0; i < nDad; ++i) { // for all genotype for the Dad (including duplicates)
        if ( dadDPs[i] < minDP ) dadGTs[i] = 0; // if the depth is less than the 
      }
      famGTs.push_back(nDad > 0 ? dadGTs[0] : 0);
      duoGTs[0] = (nDad > 0 ? dadGTs[0] : 0);
      
      if ( nDad > 1 ) {
        dupConc[pFam->pDad].addGenotype(dadGTs);
        for( i=1; i < nDad; ++i )
          for( k=0; k < i; ++k) 
            varDupConc.addGenotype(dadGTs[k],dadGTs[i]);
      }
      if ( !pFam->pKids.empty() ) { // the family 
        nSons.clear();
        sonGTs.clear();
        sonDPs.clear();
        for(j=0, k = 0; j < (int32_t)pFam->pKids.size(); ++j) {
          if ( pFam->pKids[j]->sex == 1 ) {
            nSons.resize(k+1);
            sonGTs.resize(k+1);
            sonDPs.resize(k+1);
            nSons[k] = getPersonGenoDepth( pCaller->genos, pFam->pKids[j], sonGTs[k], sonDPs[k]);
            for(i=0; i < nSons[k]; ++i) {
              if ( sonDPs[k][i] < minDP )
                sonGTs[k][i] = 0;
              duoGTs[1] = sonGTs[k][i];
              varFamConc.addGenotype(duoGTs);
            }
            famGTs.push_back( (nSons[k] > 0) ? sonGTs[k][0] : 0);            
            if ( nSons[k] > 1 ) {
              dupConc[pFam->pKids[j]].addGenotype(sonGTs[k]);
              for( i=1; i < nSons[k]; ++i ) {
                for( m=0; m < i; ++m) 
                  varDupConc.addGenotype(sonGTs[k][m],sonGTs[k][i]);
              }
            }
            k++;
          }
        }
        // get the duplicate concordance and trio concordance
        famConc[pFam].addGenotype(famGTs);
      }
    }
    return 0;
  }
};

/////////////////////////////////////////////////////////////////////////
// dge-barcode-match : Match a pair of barcode libraries
////////////////////////////////////////////////////////////////////////
int32_t qbcf_callY(int32_t argc, char** argv) {
  std::string inbcf;  
  std::string outbcf;
  std::string sexmapf;
  std::string pedf;
  std::string region;
  double alpha = 5;
  double beta = 5;

  int32_t xStart = 2781479;
  int32_t xStop = 15570138;
  std::string smID;
  std::string xLabel("chrX");
  std::string yLabel("chrY");
  std::string mtLabel("chrM");
  int32_t mendelMinDP = 0;

  paramList pl;

  BEGIN_LONG_PARAMS(longParameters)
    LONG_PARAM_GROUP("Input options", NULL)
    LONG_STRING_PARAM("bcf", &inbcf, "Input BCF/VCF file")
    LONG_STRING_PARAM("region", &region, "Region in [chrom]:[beg]-[end] format (1-based, inclusive)")
    LONG_STRING_PARAM("sex-map", &sexmapf, "Sex map file (1 male, 2 female) for all samples")
    LONG_STRING_PARAM("ped", &pedf, "PED file for nuclear families")
    LONG_STRING_PARAM("y-label", &yLabel, "Chromosome Y label")
    LONG_INT_PARAM("mendel-min-dp", &mendelMinDP, "Minimum depth for Mendelian inconsistency check")

    LONG_PARAM_GROUP("Model parameters", NULL)    
    LONG_DOUBLE_PARAM("alpha", &alpha, "Alpha parameter for beta-binomial distribution")
    LONG_DOUBLE_PARAM("beta", &beta, "Beta parameter for beta-binomial distribution")

    LONG_PARAM_GROUP("Output Options", NULL)
    LONG_STRING_PARAM("out",&outbcf,"Output file name")
  END_LONG_PARAMS();

  pl.Add(new longParams("Available Options", longParameters));
  pl.Read(argc, argv);
  pl.Status();

  // check input file sanity
  if ( inbcf.empty() || sexmapf.empty() || pedf.empty() || outbcf.empty() )
    error("--bcf, --sex-map, --ped, --out is required but at least one is missing");

  bcf_srs_t* sr = bcf_sr_init();

  if ( !region.empty() ) {
    if ( bcf_sr_set_regions(sr, region.c_str(), 0) < 0 )
      error("Cannot apply --region %s parameter", region.c_str());
  }  
  
  if ( !bcf_sr_add_reader(sr, inbcf.c_str()) )
    error("Failed to open %s", inbcf.c_str());

  // load the pedigree
  notice("Loading pedigree file %s", pedf.c_str());
  NuclearPedigree* ped = new NuclearPedigree(pedf.c_str());

  bcf_hdr_t* hdr = sr->readers[0].header;
  
  int32_t nsamples = bcf_hdr_nsamples(hdr);
  for(int32_t i=0; i < nsamples; ++i)
    ped->setSampleIndex(hdr->samples[i], i);
  notice("In the Original Pedigree: %d families, %d samples, %d unique individuals", (int32_t)ped->famIDmap.size(), (int32_t)ped->smIDmap.size(), ped->numPeople());
  
  int32_t nremoved = ped->removeSamplesWithoutIndex();
  notice("Removed %d samples not in the VCF file from the pedigree",nremoved);
  notice("Remaining samples : %d families, %d individuals, %d sequenced samples", (int32_t)ped->famIDmap.size(), ped->numPeople(), ped->numSamplesWithIndex());

  // load the sex map file
  sex_ploidy_map spmap(xLabel, yLabel, mtLabel, xStart, xStop);
  spmap.load_sex_map_file(sexmapf.empty() ? NULL : sexmapf.c_str(), hdr);
  notice("Loaded sex map of %zu samples, and identified %d males", spmap.vSex.size(), spmap.n_males);

  // iterate over the BCF file
  bcf_record rec(hdr);  
  int32_t nrecs = 0;
  int8_t* ploidies = NULL;

  // open VCF for writing and write headers
  const char* mode = "w";
  if ( outbcf.compare(outbcf.size() - 4, 4, ".bcf") == 0 ) mode = "wb";
  else if ( outbcf.compare(outbcf.size() - 7, 7, ".vcf.gz") == 0 ) mode = "wz";
  vcfFile* outVcf = bcf_open(outbcf.c_str(), mode);
  bcf_hdr_t* ohdr = bcf_hdr_init("w");
  bcf_hdr_transfer_contigs(hdr, ohdr); // transfer the contig
  
  const char* fmt_hdrs[3] = {
    "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">",
    "##FORMAT=<ID=AD,Number=A,Type=Integer,Description=\"Allele Depth\">",
    "##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Total Allele Depth, including unidentifiable alleles\">"
  };
  for(int i=0; i < 3; ++i) bcf_hdr_append(ohdr, fmt_hdrs[i]);
  
  const char* info_hdrs[13] = {
    "##INFO=<ID=AC,Number=A,Type=Integer,Description=\"Alternate Allele Counts\">",
    "##INFO=<ID=AN,Number=1,Type=Integer,Description=\"Total Number Allele Counts\">",
    "##INFO=<ID=AF,Number=A,Type=Float,Description=\"Alternate Allele Frequency from Best-guess Genotypes\">",
    "##INFO=<ID=RDP,Number=1,Type=Float,Description=\"Relative depth of females to males\">",
    "##INFO=<ID=MZD,Number=1,Type=Float,Description=\"Proportion of males with zero depth\">",
    "##INFO=<ID=FZD,Number=1,Type=Float,Description=\"Proportion of females with zero depth\">",
    "##INFO=<ID=MIS,Number=1,Type=Float,Description=\"Proportion of missing genotypes in males\">",
    "##INFO=<ID=FED,Number=1,Type=Float,Description=\"Proportion of putatively incorrect (switched) alleles\">",
    "##INFO=<ID=FOD,Number=1,Type=Float,Description=\"Proportion of alleles neither REF nor ALT alleles\">",
    "##INFO=<ID=DUPC,Number=4,Type=Integer,Description=\"Genotype count for each dup pair of males in the order of RR, RA, AR, AA\">",
    "##INFO=<ID=DUPR,Number=1,Type=Float,Description=\"Correlation coefficient for each dup pair with soft offset 0.01\">",    
    "##INFO=<ID=MENC,Number=4,Type=Integer,Description=\"Genotype count for each father-offspring pair in the order of RR, RA, AR, AA\">",
    "##INFO=<ID=MENR,Number=1,Type=Float,Description=\"Correlation coefficient for each father-offspring pair with soft offset 0.01\">"
  };
  for(int i=0; i < 13; ++i) bcf_hdr_append(ohdr, info_hdrs[i]);
  
  const char* filt_hdrs[5] = {
    "##FILTER=<ID=RDP005,Description=\"Relative depth (RDP) is less than .005\">",
    "##FILTER=<ID=FZD05,Description=\"Fraction of females with zero depth is greater than .05\">",
    "##FILTER=<ID=MENR95,Description=\"Mendelian correlation is lower than .95 with 2 or more Mendelian errors\">",
    "##FILTER=<ID=DUPR95,Description=\"Duplicate correlation is lower than .95 with 1 or more errors\">",    
    "##FILTER=<ID=MIS05,Description=\"Missingness in males is greater than .05\">"
  };
  for(int i=0; i < 5; ++i) bcf_hdr_append(ohdr, filt_hdrs[i]);

  for(int32_t i=0; i < nsamples; ++i) {
    if ( bcf_hdr_add_sample(ohdr, rec.get_sample_id(i).c_str()) < 0 )
      error("Cannot add sample ID %s", rec.get_sample_id(i).c_str());
  }  
  
  //if ( updateFMIS && ( bcf_hdr_id2int(odr->hdr, BCF_DT_ID, "FMIS" ) < 0 ) ) {
  //  sprintf(buffer,"##INFO=<ID=FMIS,Number=1,Type=Float,Description=\"Fraction of missing genotype\">\n");
  //  bcf_hdr_append(odw.hdr, buffer);
  //}
  if ( bcf_hdr_sync(ohdr) < 0 )
    error("Cannot sync BCF header for writing");
  if ( bcf_hdr_write(outVcf, ohdr) < 0 )
    error("Cannot write header to %s", outbcf.c_str());  

  adCaller caller;
  mendelY mendY(ped, &caller, mendelMinDP);
  std::vector<int32_t> c9;
  std::string c9str;
  while( bcf_sr_next_line(sr) ) {
    if ( ++nrecs % 1000 == 0 )
      notice("Processing %d variants at %s", nrecs, rec.get_var_ID().c_str());

    rec.set_variant(sr->readers[0].buffer[0]);
    ploidies = spmap.get_ploidies(rec.rec);
    
    rec.parse_allele_depths();
    rec.parse_total_depths();    
    
    int32_t geno_counts[6] = {0, 0, 0, 0, 0, 0};
    caller.callY(nsamples, rec.ads, rec.dps, ploidies);
    for(int32_t i=0; i < nsamples; ++i) {
      ++geno_counts[caller.genos[i]->geno_filt + ploidies[i]* 3];
    }

    /*
    if ( geno_counts[5] >= 5 ) {
      notice("%s\t%d %d %d\t%d %d %d\t%.5lf %.5lf %.5lf %.5lf %.5lf %.5lf", rec.get_var_ID().c_str(), geno_counts[0], geno_counts[1], geno_counts[2], geno_counts[3], geno_counts[4], geno_counts[5], caller.rDP, caller.fMNZ, caller.fFNZ, caller.fMIS, caller.fED, caller.fOD);

      mendY.mendelCheckY();
      mendY.varDupConc.fillDupCount(0, 1, c9);
      c9str.clear();
      cat_join_int32(c9str, c9, " ");
      notice("%s", c9str.c_str());
      mendY.varFamConc.fillDupCount(0, 1, c9);
      c9str.clear();
      cat_join_int32(c9str, c9, " ");
      notice("%s", c9str.c_str());
    }
    */

    // write a new variant into the output VCF
    bcf1_t* nv = bcf_init();
    bcf_clear(nv);
    bcf_set_rid(nv, rec.rec->rid);  // copy CHROM
    bcf_set_pos0(nv, rec.rec->pos); // copy POS
    bcf_update_alleles(ohdr, nv, (const char**) rec.rec->d.allele, rec.rec->n_allele); // copy REF and ALT
    bcf_set_n_sample(nv, nsamples);
    int32_t* gts = (int32_t*) calloc(2*nsamples, sizeof(int32_t));
    
    for(int32_t i=0; i < nsamples; ++i) {
      if ( ploidies[i] == 0 ) {
        gts[2*i] = gts[2*i+1] = bcf_gt_missing;
      }
      else if ( ploidies[i] == 1 ) {
        if ( caller.genos[i]->geno_filt == 0 ) {
          gts[2*i] = gts[2*i+1] = bcf_gt_missing;          
        }
        else {
          gts[2*i] = gts[2*i+1] = bcf_gt_unphased(caller.genos[i]->geno_filt-1);
        }
      }
    }
    bcf_update_format_int32(ohdr, nv, "GT", gts, nsamples * 2);
    bcf_update_format_int32(ohdr, nv, "AD", rec.ads, nsamples * 2);
    bcf_update_format_int32(ohdr, nv, "DP", rec.dps, nsamples);

    // update info field
    int32_t an = geno_counts[4]+geno_counts[5];
    int32_t ac = geno_counts[5];
    float rdp = (float)caller.rDP;
    float mzd = (float)caller.fMNZ;
    float fzd = (float)caller.fFNZ;
    float mis = (float)caller.fMIS;
    float fed = (float)caller.fED;
    float fod = (float)caller.fOD;
    
    bcf_update_info_int32(ohdr, nv, "AN", &an, 1);
    bcf_update_info_int32(ohdr, nv, "AC", &ac, 1);
    bcf_update_info_float(ohdr, nv, "RDP", &rdp, 1);
    bcf_update_info_float(ohdr, nv, "MZD", &mzd, 1);
    bcf_update_info_float(ohdr, nv, "FZD", &fzd, 1);
    bcf_update_info_float(ohdr, nv, "MIS", &mis, 1);
    bcf_update_info_float(ohdr, nv, "FED", &fed, 1);
    bcf_update_info_float(ohdr, nv, "FOD", &fod, 1);

    // update FILTER
    std::vector<std::string> flts;
    if ( rdp > 0.005 )  flts.push_back("RDP005");
    if ( fzd > 0.05 )   flts.push_back("FZD05");
    if ( mis > 0.05 )   flts.push_back("MIS05");

    // Evaluation of Mendelian errors
    mendY.mendelCheckY();
    mendY.varDupConc.fillDupCount(0, 1, c9);
    int32_t c4[4];
    c4[0] = c9[4]; c4[1] = c9[5]; c4[2] = c9[7]; c4[3] = c9[8];
    float dupDis = c4[1]+c4[2]; // 4 - R/R, 5 - R/A, 7 - A/R, 8 - A/A
    float dupR   = (float)softR(c9[4], c9[5], c9[7], c9[8], 0.01);
    bcf_update_info_int32(ohdr, nv, "DUPC", c4, 4);
    bcf_update_info_float(ohdr, nv, "DUPR", &dupR, 1);
    if ( dupR < 0.95 && dupDis > 0 ) {
      flts.push_back("DUPR95");
    }
    mendY.varFamConc.fillDupCount(0, 1, c9);
    c4[0] = c9[4]; c4[1] = c9[5]; c4[2] = c9[7]; c4[3] = c9[8];
    float famDis = c4[1]+c4[2]; // 4 - R/R, 5 - R/A, 7 - A/R, 8 - A/A
    float famR   = (float)softR(c9[4], c9[5], c9[7], c9[8], 0.01);
    bcf_update_info_int32(ohdr, nv, "MENC", c4, 4);
    bcf_update_info_float(ohdr, nv, "MENR", &famR, 1);
    if ( famR < 0.95 && famDis > 1 ) {
      flts.push_back("MENR95");
    }    

    if ( flts.empty() ) flts.push_back("PASS");
    int32_t n_flts = (int32_t)flts.size();
    int32_t flt_ids[15];
    for(int i=0; i < (int32_t)flts.size(); ++i) {
      flt_ids[i] = bcf_hdr_id2int(ohdr, BCF_DT_ID, flts[i].c_str());
    }
    bcf_update_filter(ohdr, nv, flt_ids, n_flts);    

    if ( bcf_write(outVcf, ohdr, nv) < 0 )
      error("Cannot write a variant at position %d", nv->pos+1);
    
    free(gts);
    bcf_destroy(nv);    
  }
  notice("Finished processing %d variants", nrecs);
  hts_close(outVcf);
  return 0;
}
