#ifndef NUCLEAR_PEDIGREE_H
#define NUCLEAR_PEDIGREE_H

#include <cstdlib>
#include <stdint.h>
#include <cstring>
#include <cmath>
#include <cfloat>
#include <vector>
#include <list>
#include <map>
#include <queue>

#include "qgenlib/hts_utils.h"


/**
 * A class for storing nuclear family allowing duplicates
 *
 */
class NuclearFamilySample;
class NuclearFamilyPerson;
class NuclearFamily;
class NuclearPedigree;

class NuclearFamilySample {  // Represents each sequenced/genotyped sample.
 public:
  std::string sampleID; // sample ID
  int index;      // sample index in the BCF file
  NuclearFamilySample(const char* id) : sampleID(id), index(-1) {};
};

class NuclearFamilyPerson {  // represents each individuals
 public:
  std::vector<NuclearFamilySample*> samples; // could be multiple samples if duplicated
  int sex;  // 0 : unknown, 1 : male, 2 : female
  NuclearFamilyPerson(int _sex) : sex(_sex) {}
  const char* getUID() { return samples.front()->sampleID.c_str(); }
  bool hasSample(const char* smID) {
    for(size_t i=0; i < samples.size(); ++i) {
      if ( samples[i]->sampleID.compare(smID) == 0 )
        return true;
    }
    return false;
  }
  void addSample(NuclearFamilySample* pSample) { samples.push_back(pSample); }
  int removeSamplesWithoutIndex(std::map<std::string,NuclearFamilySample*>& smIDmap) {
    int nRemoved = 0;
    for(int i=(int)samples.size()-1; i >= 0; --i) {
      if ( samples[i]->index < 0 ) {
	smIDmap.erase(samples[i]->sampleID);
	delete(samples[i]);               // call desctuctor
	samples.erase(samples.begin()+i); // remove the pointer from the vector
	++nRemoved;
      }
    }
    return nRemoved;
  }
};

class NuclearFamily {
 public:
  std::string famID;
  NuclearFamilyPerson* pDad;
  NuclearFamilyPerson* pMom;
  std::vector<NuclearFamilyPerson*> pKids;

  NuclearFamily(const char* fid) : famID(fid), pDad(NULL), pMom(NULL) {}
};

class NuclearPedigree
{
 public:
  std::map<std::string, NuclearFamily*> famIDmap;
  std::map<std::string, NuclearFamilySample*> smIDmap;

  NuclearPedigree(const char* pedFile);
  void addPerson(const char* famID, const std::vector<std::string>& smIDs, const char* dadID, const char* momID, int sex);
  bool setSampleIndex(const char* smID, int index);
  int numPeople();
  int numSamplesWithIndex();
  int removeSamplesWithoutIndex();
};

#endif
