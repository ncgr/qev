#include <stdio.h>
#include <google/sparse_hash_map>
#define _cplusplus
#include "bam.h"
#include <list>
#include<ostream>

#ifndef BAM_COVERAGE
#define BAM_COVERAGE

using namespace std;
using google::sparse_hash_map;

using __gnu_cxx::hash;


typedef struct{
        bam1_t *first;
        bam1_t *second;
        unsigned int depth; //depth to draw pair in pileup.
} pair_t;

struct StrHash {
  size_t operator()(const string& s) const {
    return HASH_NAMESPACE::hash<const char*>()(s.c_str());
  }
};


typedef sparse_hash_map< const string,pair_t *, StrHash> phash_t;

class hash_list_t{
        public:
        list<pair_t *> *plist;
        phash_t *hash;
        list<pair_t *> *bad_plist;
        phash_t *bad_hash;
        ~hash_list_t();
        hash_list_t();
};

int fetch_func( const bam1_t *b, void *data);
void remove_singlets(list<pair_t *> *plist);

#endif


