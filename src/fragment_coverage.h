#include <stdio.h>
#include <google/sparse_hash_map>
#define _cplusplus
#include "bam.h"
#include <list>


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


//function to calculate the span of a pair.
void calc_span(unsigned int *start, unsigned int *stop, pair_t *pair);

//call back for fetching pairs
int fetch_func( const bam1_t *b, void *data);
void remove_singlets(list<pair_t *> *plist);

/*Function to calculate the coverage from a list of pairs.  Also returns mean_insert, read pair counts.
	@plist is a stl list contain pair objects to be tallied.
	@pcov 
	@main_inert returns the mean of the pair lengths
	@counts returns the number of pairs used to calculate coverage
	@length input for the length of the contig
*/
//void calculate_fragment_coverage(list<pair_t *> *plist, unsigned long *pcov, double *mean_insert, unsigned long * counts, unsigned int length);
void calculate_fragment_coverage(list<pair_t *> *plist, unsigned long *pcov, double *mean_insert, unsigned long * counts, unsigned long *fragment_hist , unsigned int length);

void calculate_fragment_hist(list<pair_t *> *plist, unsigned long *fragment_hist, unsigned int length);
void calculate_fragment_coverage(list<pair_t *> *plist, unsigned long *pcov, double *mean_insert, unsigned long * counts, unsigned int length);

