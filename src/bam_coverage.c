#include <stdio.h>
#include <google/sparse_hash_map>
#define _cplusplus
#include "bam.h"
#include <list>
#include<ostream>
#include "bam_coverage.h"

hash_list_t::hash_list_t(){
        plist = new list<pair_t *>();
        bad_plist = new list<pair_t *>();
        hash = new phash_t();
        bad_hash = new phash_t();
}

hash_list_t::~hash_list_t(){
        pair_t *pair;
        for(  phash_t::iterator i=hash->begin(); i != hash->end(); ++i){
                pair = i->second;
                bam_destroy1(pair->first);
                bam_destroy1(pair->second);
                free(pair);
        }
       for(  phash_t::iterator i=bad_hash->begin(); i != bad_hash->end(); ++i){
                pair = i->second;
                bam_destroy1(pair->first);
                bam_destroy1(pair->second);
                free(pair);
        }
        delete plist;
        delete bad_plist;
        delete hash;
        delete bad_hash;
}
void remove_singlets(list<pair_t *> *plist){
         //remove singlets
        pair_t *pair;
        for( list<pair_t *>::iterator i=plist->begin(); i != plist->end(); ++i){
                pair= *i;
                if(!(pair->first && pair->second)){
                        plist->erase(i);
                        i--;
                }
        }
}

int fetch_func( const bam1_t *b, void *data)
{
        uint32_t end;
        uint32_t end2;
        hash_list_t *hash_list = (hash_list_t *)data;
        phash_t *buf = hash_list->hash;
        list<pair_t *> *plist = hash_list->plist;
        phash_t *bad_buf = hash_list->bad_hash;
        list<pair_t *> *bad_plist = hash_list->bad_plist;


        pair_t *pair=NULL;
        //calculate the end of the alignment from b
        end =  bam_calend(&b->core,  bam1_cigar(b));
        end2 = b->core.isize + b->core.pos - 1;
        if(b->core.flag & BAM_FPROPER_PAIR){
                //printf("proper pair qname %s pos1 %d end1 %d pos2 %d end2 %d isize %d \n", bam1_qname(b), b->core.pos, end, b->core.mpos, end2, b->core.isize);
                string name = bam1_qname(b);
                //check to see if a read has been seen
                if(buf->find(name)==buf->end()){
                        //add alignment
                        pair=new pair_t();
                        pair->first=NULL;
                        pair->second=NULL;
                        (*buf)[name]=pair;
                        (*plist).push_back(pair);

                }
                if(b->core.flag & BAM_FREAD1){
                        (*buf)[name]->first = new bam1_t();
                        (*buf)[name]->first = bam_copy1((*buf)[name]->first, b);
                }else{
                        (*buf)[name]->second = new bam1_t();
                        (*buf)[name]->second = bam_copy1((*buf)[name]->second, b);
                }
        }else if(b->core.flag & BAM_FPAIRED && ! (b->core.flag & BAM_FUNMAP) && ! (b->core.flag & BAM_FMUNMAP) && !(b->core.flag &  BAM_FSECONDARY) ){
                // not proper pair, but both pairs are mapped and alignment is not a secondary alignment.
                       //printf("proper pair qname %s pos1 %d end1 %d pos2 %d end2 %d isize %d \n", bam1_qname(b), b->core.pos, end, b->core.mpos, end2, b->core.isize);
                string name = bam1_qname(b);
                //check to see if a read has been seen
                if(bad_buf->find(name)==bad_buf->end()){
                        //add alignment
                        pair=new pair_t();
                        pair->first=NULL;
                        pair->second=NULL;
                        (*bad_buf)[name]=pair;
                        (*bad_plist).push_back(pair);

                }
                if(b->core.flag & BAM_FREAD1){
                        (*bad_buf)[name]->first = new bam1_t();
                        (*bad_buf)[name]->first = bam_copy1((*bad_buf)[name]->first, b);
                }else{
                        (*bad_buf)[name]->second = new bam1_t();
                        (*bad_buf)[name]->second = bam_copy1((*bad_buf)[name]->second, b);
                }

        }
        return 0;
}


