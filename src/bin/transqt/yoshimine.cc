/*
 *@BEGIN LICENSE
 *
 * PSI4: an ab initio quantum chemistry software package
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along
 * with this program; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 *
 *@END LICENSE
 */

/*! \file
    \ingroup TRANSQT
    \brief Enter brief description of file here
*/
/*
** YOSHIMINE.C: Functions for the Yoshimine Sort Object
**
** David Sherrill
** Center for Computational Quantum Chemistry, UGA
** February 1995
**
** Made slightly more general to handle MP2 restricted transformations by
** Daniel Crawford
** September 1995
**
*/

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <cstring>
#include <libciomr/libciomr.h>
#include <libqt/qt.h>
#include <libiwl/iwl.h>
#include <libpsio/psio.hpp>
#include <libpsio/psio.h>
#include <psifiles.h>
#define YEXTERN
#include "yoshimine.h"
#include "MOInfo.h"

#include <boost/shared_ptr.hpp>

namespace psi { namespace transqt {

extern struct MOInfo moinfo;

#define MIN0(a,b) (((a)<(b)) ? (a) : (b))
#define MAX0(a,b) (((a)>(b)) ? (a) : (b))
#define INDEX(i,j) ((i>j) ? (ioff[(i)]+(j)) : (ioff[(j)]+(i)))

/*
** YOSH_INIT(): This function initializes a Yoshimine sort object.
**    The data is contained in a structure YBuff.
**
** Parameters:
**    YBuff          =  pointer to struct which will hold data for the object
**    bra_indices    =  the number of bra_index pairs for the two-electron
**                      integrals being sorted
**    ket_indices    =  the number of ket_index pairs for the two-electron
**                      integrals being sorted
**    maxcor         =  the core memory, in bytes
**    maxcord        =  the core memory, in doubles
**    max_buckets    =  the max number of buckets to use (may be limited due
**                      to the fact that each bucket requires a consecutively
**                      numbered binary temp file).
**    first_tmp_file =  the number of the first tmp file used in the
**                      Yoshimine sort (e.g. 80 for file80).
**    cutoff         =  minimum value to be kept for any value during the sort
**    outfile        =  the text output file
**
** Returns: none
**
** Note:  bra_indices and ket_indices replace nbstri in an attempt to somewhat
** generalize the sort for four-index quantities which index pairs may or
** may not be canonicalizable, e.g. integrals of (ov|ov) type, as may be
** found in MP2 energy calculations...the first and second (or third and
** fourth) indices are not necessarily interchangeable.  Additionally,
** this modification has the added benefit that it will work when the
** number of left and right basis functions are not equal (e.g., when
** one has half-backtransformed integrals of the type (ao ao | mo mo)
** where nbfao != nbfmo.
*/
void yosh_init(struct yoshimine *YBuff, unsigned bra_indices,
               unsigned ket_indices, long maxcor,
               long maxcord, const int max_buckets,
               unsigned int first_tmp_file,
               double cutoff, std::string OutFileRMR)
{
   unsigned long long int twoel_array_size;              /*--- Although on 32-bit systems one can only allocate 2GB arrays
                                                           in 1 process space, one can store much bigger integrals files on disk ---*/
   unsigned int nbuckets;
   int i, j, pq;
   unsigned long int bytes_per_bucket, free_bytes_per_bucket;

   YBuff->first_tmp_file = first_tmp_file;
   twoel_array_size = bra_indices; twoel_array_size *= ket_indices;
   YBuff->core_loads = (twoel_array_size - 1) / maxcord + 1 ;
   nbuckets = YBuff->core_loads ;
   outfile->Printf("nbuckets is %i\n", nbuckets );
   YBuff->nbuckets = nbuckets;
   YBuff->cutoff = cutoff;
   YBuff->bra_indices = bra_indices;
   YBuff->ket_indices = ket_indices;
   if (nbuckets > max_buckets) {
      outfile->Printf( "(yosh_init): maximum number of buckets exceeded\n") ;
      outfile->Printf( "(yosh_init): maximum number of buckets exceeded\n") ;
      outfile->Printf( "   wanted %d buckets\n", nbuckets) ;
      tstop() ;
      exit(PSI_RETURN_FAILURE) ;
      }

   /* if the number of pq does not divide evenly among the buckets, then
    * the last bucket will have the remainder of the pq's.
    */
   YBuff->pq_per_bucket = bra_indices / nbuckets ;

   if (nbuckets == 1) {
      bytes_per_bucket = ((unsigned long int) (4*sizeof(int) + sizeof(double))) *
       ((unsigned long int) twoel_array_size) + (unsigned long int) (sizeof(struct iwlbuf)
       + IWL_INTS_PER_BUF * (4*sizeof(Label) + sizeof(Value)));
    if (bytes_per_bucket > (unsigned long int) (maxcor/nbuckets))
      bytes_per_bucket = (unsigned long int) (maxcor / nbuckets);
   }
   else
      bytes_per_bucket = (unsigned long int) (maxcor / nbuckets);
   outfile->Printf("There are %i bytes per buckets\n", bytes_per_bucket);

   free_bytes_per_bucket = bytes_per_bucket -
     (unsigned long int) (sizeof(struct iwlbuf) + IWL_INTS_PER_BUF * (4*sizeof(Label) + sizeof(Value)));
   YBuff->bucketsize = free_bytes_per_bucket / (4 * sizeof(int) +
      sizeof(double));
   YBuff->buckets = (struct bucket *) malloc(nbuckets * sizeof(struct bucket));
   YBuff->bucket_for_pq = (int *) malloc (bra_indices * sizeof(int));
   for (i=0,pq=0; i<nbuckets; i++) {
      if (i != (nbuckets - 1)) {
         for (j=0; j<YBuff->pq_per_bucket; j++)
            YBuff->bucket_for_pq[pq++] = i;
            }
      else for (pq=pq; pq<bra_indices; pq++) YBuff->bucket_for_pq[pq] = i;
      }

   for (i=0,j=0; i<(nbuckets-1); i++) {
      YBuff->buckets[i].in_bucket = 0;
      YBuff->buckets[i].lo = j;
      YBuff->buckets[i].hi = YBuff->buckets[i].lo + YBuff->pq_per_bucket - 1;
      j += YBuff->pq_per_bucket;
      }
   YBuff->buckets[i].in_bucket = 0;
   YBuff->buckets[i].lo = j;
   YBuff->buckets[i].hi = bra_indices - 1;
}

/*
** YOSH_INIT_PK(): This function initializes a Yoshimine sort object for
**                 sorting PK integrals.
**    The data is contained in a structure YBuff.
**
** Parameters:
**    YBuff          =  pointer to struct which will hold data for the object
**    bra_indices    =  the number of bra_index pairs for the two-electron
**                      integrals being sorted
**    maxcor         =  the core memory, in bytes, available for all buckets
**                       simultaneously during presorting.
**    maxcord        =  the core memory, in doubles, available for reading in
**                      a single presorted bucket file during sorting.
**    max_buckets    =  the max number of buckets to use (may be limited due
**                      to the fact that each bucket requires a consecutively
**                      numbered binary temp file).
**    first_tmp_file =  the number of the first tmp file used in the
**                      Yoshimine sort (e.g. 80 for file80).
**    cutoff         =  minimum value to be kept for any value during the sort
**    outfile        =  the text output file
**
** Returns: none
**
** Note:  bra_indices indicates the number of pq pairs. In PK integrals, we know
** that we have pq >= rs, p >= q and r >= s. YBuff contains buckets with
** varying numbers of pq indices to take into account this triangular structure.
*/
void yosh_init_pk(struct yoshimine *YBuff, unsigned bra_indices,
               long maxcor, long maxcord, const int max_buckets,
               unsigned int first_tmp_file,
               double cutoff, std::string OutFileRMR)
{
   unsigned long long int twoel_array_size;              /*--- Although on 32-bit systems one can only allocate 2GB arrays
                                                           in 1 process space, one can store much bigger integrals files on disk ---*/
   unsigned int nbuckets;
   int i, j, pq;
   unsigned long int bytes_per_bucket, free_bytes_per_bucket;

   YBuff->first_tmp_file = first_tmp_file;
   twoel_array_size = bra_indices * (bra_indices + 1) / 2;
   YBuff->core_loads = (twoel_array_size - 1) / maxcord + 1 ;
   nbuckets = YBuff->core_loads ;
   outfile->Printf("nbuckets is %i\n", nbuckets );
   YBuff->nbuckets = nbuckets;
   YBuff->cutoff = cutoff;
   YBuff->bra_indices = bra_indices;
   // -1 because this should not be used here.
   YBuff->ket_indices = -1;
   if (nbuckets > max_buckets) {
      outfile->Printf( "(yosh_init): maximum number of buckets exceeded\n") ;
      outfile->Printf( "(yosh_init): maximum number of buckets exceeded\n") ;
      outfile->Printf( "   wanted %d buckets\n", nbuckets) ;
      tstop() ;
      exit(PSI_RETURN_FAILURE) ;
      }

   if (nbuckets == 1) {
      bytes_per_bucket = ((unsigned long int) (4*sizeof(int) + sizeof(double))) *
       ((unsigned long int) twoel_array_size) + (unsigned long int) (sizeof(struct iwlbuf)
       + IWL_INTS_PER_BUF * (4*sizeof(Label) + sizeof(Value)));
    if (bytes_per_bucket > (unsigned long int) (maxcor/nbuckets))
      bytes_per_bucket = (unsigned long int) (maxcor / nbuckets);
   }
   else
      bytes_per_bucket = (unsigned long int) (maxcor / nbuckets);
   outfile->Printf("There are %i bytes per buckets\n", bytes_per_bucket);

   free_bytes_per_bucket = bytes_per_bucket -
     (unsigned long int) (sizeof(struct iwlbuf) + IWL_INTS_PER_BUF * (4*sizeof(Label) + sizeof(Value)));
   YBuff->bucketsize = free_bytes_per_bucket / (4 * sizeof(int) +
      sizeof(double));
   YBuff->buckets = (struct bucket *) malloc(nbuckets * sizeof(struct bucket));
   YBuff->bucket_for_pq = (int *) malloc (bra_indices * sizeof(int));

   /* Now we determine which buckets have which pq indices. Since we have
    * integrals in triangular form, low pq indices have very few rs indices,
    * thus buckets will have varying numbers of pq indices to handle.
    */
   // pq_per_bucket should not be used here.
   YBuff->pq_per_bucket = -1;

   // I *think* all of this should work the same for J and K
   // sorting
   unsigned long long int pq_incore = 0;
   i = 0;
   YBuff->buckets[i].lo = 0;
   YBuff->buckets[i].in_bucket = 0;
   for(unsigned pq = 0; pq < bra_indices; ++pq) {
       // Increment counters
       if(pq_incore + pq + 1 > maxcord) {
           //The batch is full. Save info.
           YBuff->buckets[i++].hi = pq - 1;
           YBuff->buckets[i].in_bucket = 0;
           YBuff->buckets[i].lo = pq;
           pq_incore = 0;
       }
       YBuff->bucket_for_pq[pq] = i;
       pq_incore += pq + 1;
   }
   YBuff->buckets[i].hi = bra_indices - 1;

}

/* YOSH_DONE(): Free allocated memory and reset all options for a
** given Yoshimine sorting structure.
*/
void yosh_done(struct yoshimine *YBuff)
{
  YBuff->core_loads = 0;
  YBuff->nbuckets = 0;
  free(YBuff->bucket_for_pq);
  YBuff->bucketsize = 0;
  free(YBuff->buckets);
  YBuff->first_tmp_file = 0;
  YBuff->pq_per_bucket = 0;
  YBuff->bra_indices = 0;
  YBuff->ket_indices = 0;
  YBuff->cutoff = 0;
}

void yosh_init_buckets(struct yoshimine *YBuff)
{
  int i;

   for (i=0; i<(YBuff->nbuckets); i++) {
      YBuff->buckets[i].p = init_int_array(YBuff->bucketsize);
      YBuff->buckets[i].q = init_int_array(YBuff->bucketsize);
      YBuff->buckets[i].r = init_int_array(YBuff->bucketsize);
      YBuff->buckets[i].s = init_int_array(YBuff->bucketsize);
      YBuff->buckets[i].val = init_array(YBuff->bucketsize);
      iwl_buf_init(&(YBuff->buckets[i].IWLBuf),YBuff->first_tmp_file+i,
                   YBuff->cutoff,0, 0);
      }
}

/*
** YOSH_PRINT(): Print out the Yoshimine structure
**
*/
void yosh_print(struct yoshimine *YBuff, std::string OutFileRMR)
{
   outfile->Printf( " Yoshimine structure:\n");
   outfile->Printf( "\tbra_indices  = %10d\n", YBuff->bra_indices);
   outfile->Printf( "\tket_indices  = %10d\n", YBuff->ket_indices);
   outfile->Printf( "\tbin size   = %10lu\n", YBuff->bucketsize) ;
   outfile->Printf( "\tbins       = %10d\n", YBuff->nbuckets) ;
   outfile->Printf( "\tcore loads = %10d\n", YBuff->core_loads) ;
   outfile->Printf( "\tpq/bin     = %10d\n", YBuff->pq_per_bucket) ;
   outfile->Printf( "\tcutoff     = %10.2E\n", YBuff->cutoff) ;
   outfile->Printf( " Structure for each bucket:\n");
   for(int i = 0; i < YBuff->nbuckets; ++i ) {
       outfile->Printf(" Bucket number %i\n", i);
       outfile->Printf("\tLow pq = %10d\n", YBuff->buckets[i].lo);
       outfile->Printf("\tHigh pq = %10d\n", YBuff->buckets[i].hi);
   }
}


/*
** YOSH_CLOSE_BUCKETS(): Close the temporary binary files used for the
** Yoshimine sort (not the final output file).
**
** Arguments:
**    YBuff   =  pointer to yoshimine object
**    erase   =  1 to erase tmp files, else 0
**
** Returns: none
**
** This function was formerly called yosh_close_tmp_files().
*/
void yosh_close_buckets(struct yoshimine *YBuff, int erase)
{
   int i;

   for (i=0; i<YBuff->nbuckets; i++) { /* close but keep */
      iwl_buf_close(&(YBuff->buckets[i].IWLBuf), !erase);
      free(YBuff->buckets[i].p);
      free(YBuff->buckets[i].q);
      free(YBuff->buckets[i].r);
      free(YBuff->buckets[i].s);
      free(YBuff->buckets[i].val);
      }
}



/*
** YOSH_RDTWO() : Read two-electron integrals from
**    file33 (in IWL format) and prepare them for Yoshimine sorting.
**
** Adapted from Ed Seidl's CSCF rdtwo.c routine
** David Sherrill
** Center for Computational Quantum Chemistry, UGA
**
** Created 1993
** Modified February 1995 to use new Yoshimine data structure
** Modified March 1995 to use nsoff array instead of call to abs_orb()
**
** Arguments:
**   YBuff        = Yoshimine object pointer
**   itapERI      = unit number for two el. file (33)
**   num_so       = array of number of symm orbs in each irrep (for reindex)
**   nirreps      = number of irreps
**   ioff         = standard lexical index array
**   elbert       = 1 for Elbert ordering, 0 for canonical ordering
**   P            = frozen core density matrix (lower triangle)
**   Hc           = frozen core operator (lower triangle)
**   matrix       = 1 for all rs for given pq, 0 otherwise
**                  (for matrix multiplication algorithm)
**   del_tei_file = 1 to delete the tei file (33), 0 otherwise
**   printflag    = 1 for printing (for debugging only!) else 0
**   outfile      = file to print integrals to (if printflag is set)
*/
void yosh_rdtwo(struct yoshimine *YBuff, int itapERI, int del_tei_file, int *num_so,
      int nirreps, int *ioff, int elbert, int fzcflag, double *P, double *Hc,
      int matrix, int printflag, std::string OutFileRMR)
{
  int ilsti, nbuf;
  int i, ij, kl, ijkl;
  int ior, ism, jor, jsm;
  int kor, ksm, lor, lsm;
  int iabs, jabs, kabs, labs ;
  int d2i ;
  double value;
  int *tmp;
  struct bucket *bptr ;
  long int tmpi;
  int whichbucket, lastflag = 0, firstfile;
  int *nsoff;
  int a,b,c,d,ab,cd,ad,bc,dum,found=0;
  int al[8], bl[8], cl[8], dl[8];
  int fi;
  struct iwlbuf ERIIN;

  if (printflag) {
    outfile->Printf( "Yoshimine rdtwo routine entered\n");
    outfile->Printf( "Two-electron integrals from file%d:\n",itapERI);
  }

  firstfile = YBuff->first_tmp_file;

  iwl_buf_init(&ERIIN,itapERI,0.0,1,1);

  nsoff = init_int_array(nirreps);
  nsoff[0] = 0;
  for (i=1; i<nirreps; i++) {
    nsoff[i] = nsoff[i-1] + num_so[i-1];
  }

  do {
    /* read a buffer full */
    ilsti = ERIIN.lastbuf;
    nbuf = ERIIN.inbuf;

    fi = 0;
    for (i=0; i < nbuf ; i++,tmp += 2) { /* do funky stuff to unpack ints */
      iabs = abs(ERIIN.labels[fi]);
      jabs = ERIIN.labels[fi+1];
      kabs = ERIIN.labels[fi+2];
      labs = ERIIN.labels[fi+3];
      value = ERIIN.values[i];
      fi += 4;

      /* calculate ijkl lexical index */
      ij = ioff[iabs] + jabs;
      kl = ioff[kabs] + labs;
      ijkl = ioff[ij] + kl;

      /* newly added March 1995: construct the frozen core operator */
      if (fzcflag) {
        a = al[0] = iabs;
        b = bl[0] = jabs;
        c = cl[0] = kabs;
        d = dl[0] = labs;
        ab = ioff[MAX0(a,b)] + MIN0(a,b);
        cd = ioff[MAX0(c,d)] + MIN0(c,d);
        bc = ioff[MAX0(b,c)] + MIN0(b,c);
        ad = ioff[MAX0(a,d)] + MIN0(a,d);
        Hc[cd] += 2.0 * P[ab] * value;
        if (b >= c) Hc[bc] -= P[ad] * value;

        a = al[1] = jabs;
        b = bl[1] = iabs;
        c = cl[1] = kabs;
        d = dl[1] = labs;
        if (!(a == al[0] && b == bl[0] && c == cl[0] && d == dl[0])) {
          ab = ioff[MAX0(a,b)] + MIN0(a,b);
          cd = ioff[MAX0(c,d)] + MIN0(c,d);
          bc = ioff[MAX0(b,c)] + MIN0(b,c);
          ad = ioff[MAX0(a,d)] + MIN0(a,d);
          if (c >= d) Hc[cd] += 2.0 * P[ab] * value;
          if (b >= c) Hc[bc] -= P[ad] * value;
        }

        a = al[2] = iabs;
        b = bl[2] = jabs;
        c = cl[2] = labs;
        d = dl[2] = kabs;
        for (dum=0, found=0; dum < 2 && !found; dum++) {
          if (a==al[dum] && b==bl[dum] && c==cl[dum] && d==dl[dum]) found=1;
        }
        if (!found) {
          ab = ioff[MAX0(a,b)] + MIN0(a,b);
          cd = ioff[MAX0(c,d)] + MIN0(c,d);
          bc = ioff[MAX0(b,c)] + MIN0(b,c);
          ad = ioff[MAX0(a,d)] + MIN0(a,d);
          if (c >= d) Hc[cd] += 2.0 * P[ab] * value;
          if (b >= c) Hc[bc] -= P[ad] * value;
        }

        a = al[3] = jabs;
        b = bl[3] = iabs;
        c = cl[3] = labs;
        d = dl[3] = kabs;
        for (dum=0, found=0; dum < 3 && !found; dum++) {
          if(a==al[dum] && b==bl[dum] && c==cl[dum] && d==dl[dum]) found=1;
        }
        if (!found) {
          ab = ioff[MAX0(a,b)] + MIN0(a,b);
          cd = ioff[MAX0(c,d)] + MIN0(c,d);
          bc = ioff[MAX0(b,c)] + MIN0(b,c);
          ad = ioff[MAX0(a,d)] + MIN0(a,d);
          if (c >= d) Hc[cd] += 2.0 * P[ab] * value;
          if (b >= c) Hc[bc] -= P[ad] * value;
        }

        a = al[4] = kabs;
        b = bl[4] = labs;
        c = cl[4] = iabs;
        d = dl[4] = jabs;
        for (dum=0, found=0; dum < 4 && !found; dum++) {
          if(a==al[dum] && b==bl[dum] && c==cl[dum] && d==dl[dum]) found=1;
        }
        if (!found) {
          ab = ioff[MAX0(a,b)] + MIN0(a,b);
          cd = ioff[MAX0(c,d)] + MIN0(c,d);
          bc = ioff[MAX0(b,c)] + MIN0(b,c);
          ad = ioff[MAX0(a,d)] + MIN0(a,d);
          if (c >= d) Hc[cd] += 2.0 * P[ab] * value;
          if (b >= c) Hc[bc] -= P[ad] * value;
        }

        a = al[5] = kabs;
        b = bl[5] = labs;
        c = cl[5] = jabs;
        d = dl[5] = iabs;
        for (dum=0, found=0; dum < 5 && !found; dum++) {
          if (a==al[dum] && b==bl[dum] && c==cl[dum] && d==dl[dum]) found=1;
        }
        if (!found) {
          ab = ioff[MAX0(a,b)] + MIN0(a,b);
          cd = ioff[MAX0(c,d)] + MIN0(c,d);
          bc = ioff[MAX0(b,c)] + MIN0(b,c);
          ad = ioff[MAX0(a,d)] + MIN0(a,d);
          if (c >= d) Hc[cd] += 2.0 * P[ab] * value;
          if (b >= c) Hc[bc] -= P[ad] * value;
        }

        a = al[6] = labs;
        b = bl[6] = kabs;
        c = cl[6] = iabs;
        d = dl[6] = jabs;
        for (dum=0, found=0; dum < 6 && !found; dum++) {
          if (a==al[dum] && b==bl[dum] && c==cl[dum] && d==dl[dum]) found=1;
        }
        if (!found) {
          ab = ioff[MAX0(a,b)] + MIN0(a,b);
          cd = ioff[MAX0(c,d)] + MIN0(c,d);
          bc = ioff[MAX0(b,c)] + MIN0(b,c);
          ad = ioff[MAX0(a,d)] + MIN0(a,d);
          if (c >= d) Hc[cd] += 2.0 * P[ab] * value;
          if (b >= c) Hc[bc] -= P[ad] * value;
        }

        a = al[7] = labs;
        b = bl[7] = kabs;
        c = cl[7] = jabs;
        d = dl[7] = iabs;
        for (dum=0, found=0; dum < 7 && !found; dum++) {
          if(a==al[dum] && b==bl[dum] && c==cl[dum] && d==dl[dum]) found=1;
        }
        if (!found) {
          ab = ioff[MAX0(a,b)] + MIN0(a,b);
          cd = ioff[MAX0(c,d)] + MIN0(c,d);
          bc = ioff[MAX0(b,c)] + MIN0(b,c);
          ad = ioff[MAX0(a,d)] + MIN0(a,d);
          if (c >= d) Hc[cd] += 2.0 * P[ab] * value;
          if (b >= c) Hc[bc] -= P[ad] * value;
        }
      } /* end construction of frozen core operator */

      /* figure out what bucket to put it in, and do so
       *
       * Elbert wants us to sort by the lower index (kl)
       * i.e. for us, ij > kl (guaranteed in 33), but for them kl > ij
       *
       */

      if (elbert) whichbucket = YBuff->bucket_for_pq[kl] ;
      else whichbucket = YBuff->bucket_for_pq[ij] ;

      bptr = YBuff->buckets+whichbucket ;
      tmpi = (bptr->in_bucket)++ ;

      /* switch things around here for Elbert (k->p, l->q, i->r, j->s) */
      if (elbert) {
        bptr->p[tmpi] = kabs;
        bptr->q[tmpi] = labs;
        bptr->r[tmpi] = iabs;
        bptr->s[tmpi] = jabs;
      }
      else {
        bptr->p[tmpi] = iabs;
        bptr->q[tmpi] = jabs;
        bptr->r[tmpi] = kabs;
        bptr->s[tmpi] = labs;
      }

      bptr->val[tmpi] = value;

      if (printflag)
        outfile->Printf( "%4d %4d %4d %4d  %4d   %10.6lf\n",
                iabs, jabs, kabs, labs, ijkl, value) ;
      if ((tmpi+1) == YBuff->bucketsize) { /* need to flush bucket to disk */
        flush_bucket(bptr, 0);
        bptr->in_bucket = 0;
      }

      if(matrix && ij != kl) {
        whichbucket = YBuff->bucket_for_pq[kl] ;
        bptr = YBuff->buckets+whichbucket ;
        tmpi = (bptr->in_bucket)++;
        bptr->p[tmpi] = kabs;
        bptr->q[tmpi] = labs;
        bptr->r[tmpi] = iabs;
        bptr->s[tmpi] = jabs;
        bptr->val[tmpi] = value;
        if ((tmpi+1) == YBuff->bucketsize) {
          flush_bucket(bptr, 0);
          bptr->in_bucket = 0;
        }
      }

    }
    if (!ilsti)
      iwl_buf_fetch(&ERIIN);
  } while(!ilsti);

  /* flush partially filled buckets */
  /* Ok, after "matrix" was added above, we ran into the possibility of
   * flushing TWO buffers with the lastflag set.  This would be bad,
   * because the second buffer would never be read.  Therefore, I have
   * always passed a lastflag of 0 to flush_bucket() in the code above,
   * and now I flush all buckets here with lastflag set to 1.  There
   * is a small possibility that I will write a buffer of all zeroes.
   * This should not actually cause a problem, the way the iwl buf reads
   * currently work.  Make sure to be careful if rewriting iwl routines!
   */
  for (i=0; i<YBuff->nbuckets; i++) {
    flush_bucket((YBuff->buckets)+i, 1);
  }

  free(nsoff);

  iwl_buf_close(&ERIIN, !del_tei_file);
}

/*
** YOSH_RDTWO_PK() : Read two-electron integrals from
**    file33 (in IWL format) and prepare them for Yoshimine sorting
**    for PK integrals.
**
**    We need to sort Coulomb using ij indices, but exchange using
**    ik and il indices, thus we need two Yoshimine objects, and we
**    write to two sets of temporary files, while only reading the
**    integral file once.
**
** Arguments:
**   YBuffJ       = Yoshimine object pointer for Coulomb integrals
**   YBuffK       = Yoshimine object pointer for exchange integrals
**   itapERI      = unit number for two el. file (33)
**   nirreps      = number of irreps
**   so2rel       = array mapping absolute basis function index to relative
**                  basis function index within an irrep, so2rel[abs] = rel
**   so2sym       = array mapping absolute basis function index to irrep
**                  number, so2sym[abs] = sym
**   pksymoff     = array containing the offset in each irrep to convert a
**                  pq index computed with relative indices to an absolute
**                  pq index, pqrel = ioff[prel] + qrel, pqabs = pqrel + pksymoff[psym]
**   ioff         = standard lexical index array
**   del_tei_file = 1 to delete the tei file (33), 0 otherwise
**   printflag    = 1 for printing (for debugging only!) else 0
*/
void yosh_rdtwo_pk(struct yoshimine *YBuffJ, struct yoshimine *YBuffK, int itapERI,
      int del_tei_file, int nirreps, int* so2rel, int* so2sym, int* pksymoff, int *ioff, int printflag)
{
  int ilsti, nbuf;
  int i, ij, kl, ijkl;
  int ik, il;
  int iabs, jabs, kabs, labs ;
  int irel, jrel, krel, lrel;
  int isym, jsym, ksym, lsym;
  double value;
  int *tmp;
  struct bucket *bptr_J, *bptr_K1, *bptr_K2 ;
  long int tmpi_J, tmpi_K1, tmpi_K2;
  int whichbucket_J, whichbucket_K1, whichbucket_K2, firstfile_J, firstfile_K;
  int fi;
  struct iwlbuf ERIIN;
  int* num_int_J;
  int* num_int_K;

  num_int_J = init_int_array(YBuffJ->nbuckets);
  num_int_K = init_int_array(YBuffK->nbuckets);

  for(int h  = 0; h < YBuffJ->nbuckets; ++h) {
      num_int_J[h] = 0;
  }
  for(int h  = 0; h < YBuffK->nbuckets; ++h) {
      num_int_K[h] = 0;
  }

  if (printflag) {
    outfile->Printf( "Yoshimine rdtwo routine entered\n");
    outfile->Printf( "Two-electron integrals from file%d:\n",itapERI);
  }

  firstfile_J = YBuffJ->first_tmp_file;
  firstfile_K = YBuffK->first_tmp_file;

  iwl_buf_init(&ERIIN,itapERI,0.0,1,1);

  do {
    /* read a buffer full */
    ilsti = ERIIN.lastbuf;
    nbuf = ERIIN.inbuf;

    fi = 0;
    for (i=0; i < nbuf ; i++,tmp += 2) { /* do funky stuff to unpack ints */
      iabs = abs(ERIIN.labels[fi]);
      jabs = ERIIN.labels[fi+1];
      kabs = ERIIN.labels[fi+2];
      labs = ERIIN.labels[fi+3];

      irel = so2rel[iabs];
      jrel = so2rel[jabs];
      krel = so2rel[kabs];
      lrel = so2rel[labs];

      isym = so2sym[iabs];
      jsym = so2sym[jabs];
      ksym = so2sym[kabs];
      lsym = so2sym[labs];

      value = ERIIN.values[i];
      fi += 4;

      // Need to add the logic for symmetry of the indices here

      // K for first sort IKJL
      if ((isym == ksym) && (jsym == lsym)) {
          /* Calculate ik for exchange buckets */
          ik = ioff[irel] + krel;
          ik += pksymoff[isym];

          /* figure out what bucket to put it in, and do so
           */

          whichbucket_K1 = YBuffK->bucket_for_pq[ik] ;

          bptr_K1= YBuffK->buckets + whichbucket_K1 ;
          tmpi_K1 = (bptr_K1->in_bucket)++ ;

          // Fill the first exchange bucket.
          bptr_K1->p[tmpi_K1] = iabs;
          bptr_K1->q[tmpi_K1] = jabs;
          bptr_K1->r[tmpi_K1] = kabs;
          bptr_K1->s[tmpi_K1] = labs;

          bptr_K1->val[tmpi_K1] = value;
          num_int_K[whichbucket_K1]++;
      }
      if ((tmpi_K1 + 1) == YBuffK->bucketsize) { /* need to flush bucket to disk */
          flush_bucket(bptr_K1, 0);
          bptr_K1->in_bucket = 0;
      }


      // We need the symmetry offset to get pqrs properly.
      if((isym == jsym) && (ksym == lsym)) {
          /* calculate ijkl lexical index */
          ij = ioff[irel] + jrel;
          ij += pksymoff[isym];
          kl = ioff[krel] + lrel;
          kl += pksymoff[ksym];
          // ijkl only here for debug printing
          ijkl = ioff[ij] + kl;
          whichbucket_J = YBuffJ->bucket_for_pq[ij] ;
          bptr_J = YBuffJ->buckets + whichbucket_J ;
          tmpi_J = (bptr_J->in_bucket)++ ;

          // Fill the Coulomb bucket
          bptr_J->p[tmpi_J] = iabs;
          bptr_J->q[tmpi_J] = jabs;
          bptr_J->r[tmpi_J] = kabs;
          bptr_J->s[tmpi_J] = labs;

          bptr_J->val[tmpi_J] = value;
          num_int_J[whichbucket_J]++;

          if (printflag)
            outfile->Printf( "%4d %4d %4d %4d  %4d   %10.6lf\n",
                    iabs, jabs, kabs, labs, ijkl, value) ;

          // Now we do the second sort for K (ILJK), which apparently should
          // be in there.

          if((irel != jrel) && (krel != lrel)) {
              if((isym == lsym) && (jsym == ksym)) {
                  /* Calculate il for exchange buckets */
                  il = ioff[irel] + lrel;
                  il += pksymoff[isym];

                  // outfile->Printf("ik is %i and il is %i\n", ik, il);
                  // outfile->Printf("for integral <%i %i |%i %i>", iabs, jabs, kabs, labs);
                  whichbucket_K2 = YBuffK->bucket_for_pq[il] ;
                  if (whichbucket_K1 != whichbucket_K2) {
                      bptr_K2 = YBuffK->buckets + whichbucket_K2 ;
                      tmpi_K2 = (bptr_K2->in_bucket)++ ;

                      //outfile->Printf("We do not skip\n");
                      bptr_K2->p[tmpi_K2] = iabs;
                      bptr_K2->q[tmpi_K2] = jabs;
                      bptr_K2->r[tmpi_K2] = kabs;
                      bptr_K2->s[tmpi_K2] = labs;

                      bptr_K2->val[tmpi_K2] = value;
                      num_int_K[whichbucket_K2]++;

                      if ((tmpi_K2 + 1) == YBuffK->bucketsize) { /* need to flush bucket to disk */
                          flush_bucket(bptr_K2, 0);
                          bptr_K2->in_bucket = 0;
                      }
                  }

              }
          }
      } else if((isym == lsym) && (jsym == ksym)) {
          outfile->Printf("This integral should be in K 2nd sort\n");
      }
      if ((tmpi_J + 1) == YBuffJ->bucketsize) { /* need to flush bucket to disk */
        flush_bucket(bptr_J, 0);
        bptr_J->in_bucket = 0;
      }
    }
    if (!ilsti)
      iwl_buf_fetch(&ERIIN);
  } while(!ilsti);

  /* flush partially filled buckets */
  // TODO: For PK, we never use matrix so we could fix that.
  /* Ok, after "matrix" was added above, we ran into the possibility of
   * flushing TWO buffers with the lastflag set.  This would be bad,
   * because the second buffer would never be read.  Therefore, I have
   * always passed a lastflag of 0 to flush_bucket() in the code above,
   * and now I flush all buckets here with lastflag set to 1.  There
   * is a small possibility that I will write a buffer of all zeroes.
   * This should not actually cause a problem, the way the iwl buf reads
   * currently work.  Make sure to be careful if rewriting iwl routines!
   */
  for (i=0; i<YBuffJ->nbuckets; i++) {
    flush_bucket((YBuffJ->buckets)+i, 1);
  }

  for (i=0; i<YBuffK->nbuckets; i++) {
    flush_bucket((YBuffK->buckets)+i, 1);
  }

  for (int h = 0; h < YBuffJ->nbuckets; ++h) {
      outfile->Printf("We wrote %i integrals in J bucket %i\n", num_int_J[h], h);
  }
  for (int h = 0; h < YBuffK->nbuckets; ++h) {
      outfile->Printf("We wrote %i integrals in K bucket %i\n", num_int_K[h], h);
  }

  free(num_int_J);
  free(num_int_K);

  iwl_buf_close(&ERIIN, !del_tei_file);
}

/*
** YOSH_RDTWO_UHF() : Read two-electron integrals from file33 (in IWL
** format) and prepare them for Yoshimine sorting.
**
** Arguments:
**   YBuff        = Yoshimine object pointer
**   itapERI      = unit number for two el. file (33)
**   num_so       = array of number of symm orbs in each irrep (for reindex)
**   nirreps      = number of irreps
**   ioff         = standard lexical index array
**   elbert       = 1 for Elbert ordering, 0 for canonical ordering
**   Pa           = alpha frozen core density matrix (lower triangle)
**   Pb           = beta frozen core density matrix (lower triangle)
**   Hca          = alpha frozen core operator (lower triangle)
**   Hcb          = beta frozen core operator (lower triangle)
**   matrix       = 1 for all rs for given pq, 0 otherwise
**                  (for matrix multiplication algorithm)
**   del_tei_file = 1 to delete the tei file (33), 0 otherwise
**   printflag    = 1 for printing (for debugging only!) else 0
**   outfile      = file to print integrals to (if printflag is set)
*/
void yosh_rdtwo_uhf(struct yoshimine *YBuff, int itapERI, int del_tei_file, int *num_so,
                    int nirreps, int *ioff, int elbert, int fzcflag, double *Pa, double *Pb,
                    double *Hca, double *Hcb, int matrix, int printflag, std::string OutFileRMR)
{
  int ilsti, nbuf;
  int i, ij, kl, ijkl;
  int ior, ism, jor, jsm;
  int kor, ksm, lor, lsm;
  int iabs, jabs, kabs, labs ;
  int d2i ;
  double value;
  int *tmp;
  struct bucket *bptr ;
  long int tmpi;
  int whichbucket, lastflag = 0, firstfile;
  int *nsoff;
  int a,b,c,d,ab,cd,ad,bc,dum,found=0;
  int al[8], bl[8], cl[8], dl[8];
  int fi;
  struct iwlbuf ERIIN;

  if (printflag) {
    outfile->Printf( "Yoshimine rdtwo routine entered\n");
    outfile->Printf( "Two-electron integrals from file%d:\n",itapERI);
  }

  firstfile = YBuff->first_tmp_file;

  iwl_buf_init(&ERIIN,itapERI,0.0,1,1);

  nsoff = init_int_array(nirreps);
  nsoff[0] = 0;
  for (i=1; i<nirreps; i++) {
    nsoff[i] = nsoff[i-1] + num_so[i-1];
  }

  do {
    /* read a buffer full */
    ilsti = ERIIN.lastbuf;
    nbuf = ERIIN.inbuf;

    fi = 0;
    for (i=0; i < nbuf ; i++,tmp += 2) { /* do funky stuff to unpack ints */
      iabs = abs(ERIIN.labels[fi]);
      jabs = ERIIN.labels[fi+1];
      kabs = ERIIN.labels[fi+2];
      labs = ERIIN.labels[fi+3];
      value = ERIIN.values[i];
      fi += 4;

      /* calculate ijkl lexical index */
      ij = ioff[iabs] + jabs;
      kl = ioff[kabs] + labs;
      ijkl = ioff[ij] + kl;

      /* construct the UHF frozen core operator */
      if (fzcflag) {
        a = al[0] = iabs;
        b = bl[0] = jabs;
        c = cl[0] = kabs;
        d = dl[0] = labs;
        ab = ioff[MAX0(a,b)] + MIN0(a,b);
        cd = ioff[MAX0(c,d)] + MIN0(c,d);
        bc = ioff[MAX0(b,c)] + MIN0(b,c);
        ad = ioff[MAX0(a,d)] + MIN0(a,d);
        Hca[cd] += (Pa[ab] + Pb[ab]) * value;
        Hcb[cd] += (Pa[ab] + Pb[ab]) * value;
        if (b >= c) {
          Hca[bc] -= Pa[ad] * value;
          Hcb[bc] -= Pb[ad] * value;
        }

        a = al[1] = jabs;
        b = bl[1] = iabs;
        c = cl[1] = kabs;
        d = dl[1] = labs;
        if (!(a == al[0] && b == bl[0] && c == cl[0] && d == dl[0])) {
          ab = ioff[MAX0(a,b)] + MIN0(a,b);
          cd = ioff[MAX0(c,d)] + MIN0(c,d);
          bc = ioff[MAX0(b,c)] + MIN0(b,c);
          ad = ioff[MAX0(a,d)] + MIN0(a,d);
          if (c >= d) {
            Hca[cd] += (Pa[ab] + Pb[ab]) * value;
            Hcb[cd] += (Pa[ab] + Pb[ab]) * value;
          }
          if (b >= c) {
            Hca[bc] -= Pa[ad] * value;
            Hcb[bc] -= Pb[ad] * value;
          }
        }

        a = al[2] = iabs;
        b = bl[2] = jabs;
        c = cl[2] = labs;
        d = dl[2] = kabs;
        for (dum=0, found=0; dum < 2 && !found; dum++) {
          if (a==al[dum] && b==bl[dum] && c==cl[dum] && d==dl[dum]) found=1;
        }
        if (!found) {
          ab = ioff[MAX0(a,b)] + MIN0(a,b);
          cd = ioff[MAX0(c,d)] + MIN0(c,d);
          bc = ioff[MAX0(b,c)] + MIN0(b,c);
          ad = ioff[MAX0(a,d)] + MIN0(a,d);
          if (c >= d) {
            Hca[cd] += (Pa[ab] + Pb[ab]) * value;
            Hcb[cd] += (Pa[ab] + Pb[ab]) * value;
          }
          if (b >= c) {
            Hca[bc] -= Pa[ad] * value;
            Hcb[bc] -= Pb[ad] * value;
          }
        }

        a = al[3] = jabs;
        b = bl[3] = iabs;
        c = cl[3] = labs;
        d = dl[3] = kabs;
        for (dum=0, found=0; dum < 3 && !found; dum++) {
          if(a==al[dum] && b==bl[dum] && c==cl[dum] && d==dl[dum]) found=1;
        }
        if (!found) {
          ab = ioff[MAX0(a,b)] + MIN0(a,b);
          cd = ioff[MAX0(c,d)] + MIN0(c,d);
          bc = ioff[MAX0(b,c)] + MIN0(b,c);
          ad = ioff[MAX0(a,d)] + MIN0(a,d);
          if (c >= d) {
            Hca[cd] += (Pa[ab] + Pb[ab]) * value;
            Hcb[cd] += (Pa[ab] + Pb[ab]) * value;
          }
          if (b >= c) {
            Hca[bc] -= Pa[ad] * value;
            Hcb[bc] -= Pb[ad] * value;
          }
        }

        a = al[4] = kabs;
        b = bl[4] = labs;
        c = cl[4] = iabs;
        d = dl[4] = jabs;
        for (dum=0, found=0; dum < 4 && !found; dum++) {
          if(a==al[dum] && b==bl[dum] && c==cl[dum] && d==dl[dum]) found=1;
        }
        if (!found) {
          ab = ioff[MAX0(a,b)] + MIN0(a,b);
          cd = ioff[MAX0(c,d)] + MIN0(c,d);
          bc = ioff[MAX0(b,c)] + MIN0(b,c);
          ad = ioff[MAX0(a,d)] + MIN0(a,d);
          if (c >= d) {
            Hca[cd] += (Pa[ab] + Pb[ab]) * value;
            Hcb[cd] += (Pa[ab] + Pb[ab]) * value;
          }
          if (b >= c) {
            Hca[bc] -= Pa[ad] * value;
            Hcb[bc] -= Pb[ad] * value;
          }
        }

        a = al[5] = kabs;
        b = bl[5] = labs;
        c = cl[5] = jabs;
        d = dl[5] = iabs;
        for (dum=0, found=0; dum < 5 && !found; dum++) {
          if (a==al[dum] && b==bl[dum] && c==cl[dum] && d==dl[dum]) found=1;
        }
        if (!found) {
          ab = ioff[MAX0(a,b)] + MIN0(a,b);
          cd = ioff[MAX0(c,d)] + MIN0(c,d);
          bc = ioff[MAX0(b,c)] + MIN0(b,c);
          ad = ioff[MAX0(a,d)] + MIN0(a,d);
          if (c >= d) {
            Hca[cd] += (Pa[ab] + Pb[ab]) * value;
            Hcb[cd] += (Pa[ab] + Pb[ab]) * value;
          }
          if (b >= c) {
            Hca[bc] -= Pa[ad] * value;
            Hcb[bc] -= Pb[ad] * value;
          }
        }

        a = al[6] = labs;
        b = bl[6] = kabs;
        c = cl[6] = iabs;
        d = dl[6] = jabs;
        for (dum=0, found=0; dum < 6 && !found; dum++) {
          if (a==al[dum] && b==bl[dum] && c==cl[dum] && d==dl[dum]) found=1;
        }
        if (!found) {
          ab = ioff[MAX0(a,b)] + MIN0(a,b);
          cd = ioff[MAX0(c,d)] + MIN0(c,d);
          bc = ioff[MAX0(b,c)] + MIN0(b,c);
          ad = ioff[MAX0(a,d)] + MIN0(a,d);
          if (c >= d) {
            Hca[cd] += (Pa[ab] + Pb[ab]) * value;
            Hcb[cd] += (Pa[ab] + Pb[ab]) * value;
          }
          if (b >= c) {
            Hca[bc] -= Pa[ad] * value;
            Hcb[bc] -= Pb[ad] * value;
          }
        }

        a = al[7] = labs;
        b = bl[7] = kabs;
        c = cl[7] = jabs;
        d = dl[7] = iabs;
        for (dum=0, found=0; dum < 7 && !found; dum++) {
          if(a==al[dum] && b==bl[dum] && c==cl[dum] && d==dl[dum]) found=1;
        }
        if (!found) {
          ab = ioff[MAX0(a,b)] + MIN0(a,b);
          cd = ioff[MAX0(c,d)] + MIN0(c,d);
          bc = ioff[MAX0(b,c)] + MIN0(b,c);
          ad = ioff[MAX0(a,d)] + MIN0(a,d);
          if (c >= d) {
            Hca[cd] += (Pa[ab] + Pb[ab]) * value;
            Hcb[cd] += (Pa[ab] + Pb[ab]) * value;
          }
          if (b >= c) {
            Hca[bc] -= Pa[ad] * value;
            Hcb[bc] -= Pb[ad] * value;
          }
        }
      } /* end construction of frozen core operator */

      /* figure out what bucket to put it in, and do so
       *
       * Elbert wants us to sort by the lower index (kl)
       * i.e. for us, ij > kl (guaranteed in 33), but for them kl > ij
       *
       */

      if (elbert) whichbucket = YBuff->bucket_for_pq[kl] ;
      else whichbucket = YBuff->bucket_for_pq[ij] ;

      bptr = YBuff->buckets+whichbucket ;
      tmpi = (bptr->in_bucket)++ ;

      /* switch things around here for Elbert (k->p, l->q, i->r, j->s) */
      if (elbert) {
        bptr->p[tmpi] = kabs;
        bptr->q[tmpi] = labs;
        bptr->r[tmpi] = iabs;
        bptr->s[tmpi] = jabs;
      }
      else {
        bptr->p[tmpi] = iabs;
        bptr->q[tmpi] = jabs;
        bptr->r[tmpi] = kabs;
        bptr->s[tmpi] = labs;
      }

      bptr->val[tmpi] = value;

      if (printflag)
        outfile->Printf( "%4d %4d %4d %4d  %4d   %10.6lf\n",
                iabs, jabs, kabs, labs, ijkl, value) ;
      if ((tmpi+1) == YBuff->bucketsize) { /* need to flush bucket to disk */
        flush_bucket(bptr, 0);
        bptr->in_bucket = 0;
      }

      if(matrix && ij != kl) {
        whichbucket = YBuff->bucket_for_pq[kl] ;
        bptr = YBuff->buckets+whichbucket ;
        tmpi = (bptr->in_bucket)++;
        bptr->p[tmpi] = kabs;
        bptr->q[tmpi] = labs;
        bptr->r[tmpi] = iabs;
        bptr->s[tmpi] = jabs;
        bptr->val[tmpi] = value;
        if ((tmpi+1) == YBuff->bucketsize) {
          flush_bucket(bptr, 0);
          bptr->in_bucket = 0;
        }
      }

    }
    if (!ilsti)
      iwl_buf_fetch(&ERIIN);
  } while(!ilsti);

  /* flush partially filled buckets */
  /* Ok, after "matrix" was added above, we ran into the possibility of
   * flushing TWO buffers with the lastflag set.  This would be bad,
   * because the second buffer would never be read.  Therefore, I have
   * always passed a lastflag of 0 to flush_bucket() in the code above,
   * and now I flush all buckets here with lastflag set to 1.  There
   * is a small possibility that I will write a buffer of all zeroes.
   * This should not actually cause a problem, the way the iwl buf reads
   * currently work.  Make sure to be careful if rewriting iwl routines!
   */
  for (i=0; i<YBuff->nbuckets; i++) {
    flush_bucket((YBuff->buckets)+i, 1);
  }

  free(nsoff);
  iwl_buf_close(&ERIIN, !del_tei_file);
}

/*
** YOSH_RDTWO_BACKTR() : Read two-electron integrals from an IWL file and
**    prepare them for Yoshimine sorting.   We have removed support for
**    Elbert loops and frozen core, since the former is not currently
**    being used and the latter should not apply to backtransforms.
**    The main issue is whether we need to symmetrize the twopdm, or
**    whether it has already been symmetrized.  We assume that it always
**    has the symmetry (pq|rs) = (rs|pq), but it may not have the other
**    left-pair and right-pair permutational symmetries (the CI twopdm
**    for example does not have this symmetry naturally).  The transform
**    needs this symmetry (as does the derivative program), so we will
**    enforce it here if the user specifies.  We currently assume that
**    the tpdm on disk has only unique (pq|rs) pairs, and so for our
**    purposes we need to generate (rs|pq) just like we do when reading
**    the AO integrals in the regular forwards transform.
**
** Based on the YOSH_RDTWO34() function
** C. David Sherrill
** Created August 1997
**
** Arguments:
**   YBuff        = Yoshimine object pointer
**   tei_file     = unit number for two-electron integrals
**   ioff         = standard lexical index array
**   symmetrize   = symmetrize the incoming 2pdm
**   add_ref_pt   = Add the factors arising from a reference determinant
**                  (n.b. assumes lowest MO's occupied)
**   del_tei_file = 1 to delete the tei file, 0 otherwise
**   printflag    = 1 for printing (for debugging only!) else 0
**   outfile      = file to print integrals to (if printflag is set)
*/
void yosh_rdtwo_backtr(struct yoshimine *YBuff, int tei_file, int *ioff,
                       int symmetrize, int add_ref_pt, int del_tei_file,
                       int printflag, std::string OutFileRMR)
{
  int i, ij, kl, ijkl;
  int iabs, jabs, kabs, labs;
  double value;
  struct bucket *bptr;
  int whichbucket, lastbuf, idx;
  long int tmpi;
  struct iwlbuf Inbuf;
  Value *valptr;
  Label *lblptr;

  if (printflag) {
    outfile->Printf( "Yoshimine rdtwo_backtr routine entered\n");
    outfile->Printf( "Two-particle density from file %d:\n", tei_file);
  }

  iwl_buf_init(&Inbuf, tei_file, YBuff->cutoff, 1, 0);
  lblptr = Inbuf.labels;
  valptr = Inbuf.values;

  do {
    iwl_buf_fetch(&Inbuf);
    lastbuf = Inbuf.lastbuf;
    for (idx=4*Inbuf.idx; Inbuf.idx < Inbuf.inbuf; Inbuf.idx++) {
      iabs = (int) lblptr[idx++];
      jabs = (int) lblptr[idx++];
      kabs = (int) lblptr[idx++];
      labs = (int) lblptr[idx++];

      iabs = moinfo.corr2pitz_nofzv[iabs];
      jabs = moinfo.corr2pitz_nofzv[jabs];
      kabs = moinfo.corr2pitz_nofzv[kabs];
      labs = moinfo.corr2pitz_nofzv[labs];

      value = valptr[Inbuf.idx];
      if (symmetrize) {
        if (iabs != jabs) value *= 0.5;
        if (kabs != labs) value *= 0.5;
      }

      /* calculate ijkl lexical index (make no i>=j assumptions for now) */
      ij = INDEX(iabs,jabs);
      kl = INDEX(kabs,labs);
      ijkl = INDEX(ij,kl);  /* ijkl needed only in the print routine */

      /* figure out what bucket to put it in, and do so */

      whichbucket = YBuff->bucket_for_pq[ij];
      bptr = YBuff->buckets+whichbucket;
      tmpi = (bptr->in_bucket)++;
      bptr->p[tmpi] = iabs;
      bptr->q[tmpi] = jabs;
      bptr->r[tmpi] = kabs;
      bptr->s[tmpi] = labs;
      bptr->val[tmpi] = value;

      if (printflag)
        outfile->Printf( "%4d %4d %4d %4d  %4d   %10.6lf\n",
                iabs, jabs, kabs, labs, ijkl, value) ;
      if ((tmpi+1) == YBuff->bucketsize) { /* need to flush bucket to disk */
        flush_bucket(bptr, 0);
        bptr->in_bucket = 0;
      }


      /* this generates (kl|ij) from (ij|kl) if necessary and puts it out */
      /* anaglogous to the "matrix" option in yosh_rdtwo()              */
      if (iabs != kabs || jabs != labs) {
        whichbucket = YBuff->bucket_for_pq[kl];
        bptr = YBuff->buckets+whichbucket;
        tmpi = (bptr->in_bucket)++;
        bptr->p[tmpi] = kabs;
        bptr->q[tmpi] = labs;
        bptr->r[tmpi] = iabs;
        bptr->s[tmpi] = jabs;
        bptr->val[tmpi] = value;
        if (printflag)
          outfile->Printf( "%4d %4d %4d %4d  %4d   %10.6lf\n",
                  kabs, labs, iabs, jabs, ijkl, value) ;
        if ((tmpi+1) == YBuff->bucketsize) {
          flush_bucket(bptr, 0);
          bptr->in_bucket = 0;
        }
      }

    }

  } while(!lastbuf);

  /* now add in contributions from the reference determinant if requested */
  if (add_ref_pt) add_2pdm_ref_pt(YBuff, ioff, printflag, "outfile");

  /* flush partially filled buckets */
  for (i=0; i<YBuff->nbuckets; i++) {
    flush_bucket((YBuff->buckets)+i, 1);
  }

  iwl_buf_close(&Inbuf, !del_tei_file);
}

/*
** YOSH_RDTWO_BACKTR_UHF() : Read two-particle density elements from
** an IWL file and prepare them for Yoshimine sorting.  The sorted
** twopdm elements, G(pqrs), produced by this code have unique p-q and
** r-s combinations, but no pq-rs packing.  However, the input
** twopdm's may not have this same structure, so the boolean arguments
** swap_bk and symm_pq are used to correct this problem.  There are
** two circumstances to consider:
**
** (1) If the MO twopdm lacks pq-rs symmetry (e.g., the AB twopdm),
** its input file should include all rs for each pq.  The swap_bk flag
** should be set to "0" in this case so that no "extra" G(rs,pq)
** components are written to the sorting buffers.  If the input twopdm
** has pq-rs symmetry (e.g., the AA and BB twopdms), only unique pq-rs
** combinations should be included and the swap_bk flag should be set
** to "1".
**
** (2) If the MO density lacks p-q and r-s symmetry, the input file
** should include all combinations of p,q and r,s, and the symm_pq
** flag should be set to "1".  If the MO density has p-q and r-s
** symmetry, then its input file should include only unique p,q and
** r,s combinations, and the symm_pq flag should be set to "0".
**
** Note that "intermediate" symmetry cases, where the MO twopdm has
** p-q symmetry but not r-s symmetry, for example, are not included
** here.
**
** Also, this code assumes the input indices are in QTS ordering and
** converts them automatically to Pitzer.
**
** TDC, 1/03
**
** Based on the YOSH_RDTWO_BACKTR() function above by
** C. David Sherrill
**
** Arguments:
**   YBuff        = Yoshimine object pointer
**   tei_file     = unit number for two-electron integrals
**   ioff         = standard lexical index array
**   swap_bk      = sort both G(pq,rs) and G(rs,pq) combinations
**   symm_pq      = symmetrize both p,q and r,s combinations
**   del_tei_file = 1 to delete the tei file, 0 otherwise
**   printflag    = 1 for printing (for debugging only!) else 0
**   outfile      = file to print integrals to (if printflag is set)
*/
void yosh_rdtwo_backtr_uhf(std::string spin, struct yoshimine *YBuff, int tei_file, int *ioff,
                           int swap_bk, int symm_pq, int del_tei_file,
                           int printflag, std::string OutFileRMR)
{
  int i, ij, kl, ijkl;
  int iabs, jabs, kabs, labs;
  double value;
  struct bucket *bptr;
  int whichbucket, lastbuf, idx;
  long int tmpi;
  struct iwlbuf Inbuf;
  Value *valptr;
  Label *lblptr;
  int *iorder, *jorder, *korder, *lorder;

  if (printflag) {
    outfile->Printf( "Yoshimine rdtwo_backtr routine entered\n");
    outfile->Printf( "Two-particle density from file %d:\n", tei_file);
  }

  if(spin == "AA") {
    iorder = moinfo.corr2pitz_nofzv_a;
    jorder = moinfo.corr2pitz_nofzv_a;
    korder = moinfo.corr2pitz_nofzv_a;
    lorder = moinfo.corr2pitz_nofzv_a;
  }
  else if(spin == "BB") {
    iorder = moinfo.corr2pitz_nofzv_b;
    jorder = moinfo.corr2pitz_nofzv_b;
    korder = moinfo.corr2pitz_nofzv_b;
    lorder = moinfo.corr2pitz_nofzv_b;
  }
  else if(spin == "AB") {
    iorder = moinfo.corr2pitz_nofzv_a;
    jorder = moinfo.corr2pitz_nofzv_a;
    korder = moinfo.corr2pitz_nofzv_b;
    lorder = moinfo.corr2pitz_nofzv_b;
  }
  else {
    outfile->Printf( "\n\tInvalid spin cases requested for backtransformation!\n");
    exit(PSI_RETURN_FAILURE);
  }

  iwl_buf_init(&Inbuf, tei_file, YBuff->cutoff, 1, 0);
  lblptr = Inbuf.labels;
  valptr = Inbuf.values;

  do {
    iwl_buf_fetch(&Inbuf);
    lastbuf = Inbuf.lastbuf;
    for (idx=4*Inbuf.idx; Inbuf.idx < Inbuf.inbuf; Inbuf.idx++) {
      iabs = (int) lblptr[idx++];
      jabs = (int) lblptr[idx++];
      kabs = (int) lblptr[idx++];
      labs = (int) lblptr[idx++];

      iabs = iorder[iabs];
      jabs = jorder[jabs];
      kabs = korder[kabs];
      labs = lorder[labs];

      value = valptr[Inbuf.idx];

      if (symm_pq) {
        if (iabs != jabs) value *= 0.5;
        if (kabs != labs) value *= 0.5;
      }

      /* calculate ijkl lexical index (make no i>=j assumptions for now) */
      ij = INDEX(iabs,jabs);
      kl = INDEX(kabs,labs);
      ijkl = INDEX(ij,kl);  /* ijkl needed only in the print routine */

      /* figure out what bucket to put it in, and do so */

      whichbucket = YBuff->bucket_for_pq[ij];
      bptr = YBuff->buckets+whichbucket;
      tmpi = (bptr->in_bucket)++;
      bptr->p[tmpi] = iabs;
      bptr->q[tmpi] = jabs;
      bptr->r[tmpi] = kabs;
      bptr->s[tmpi] = labs;
      bptr->val[tmpi] = value;

      if (printflag)
        outfile->Printf( "%4d %4d %4d %4d  %4d   %10.6lf\n",
                iabs, jabs, kabs, labs, ijkl, value) ;
      if ((tmpi+1) == YBuff->bucketsize) { /* need to flush bucket to disk */
        flush_bucket(bptr, 0);
        bptr->in_bucket = 0;
      }

      /* this generates (kl|ij) from (ij|kl) if necessary and puts it out */
      /* anaglogous to the "matrix" option in yosh_rdtwo()              */
      if (swap_bk && (iabs != kabs || jabs != labs)) {
        whichbucket = YBuff->bucket_for_pq[kl];
        bptr = YBuff->buckets+whichbucket;
        tmpi = (bptr->in_bucket)++;
        bptr->p[tmpi] = kabs;
        bptr->q[tmpi] = labs;
        bptr->r[tmpi] = iabs;
        bptr->s[tmpi] = jabs;
        bptr->val[tmpi] = value;
        if (printflag)
          outfile->Printf( "%4d %4d %4d %4d  %4d   %10.6lf\n",
                  kabs, labs, iabs, jabs, ijkl, value) ;
        if ((tmpi+1) == YBuff->bucketsize) {
          flush_bucket(bptr, 0);
          bptr->in_bucket = 0;
        }
      }

    }

  } while(!lastbuf);

  /* flush partially filled buckets */
  for (i=0; i<YBuff->nbuckets; i++) {
    flush_bucket((YBuff->buckets)+i, 1);
  }

  iwl_buf_close(&Inbuf, !del_tei_file);
}

/*
** ADD_2PDM_REF_PT
**
** This function adds in the contributions to the two-particle density
** matrix from the reference determinant, as might be required in MBPT
** or CC theory.  Assume the reference is made up of the lowest-lying
** orbitals as specified by DOCC and SOCC arrays.  Assume restricted
** orbitals.   Assume correlated ordering is docc, socc, virt...may
** not always be true, but order array is already wiped out by this
** point for backtransforms, would have to fetch it again.
**
** David Sherrill, Feb 1998
*/
void add_2pdm_ref_pt(struct yoshimine *YBuff, int *ioff, int pflg,
                     std::string OutFileRMR)
{
   int i, j, ii, jj, ij;
   int iabs, jabs;
   int ndocc, nsocc;

   ndocc = moinfo.ndocc;  nsocc = moinfo.nsocc;

   /* closed-shell part */
   for (i=0; i<ndocc; i++) {
     iabs = moinfo.corr2pitz_nofzv[i];
     ii = ioff[iabs] + iabs;

     for (j=0; j<i; j++) {
       jabs = moinfo.corr2pitz_nofzv[j];
       jj = ioff[jabs] + jabs;
       ij = ioff[iabs] + jabs;

       yosh_buff_put_val(YBuff,ioff,ii,iabs,iabs,jabs,jabs, 2.0,pflg,"outfile");
       yosh_buff_put_val(YBuff,ioff,jj,jabs,jabs,iabs,iabs, 2.0,pflg,"outfile");
       yosh_buff_put_val(YBuff,ioff,ij,iabs,jabs,jabs,iabs,-0.25,pflg,"outfile");
       yosh_buff_put_val(YBuff,ioff,ij,jabs,iabs,iabs,jabs,-0.25,pflg,"outfile");

     }

     jabs = moinfo.corr2pitz_nofzv[j];
     jj = ioff[jabs] + jabs;
     ij = ioff[iabs] + jabs;
     yosh_buff_put_val(YBuff,ioff,ii,iabs,iabs,jabs,jabs, 1.0,pflg,"outfile");

   }


   /* open-shell part, if any */
   for (i=ndocc; i<ndocc+nsocc; i++) {
     iabs = moinfo.corr2pitz_nofzv[i];
     ii = ioff[iabs] + iabs;

     for (j=0; j<ndocc; j++) {
       jabs = moinfo.corr2pitz_nofzv[j];
       jj = ioff[jabs] + jabs;
       ij = ioff[iabs] + jabs;

       yosh_buff_put_val(YBuff,ioff,ii,iabs,iabs,jabs,jabs, 1.0,pflg,"outfile");
       yosh_buff_put_val(YBuff,ioff,jj,jabs,jabs,iabs,iabs, 1.0,pflg,"outfile");
       yosh_buff_put_val(YBuff,ioff,ij,iabs,jabs,jabs,iabs,-0.125,pflg,"outfile");
       yosh_buff_put_val(YBuff,ioff,ij,jabs,iabs,iabs,jabs,-0.125,pflg,"outfile");
     }

     for (j=ndocc; j<i; j++) {
       jabs = moinfo.corr2pitz_nofzv[j];
       jj = ioff[jabs] + jabs;
       ij = ioff[iabs] + jabs;

       yosh_buff_put_val(YBuff,ioff,ii,iabs,iabs,jabs,jabs, 0.5,pflg,"outfile");
       yosh_buff_put_val(YBuff,ioff,jj,jabs,jabs,iabs,iabs, 0.5,pflg,"outfile");
       yosh_buff_put_val(YBuff,ioff,ij,iabs,jabs,jabs,iabs,-0.125,pflg,"outfile");
       yosh_buff_put_val(YBuff,ioff,ij,jabs,iabs,iabs,jabs,-0.125,pflg,"outfile");
     }

   }

}



/*
** YOSH_BUFF_PUT_VAL
**
** This function puts a value and its associated indices to a yoshimine
** sorting buffer.
*/
void yosh_buff_put_val(struct yoshimine *YBuff, int *ioff, int pq,
                       int p, int q, int r, int s, double value, int prtflg,
                       std::string OutFileRMR)
{
   struct bucket *bptr;
   int whichbucket;
   long int tmpi;

   whichbucket = YBuff->bucket_for_pq[pq];
   bptr = YBuff->buckets+whichbucket;
   tmpi = (bptr->in_bucket)++;
   bptr->p[tmpi] = p;
   bptr->q[tmpi] = q;
   bptr->r[tmpi] = r;
   bptr->s[tmpi] = s;
   bptr->val[tmpi] = value;

   if (prtflg)
     outfile->Printf( "%4d %4d %4d %4d         %10.6lf\n", p, q, r, s,
             value);

   if ((tmpi+1) == YBuff->bucketsize) { /* need to flush bucket to disk */
     flush_bucket(bptr, 0);
     bptr->in_bucket = 0;
   }


}



/*
** YOSH_SORT(): Sort all the buckets in the Yoshimine sorting algorithm.
**    The call to sortbuf() will cause the intermediate files to be
**    deleted unless keep_bins is set to 1.
**
** Arguments:
**    YBuff        =  pointer to Yoshimine object
**    out_tape     =  number for binary output file
**    keep_bins    =  keep the intermediate tmp files
**    ioff         =  the usual offset array unless no_pq_perm, in which
**                    case it is appropriate for the left indices, i.e.,
**                    pq = ioff[p] + q
**    ioff2        =  the Elbert ioff2 array if elbert=true.  If no_pq_perm,
**                    then this is the usual ioff array, appropriate
**                    for the right indices
**    nbfso        =  number of basis fns in symmetry orbitals
**    ket_indices  =  number of ket indices (usually ntri)
**    elbert       =  are inputs in Elbert ordering? p>=q, r>=s, rs>=pq
**    intermediate =  1 if sorting an intermediate in the transformation
**                    (argument to sortbuf()).  This implies that the full
**                    set of rs are available for a given pq.
**    no_pq_perm   =  if p and q are not interchangeable (e.g., one is
**                    occupied and one is virtual, as in MP2)
**    qdim         =  the number of possible values for index q
**    add          =  do additions of integrals during the sort?
**    print_lvl    =  verbosity level (how much to print)
**    outfile      =  text output file
**
** Returns: none
*/
void yosh_sort(struct yoshimine *YBuff, int out_tape, int keep_bins,
      int *ioff, int *ioff2, int nbfso, int ket_indices, int elbert,
      int intermediate, int no_pq_perm, int qdim, int add, int print_lvl,
      std::string OutFileRMR)
{
   double *twoel_ints;
   int i, max_pq;
   struct iwlbuf inbuf, outbuf;

   /* may be slightly more than pq_per_bucket pq's in each bucket
    * if the pq's didn't divide evenly among the buckets.  The remainder
    * will go to the last bucket.
    */
   max_pq = YBuff->buckets[YBuff->core_loads-1].hi -
            YBuff->buckets[YBuff->core_loads-1].lo
            + 1;

   twoel_ints = init_array(max_pq * ket_indices) ;
   iwl_buf_init(&outbuf, out_tape, YBuff->cutoff, 0, 0);

   for (i=0; i<YBuff->core_loads-1; i++) {
      if (print_lvl > 1) outfile->Printf( "Sorting bin %d\n", i+1);
      iwl_buf_init(&inbuf, YBuff->first_tmp_file+i, YBuff->cutoff, 1, 0);
      sortbuf(&inbuf, &outbuf, twoel_ints, (YBuff->buckets)[i].lo,
              (YBuff->buckets)[i].hi, ioff, ioff2, nbfso, elbert,
               intermediate, no_pq_perm, qdim, add, (print_lvl > 4), "outfile");
      zero_arr(twoel_ints, max_pq * ket_indices);
      /* zero_arr(twoel_ints, YBuff->pq_per_bucket * YBuff->bra_indices); */
      iwl_buf_close(&inbuf, keep_bins);
      }


   if (print_lvl > 1) outfile->Printf( "Sorting bin %d\n", i+1) ;
   iwl_buf_init(&inbuf, YBuff->first_tmp_file+i, YBuff->cutoff, 1, 0);
   sortbuf(&inbuf, &outbuf, twoel_ints, (YBuff->buckets)[i].lo,
           (YBuff->buckets)[i].hi, ioff, ioff2, nbfso, elbert,
           intermediate, no_pq_perm, qdim, add, (print_lvl > 4), "outfile");
   iwl_buf_close(&inbuf, keep_bins);

   if (print_lvl > 1) outfile->Printf( "Done sorting.\n");

   iwl_buf_flush(&outbuf, 1);
   iwl_buf_close(&outbuf, 1);
   free(twoel_ints);
}

/*
** YOSH_SORT_PK(): Sort all the buckets in the Yoshimine sorting algorithm.
**    The call to sortbuf_pk() will cause the intermediate files to be
**    deleted unless keep_bins is set to 1.
**    The sorting for Coulomb and exchange integrals is not exactly the same, thus the
**    flag is_exch has to be set to true for the exchange matrix.
**    The integrals are ordered as ij >= kl; i >= j and k >= l, with appropriate
**    factors for diagonal elements. In addition, integrals are written directly
**    to the PK file without labels.
**
** Arguments:
**    YBuff        =  pointer to Yoshimine object
**    is_exch      =  flag set to true if we are sorting for
**                    the exchange matrix
**    out_tape     =  number for binary output file
**    keep_bins    =  keep the intermediate tmp files
**    so2rel       =  array mapping absolute basis function index to relative
**                    basis function index within an irrep, so2rel[abs] = rel
**    so2sym       =  array mapping absolute basis function index to irrep
**                    number, so2sym[abs] = sym
**    pksymoff     =  array containing the offset in each irrep to convert a
**                    pq index computed with relative indices to an absolute
**                    pq index, pqrel = ioff[prel] + qrel, pqabs = pqrel + pksymoff[psym]
**    ioff         =  the usual offset array
**    num_so       =  number of basis fns per irrep
**    ket_indices  =  number of ket indices (usually ntri)
**    qdim         =  the number of possible values for index q
**    print_lvl    =  verbosity level (how much to print)
**    outfile      =  text output file
**
** Returns: none
*/
void yosh_sort_pk(struct yoshimine *YBuff, int is_exch, int out_tape, int keep_bins,
      int* so2ind, int* so2sym, int* pksymoff, int* ioff, int print_lvl)
{
   size_t batch_size = 0;
   size_t nintegrals;
   double *twoel_ints;
   int lopq, hipq, i;
   boost::shared_ptr<PSIO> psio = _default_psio_lib_;
   char* label = new char[100];
   struct iwlbuf inbuf;

   // We compute the maximum batch size
   for(int i = 0; i < YBuff->nbuckets; ++i) {
       lopq = YBuff->buckets[i].lo;
       hipq = YBuff->buckets[i].hi + 1;
       nintegrals = ioff[hipq] - ioff[lopq];
       if (nintegrals > batch_size) batch_size = nintegrals;
   }

   twoel_ints = init_array(batch_size);

   for (i = 0; i < YBuff->core_loads; i++) {
      if (print_lvl > 1) outfile->Printf( "Sorting bin %d\n", i+1);
      lopq = YBuff->buckets[i].lo;
      hipq = YBuff->buckets[i].hi;
      iwl_buf_init(&inbuf, YBuff->first_tmp_file+i, YBuff->cutoff, 1, 0);
      sortbuf_pk(&inbuf, out_tape, is_exch, twoel_ints, lopq,
              hipq, so2ind, so2sym, pksymoff, ioff, (print_lvl > 4), "outfile");
      // Since everything is in triangle form, we can totally get the size
      nintegrals = ioff[hipq + 1] - ioff[lopq];
      outfile->Printf("Batch number %i, nintegrals is %i\n", i, nintegrals);
      if (is_exch) {
        sprintf(label,"K Block (Batch %d)", i);
      } else {
        sprintf(label,"J Block (Batch %d)", i);
      }
      psio->write_entry(out_tape, label, (char*)twoel_ints, nintegrals * sizeof(double));
      if(i < YBuff->core_loads - 1) {
          outfile->Printf("Zeroing twoelints\n");
        zero_arr(twoel_ints, batch_size);
      }
      iwl_buf_close(&inbuf, keep_bins);
   }


   if (print_lvl > 1) outfile->Printf( "Done sorting.\n");

   free(twoel_ints);
   delete [] label;
}



/*
** YOSH_FREE(): Free up a Yoshimine object.  Free any dynamically-allocated
**    memory.
*/
void yosh_free(struct yoshimine *YBuff)
{
   free(YBuff->buckets);
}


/*
** YOSH_FLUSH(): Flush any of the buckets and tag them as 'last buffer'
**    for each bucket tmp file.
**
*/
void yosh_flush(struct yoshimine *YBuff)
{
   int i;

   for (i=0; i<YBuff->nbuckets; i++) {
         flush_bucket((YBuff->buckets)+i, 1);
      }
}



/*
** FLUSH_BUCKET():  This function flushes a Yoshimine bucket to a temporary
**    binary file.  Must be careful to call with lastbuf==1 when the
**    last buffer has been reached, so that the iwl buffer can set the
**    `last buffer' flag on output.
*/
void flush_bucket(struct bucket *bptr, int lastbuf)
{

   iwl_buf_wrt_arr(&(bptr->IWLBuf), bptr->val, bptr->p, bptr->q,
     bptr->r, bptr->s, bptr->in_bucket);
   iwl_buf_flush(&(bptr->IWLBuf), lastbuf);
}


/*
** YOSH_WRT_ARR(): Write an array to a Yoshimine object...useful for writing
**    intermediates to disk from a transformation program.  Write only
**    nonzero elements so that less disk space is required, and so that
**    zeroes for i,j,k,l and value for the first element in a buffer can
**    indicate an empty buffer.  (Eventually I should probably use
**    the second byte of the flags field in each buffer to indicate a
**    zero buffer, else use the second and third bytes to denote the
**    number of integrals in each buffer.  File34 had a field like
**    this, but it assumes an integral number of integer words per
**    double, which might not be true in the long run.  IWL is currently
**    free of this assumption).
**
**    The sortbuf() routine wants to have blocks containing all rs for
**    a given pq.  If instead, we want all pq for a given rs, we
**    need to reverse the order of pq and rs in the tmp files.  Set
**    sortby_rs=1.
**
** Arguments:
**    YBuff     =  pointer to Yoshimine object
**    p         =  common p value for array
**    q         =  common q value for array
**    pq        =  compound index determined from p and q
**    pqsym     =  the direct product of the symmetries of p and q
**    arr       =  the array containing the data for a given pq
**    rmax      =  the maximum value of the third index
**    ioff      =  the standard offset array
**    orbsym    =  the irrep for each orbital
**    firsti    =  the first orbital for each irrep
**    lasti     =  last orbital for each irrep
**    sortby_rs =  described above
**    printflag =  verbosity level for printing
**    outfile   =  text output file
**
*/
void yosh_wrt_arr(struct yoshimine *YBuff, int p, int q, int pq, int pqsym,
   double *arr, int rmax, int *ioff, int *orbsym, int *firsti,
   int *lasti, int sortby_rs, int printflag, std::string OutFileRMR)
{
   int r, s, rs, rsym, ssym, smax;
   long int tmpi;
   int whichbucket;
   int lastflag = 0, firstfile;
   struct bucket *bptr;
   double value;

   if (printflag) {
     outfile->Printf( "\nyosh_wrt_arr called for p=%d,q=%d\n", p, q);
   }

   firstfile = YBuff->first_tmp_file;

   for (r=0; r<rmax; r++) {
      rsym = orbsym[r];
      ssym = pqsym ^ rsym;
      if (ssym > rsym) continue;

      smax = (rsym == ssym) ? r : lasti[ssym];

      for (s=firsti[ssym]; s<=smax; s++) {
         rs = ioff[r] + s;
         value = arr[rs];

         if (fabs(value) > YBuff->cutoff) {
            /* figure out which bucket the integral should go in */
            if (sortby_rs) whichbucket = YBuff->bucket_for_pq[rs];
            else whichbucket = YBuff->bucket_for_pq[pq];

            bptr = YBuff->buckets+whichbucket ;
            tmpi = (bptr->in_bucket)++ ;

            /* sortbuf wants to sort by the first index.  If we want to
             * sort by rs, we need to map r->p, s->q, p->r, q->s
             */
            if (sortby_rs) {
               bptr->p[tmpi] = r;
               bptr->q[tmpi] = s;
               bptr->r[tmpi] = p;
               bptr->s[tmpi] = q;
               }
            else {
               bptr->p[tmpi] = p;
               bptr->q[tmpi] = q;
               bptr->r[tmpi] = r;
               bptr->s[tmpi] = s;
               }

            bptr->val[tmpi] = value;

            if (printflag)
               outfile->Printf( "%4d %4d %4d %4d  %10.6lf\n",
                  p, q, r, s, arr[rs]);

            /* if we need to flush bucket to disk */
            if ((tmpi+1) == YBuff->bucketsize) {
               flush_bucket(bptr, lastflag);
               bptr->in_bucket = 0 ;
               }

            /*
             * What if this was the last buffer, and we didn't set the lastbuf
             * flag?  Must write a buffer of zeroes using the yosh_flush()
             * routine.  Must be careful to always have the last buffer
             * set the lastbuf flag (to avoid read error) and must also be
             * sure that an exactly-filled buffer is written out.  Carelessness
             * might cause a failure to write the last buffer if it has
             * exactly 'bucketsize' elements. ---CDS
             */

            } /* end if (fabs(arr[rs])) ... */

         } /* end loop over s */

      }  /* end loop over r */



}


/*
** YOSH_WRT_ARR2(): Write an array to a Yoshimine object (nonzero
**    elements only).  The values of the integrals are in array 'arr',
**    and the indices are input as follows:  all integrals have the
**    common indices 'p' and 'q'.  The remaining indices are input
**    in arrays 'rlist' and 'slist'.  The number of integrals in the
**    current block is given by the parameter 'size'.  As currently
**    written, the routine will write out the integrals in canonical
**    ordering p>=q, r>=s, pq>=rs.  It might conceivably be useful
**    in the future to use a 'sortby_rs' flag as in yosh_wrt_arr().
**
** Arguments:
**    YBuff     =  pointer to Yoshimine object
**    size      =  number of integrals to process in the array
**    arr       =  the array containing the data for a given pq
**    p         =  common p value for array
**    q         =  common q value for array
**    rlist     =  list of indices r
**    slist     =  list of indices s (no guarantee rlist[i] >= slist[i])
**    ioff      =  the standard offset array
**    printflag =  verbosity level for printing
**    outfile   =  text output file
**
*/
void yosh_wrt_arr2(struct yoshimine *YBuff, int size, double *arr,
   int p, int q, int *rlist, int *slist, int *ioff, int printflag,
   std::string OutFileRMR)
{
   int x;
   int i1, j1, k1, l1, ij1, kl1;
   int ktmp, ltmp;
   int i2, j2, k2, l2, ij2, kl2;
   long int tmpi;
   int whichbucket;
   int lastflag = 0, firstfile;
   struct bucket *bptr;
   double value;

   firstfile = YBuff->first_tmp_file;

   i1 = MAX0(p,q);
   j1 = MIN0(p,q);
   ij1 = ioff[i1] + j1;

   for (x=0; x<size; x++) {

      value = *arr++;
      ktmp = *rlist++;
      ltmp = *slist++;

      if (fabs(value) < YBuff->cutoff) continue;

      k1 = MAX0(ktmp,ltmp);
      l1 = MIN0(ktmp,ltmp);
      kl1 = ioff[k1] + l1;

      if (ij1 < kl1) {
         i2 = k1;
         j2 = l1;
         k2 = i1;
         l2 = j1;
         ij2 = kl1;
         kl2 = ij1;
         }
      else {
         i2 = i1;
         j2 = j1;
         k2 = k1;
         l2 = l1;
         ij2 = ij1;
         kl2 = kl1;
         }

   whichbucket = YBuff->bucket_for_pq[ij2];
   bptr = YBuff->buckets+whichbucket ;
   tmpi = (bptr->in_bucket)++ ;

   bptr->p[tmpi] = i2;
   bptr->q[tmpi] = j2;
   bptr->r[tmpi] = k2;
   bptr->s[tmpi] = l2;

   bptr->val[tmpi] = value;

   if (printflag) outfile->Printf( "%4d %4d %4d %4d  %10.6lf\n",
         i2, j2, k2, l2, value);

   /* if we need to flush bucket to disk */
      if ((tmpi+1) == YBuff->bucketsize) {
         flush_bucket(bptr, lastflag);
         bptr->in_bucket = 0;
         }

   } /* end loop over x */

}

/*
** YOSH_WRT_ARR_MP2(): Write an array to a Yoshimine object...useful for
**    writing intermediates to disk from a transformation program.  Write
**    only nonzero elements so that less disk space is required, and so that
**    zeroes for i,j,k,l and value for the first element in a buffer can
**    indicate an empty buffer.
**
**    The sortbuf() routine wants to have blocks containing all rs for
**    a given pq.  If instead, we want all pq for a given rs, we
**    need to reverse the order of pq and rs in the tmp files.  Set
**    sortby_rs=1.
**
**    This routine has been altered to work _only_ for MP2 restricted
**    sorts.  It expects the r-index to correlate to occupied orbitals and
**    the s-index to virtual orbitals.  This routine may not, in general
**    be used with other orbital distributions.
**    -Daniel 9/20/95
**
** Arguments:
**    YBuff     =  pointer to Yoshimine object
**    p         =  common p value for array
**    q         =  common q value for array
**    pq        =  compound index determined from p and q
**    pqsym     =  the direct product of the symmetries of p and q
**    arr       =  the array containing the data for a given pq
**    rsym      =  the irrep of the r-indices
**    firstr    =  the first r-index for each irrep
**    lastr     =  the last r-index for each irrep
**    firsts    =  the first s-index for each irrep
**    lasts     =  the last s-index for each irrep
**    sortby_rs =  described above
**    ndocc     =  the number of doubly-occupied orbitals
**    nvirt     =  the number of virtual orbitals
**    occ       =  the Pitzer -> QTS ordering array for the occupied orbitals
**    vir       =  the Pitzer -> QTS ordering array for the virtual orbitals
**    ioff3     =  offset array for ndocc*nvirt arrays
**    printflag =  verbosity level for printing
**    outfile   =  text output file
**
*/
void yosh_wrt_arr_mp2(struct yoshimine *YBuff, int p, int q, int pq,
                      int pqsym, double **arr, int rsym, int *firstr,
                      int *lastr, int *firsts, int *lasts, int sortby_rs,
                      int ndocc, int nvirt, int *occ, int *vir, int *ioff3,
                      int printflag, std::string OutFileRMR)
{
   int r, s, rs, ssym;
   int rnew,snew;
   int R,S;
   long int tmpi;
   int whichbucket;
   int lastflag = 0, firstfile;
   struct bucket *bptr;
   double value;

   firstfile = YBuff->first_tmp_file;
   ssym = pqsym ^ rsym;
   for (r=firstr[rsym], R=0; r <= lastr[rsym]; r++, R++) {
       rnew = occ[r];
      for (s=firsts[ssym], S=0; s <=lasts[ssym]; s++, S++) {
          snew = vir[s];
          rs = ioff3[rnew] + snew;
          value = arr[R][S];

         if (fabs(value) > YBuff->cutoff) {
            /* figure out which bucket the integral should go in */
            if (sortby_rs) whichbucket = YBuff->bucket_for_pq[rs];
            else whichbucket = YBuff->bucket_for_pq[pq];

            bptr = YBuff->buckets+whichbucket ;
            tmpi = (bptr->in_bucket)++ ;

            /* sortbuf wants to sort by the first index.  If we want to
             * sort by rs, we need to map r->p, s->q, p->r, q->s
             * We also write out the QTS ordered indices for r and s, and
             * not the Pitzer indices.
             */
            if (sortby_rs) {
               bptr->p[tmpi] = rnew;
               bptr->q[tmpi] = snew;
               bptr->r[tmpi] = p;
               bptr->s[tmpi] = q;
               }
            else {
               bptr->p[tmpi] = p;
               bptr->q[tmpi] = q;
               bptr->r[tmpi] = rnew;
               bptr->s[tmpi] = snew;
               }

            bptr->val[tmpi] = value;

            if (printflag)
               outfile->Printf( "%4d %4d %4d %4d  %10.6lf\n",
                  p, q, r, s, value);

            /* if we need to flush bucket to disk */
            if ((tmpi+1) == YBuff->bucketsize) {
               flush_bucket(bptr, lastflag);
               bptr->in_bucket = 0;
               }

            } /* end if (fabs(arr[rs])) ... */

         } /* end loop over s */

      }  /* end loop over r */

}



/*
** YOSH_WRT_ARR_MP2R12A(): Write an array to a Yoshimine object...useful for
**    writing intermediates to disk from a transformation program.  Write
**    only nonzero elements so that less disk space is required, and so that
**    zeroes for i,j,k,l and value for the first element in a buffer can
**    indicate an empty buffer.
**
**    The sortbuf() routine wants to have blocks containing all rs for
**    a given pq.  If instead, we want all pq for a given rs, we
**    need to reverse the order of pq and rs in the tmp files.  Set
**    sortby_rs=1.
**
**    This routine has been altered to work _only_ for MP2R12A restricted
**    sorts.  It expects the r-index to correlate to occupied orbitals.
**    This routine may not, in general
**    be used with other orbital distributions.
**    -Edward
**
** Arguments:
**    YBuff     =  pointer to Yoshimine object
**    p         =  common p value for array
**    q         =  common q value for array
**    pq        =  compound index determined from p and q
**    pqsym     =  the direct product of the symmetries of p and q
**    arr       =  the array containing the data for a given pq
**    rsym      =  the irrep of the r-indices
**    firstr    =  the first r-index for each irrep
**    lastr     =  the last r-index for each irrep
**    firsts    =  the first s-index for each irrep
**    lasts     =  the last s-index for each irrep
**    sortby_rs =  described above
**    occ       =  the Pitzer -> QTS ordering array for the occupied orbitals
**    ioff3     =  offset array for ndocc*nmo arrays
**    printflag =  verbosity level for printing
**    outfile   =  text output file
**
*/
void yosh_wrt_arr_mp2r12a(struct yoshimine *YBuff, int p, int q, int pq,
                          int pqsym, double **arr, int rsym, int *firstr,
                          int *lastr, int *firsts, int *lasts, int sortby_rs,
                          int *occ, int *ioff3,
                          int printflag, std::string OutFileRMR)
{
   int r, s, rs, ssym;
   int rnew,snew;
   int R,S;
   long int tmpi;
   int whichbucket;
   int lastflag = 0, firstfile;
   struct bucket *bptr;
   double value;

   firstfile = YBuff->first_tmp_file;
   ssym = pqsym ^ rsym;
   for (r=firstr[rsym], R=0; r <= lastr[rsym]; r++, R++) {
       rnew = occ[r];
      for (s=firsts[ssym], S=0; s <=lasts[ssym]; s++, S++) {
          snew = s;
          rs = ioff3[rnew] + snew;
          value = arr[R][S];

         if (fabs(value) > YBuff->cutoff) {
            /* figure out which bucket the integral should go in */
            if (sortby_rs) whichbucket = YBuff->bucket_for_pq[rs];
            else whichbucket = YBuff->bucket_for_pq[pq];

            bptr = YBuff->buckets+whichbucket ;
            tmpi = (bptr->in_bucket)++ ;

            /* sortbuf wants to sort by the first index.  If we want to
             * sort by rs, we need to map r->p, s->q, p->r, q->s
             * We also write out the QTS ordered index for r, and
             * not the Pitzer index.
             */
            if (sortby_rs) {
               bptr->p[tmpi] = rnew;
               bptr->q[tmpi] = snew;
               bptr->r[tmpi] = p;
               bptr->s[tmpi] = q;
               }
            else {
               bptr->p[tmpi] = p;
               bptr->q[tmpi] = q;
               bptr->r[tmpi] = rnew;
               bptr->s[tmpi] = snew;
               }

            bptr->val[tmpi] = value;

            if (printflag)
               outfile->Printf( "%4d %4d %4d %4d  %10.6lf\n",
                  p, q, r, s, value);

            /* if we need to flush bucket to disk */
            if ((tmpi+1) == YBuff->bucketsize) {
               flush_bucket(bptr, lastflag);
               bptr->in_bucket = 0;
               }

            } /* end if (fabs(arr[rs])) ... */

         } /* end loop over s */

      }  /* end loop over r */
}

}} // end namespace psi::transqt

