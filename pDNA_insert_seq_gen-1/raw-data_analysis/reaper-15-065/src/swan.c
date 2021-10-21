
/*
 * (C) Copyright 2014, 2015 European Molecular Biology Laboratory.
 * Author: Stijn van Dongen <stijn@ebi.ac.uk>.
 * Contact: <kraken@ebi.ac.uk>
 *
 * This file is part of Reaper.   Reaper is free software: you can redistribute
 * it  and/or modify it under the terms of the  GNU  General  Public License as
 * published by the Free Software Foundation;  either version 3 of the License,
 * or  (at your option)  any later version.  This program is distributed in the
 * hope that it will be useful,  but  WITHOUT  ANY  WARRANTY;  without even the
 * implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 * See the  GNU  General  Public  License  for  more  details.  You should have
 * received a copy of the  GNU  General Public License along with this program.
 * If not, see http://www.gnu.org/licenses/.
*/


/* TODO: memory read_fasta_file, make_index
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdarg.h>
#include <ctype.h>
#include <errno.h>

#include "sw.h"
#include "slib.h"
#include "dna.h"



#define PRINT_KEY_VALUE 1
#define PRINT_EXCISE    2

static const char* g_trace = "n";
static int g_cell = 0;
static int g_mode_print = 0;
static int g_qlen = 0;
static int g_rlen = 0;
static int g_indel_allowed = 1;
static int g_dump_matrix = 0;
static int g_self = 0;


struct req
{  char* s
;  int slen
;  char* annot
;  unsigned busy
;
}  ;


struct kmer_x_req
{  unsigned  kmer
;  unsigned  req_index
;
}  ;


struct kmer_index
{  unsigned* kxr_offset
;  unsigned* kxr_length
;  struct kmer_x_req* kxr
;  unsigned k              /* defines size of the above */
;
}  ;


struct kmer_gen
{  unsigned cur_kmer
;  int last_N
;  int i
;  const char* s
;  unsigned slen
;  unsigned k
;
}  ;


void kmer_gen_init
(  struct kmer_gen* kg
,  const char* s
,  unsigned slen
,  unsigned k
)
   {  kg->cur_kmer = 0
   ;  kg->last_N   = -1
   ;  kg->i        = 0
   ;  kg->s        = s
   ;  kg->slen     = slen
   ;  kg->k        = k
;  }


unsigned* kmer_gen_next
(  struct kmer_gen* kg
)
   {  int i = kg->i
   ;  unsigned base = BASEMAP((unsigned char) kg->s[i])        /* will hit/need \0 */

   ;  if (i >= kg->slen)
      return NULL
   ;  kg->i++              /* note, i contains previous value */

   ;  kg->cur_kmer = (kg->cur_kmer << 2)
   ;  if (base == 4)
      kg->last_N = i
   ;  else if (base < 4)
      kg->cur_kmer |= base
   ;  kg->cur_kmer &= ((1 << (2*kg->k)) - 1)

   ;  if (kg->last_N + kg->k <= i)
      return &kg->cur_kmer
   ;  return kmer_gen_next(kg)
;  }


void sw_printaln5
(  struct sw_alninfo* ai
,  int accept
,  int score
,  void* fp
,  int zippit
,  int recno
)
   {  unsigned o = ai->aln_ofs
   ;  char space[512]
   ;  char buf[8192]
   ;  int n
   ;  int leftflush = ai->lft_start > ai->rgt_start
   ;  int shift = (int) ai->lft_start - (int) ai->rgt_start

   ;  memset(space, ' ', 512)
   ;  space[511] = '\0'

   ;  if (shift < 0)
      shift = -shift

   ;  n
   =  snprintf
      (  buf
      ,  8192
      ,  "%.*s%.*s%s%s : [%d,%d]\n%.*s%.*s%s : score %d\n%.*s%.*s%s%s : [%d,%d] %d/%d recno=%d\n.\n"
       /* 1 1 1 1 1 1     2  2    3 3 3 3 3          3   4 4 4 4 4 4    5  5   6  6      7 */
/* 1 */
      ,  leftflush ? 0 : shift         /* this much */
      ,  space
      ,  (int) (ai->lft_start-1)       /* this much */
      ,  ai->left
      ,  ai->aln_lft+o
      ,  ai->left + (ai->lft_end)
/* 2 */
      ,  (int) ai->lft_start
      ,  (int) ai->lft_end
/* 3 */
      ,  shift
      ,  space
      ,  (int) ((leftflush ? ai->rgt_start : ai->lft_start) - 1)
      ,  space
      ,  ai->aln_aln+o
      ,  (int) ai->data[ai->max_ij]
/* 4 */
      ,  leftflush ? shift : 0
      ,  space
      ,  (int) (ai->rgt_start -1)
      ,  ai->right
      ,  ai->aln_rgt+o
      ,  ai->right + (ai->rgt_end)
/* 5 */
      ,  (int) ai->rgt_start
      ,  (int) ai->rgt_end
/* 6 */
      ,  (int) accept
      ,  (int) score
/* 7 */
      ,  (int) recno
      )
#if WE_USE_ZLIB
   ;  if (n > 0)
      {
         if (zippit) gzwrite(fp, buf, n)
      ;  else        fputs(buf, fp)
#else
         fputs(buf, fp)
#endif
   ;  }
   }



   /* A stopgap solution.
   */
SWNUM* wrap_sw_fill2
(  struct sw_alninfo* ai
,  const char *left           /* fixme, length-encode */
,  const char *right          /* fixme, length-encode */
,  int indel_allowed
,  struct sw_param* swp
)
   {  SWNUM* allocdata  = NULL
   ;  unsigned leftlen  = strlen(left)
   ;  unsigned rightlen = strlen(right)
   ;  unsigned long thesize  = (leftlen+1) * (rightlen+1)
   ;  int ret = 0
   
   ;  allocdata = myalloc(thesize)

   ;  ret = sw_fill2(ai, allocdata, thesize, left, right, indel_allowed, swp)

   ;  if (ret)
      {  free(allocdata)
      ;  die(1, "memory allocation failure, most likely")
   ;  }
      return allocdata
;  }


int pp_cb                       /* pretty print callback */
(  const struct sw_alninfo* ai
,  char* buf
,  unsigned bufsize
)
   {  return snprintf(buf, bufsize, " ref-offset=%d query-offset=%d aln={%s}", ai->lft_start, ai->rgt_start, ai->aln_aln+ai->aln_ofs)
;  }


int do_align
(  FILE* fpo
,  const char* left
,  const char* annot_left
,  const char* right
,  const char* annot_right
,  struct sw_param* swp
,  int recno
)
   {
   ;  SWNUM* thedata = NULL
   ;  struct sw_alninfo ai = { 0 }
   ;  int cell = 0

   ;  thedata = wrap_sw_fill2(&ai, left, right, g_indel_allowed, swp)

   ;  if (g_dump_matrix)
      sw_dump(&ai)

   ;  cell = ai.max_ij
   ;  if (g_cell && (g_cell < ai.nj+1 || g_cell >= ai.ni * ai.nj))
      cell = g_cell

   ;  if (strchr(g_trace, 'o'))
      sw_trace2(&ai, swp, g_indel_allowed, g_cell)
   ;  else
      sw_trace2a
      (  &ai
      ,  swp
      ,  cell
      ,     (strchr(g_trace, 'l') ? SW_TRIMLEFT : 0)
         |  (strchr(g_trace, 'r') ? SW_TRIMRIGHT : 0)
      )

   ;  if (g_mode_print & PRINT_KEY_VALUE)
      sw_lp(annot_left, annot_right, cell, &ai, fpo, 0, recno, NULL)
   ;  else if (g_mode_print & PRINT_EXCISE)
      sw_pp_excise(annot_left, annot_right, cell, &ai, fpo, 0, recno, NULL)
   ;  else
      sw_pp2(annot_left, annot_right, cell, &ai, fpo, 0, recno, NULL)

   ;  if (thedata)
      free(thedata)

   ;  return 0
;  }


   /* returns -1 if id_threshold is != 0
   */
int search_best_match
(  FILE* fpo
,  struct kmer_index* kidx
,  struct sw_param *swp
,  struct req* theref            /* many */
,  struct req* thequery          /* one  */
,  int* overhang
,  unsigned* todo
,  unsigned n_todo
,  unsigned id_threshold
)
   {
#define MATRIX_SIZE 8192
   ;  struct sw_alninfo ai = { 0 }
   ;  struct kmer_gen kg
   ;  int best_i = 0, i
   ;  int best_identity = 0
   ;  int best_overhang = 1 << 20
   ;  unsigned* kmerp

                                    /* hierverder: check index offset/length correctness */
   ;  if (!n_todo)
      {  kmer_gen_init(&kg, thequery->s, strlen(thequery->s), kidx->k)
      ;  while ((kmerp = kmer_gen_next(&kg)))
         {  int i
         ;  unsigned thekmer = kmerp[0]
;if (0 && kg.i == kidx->k)
fprintf(stderr, "(First k-mer %d %.*s seen %d times in index)\n", (int) thekmer, (int) kidx->k, thequery->s, (int) kidx->kxr_length[thekmer])
         ;  for (i=0; i<kidx->kxr_length[thekmer]; i++)
            {  int reqid = kidx->kxr[kidx->kxr_offset[thekmer]+i].req_index
            ;  if (!theref[reqid].busy)
               {  todo[n_todo++] = reqid
               ;  theref[reqid].busy++
            ;  }
            }
         }
;if(1)fprintf(stderr, "(Query hits %d index targets)\n", (int) n_todo);
      }

      for (i=0; i<n_todo; i++)
      {  struct req* rs = theref+todo[i]
      ;  SWNUM* thedata = NULL

      ;  if (rs == thequery)
         continue

      ;  rs->busy = 0
      ;  thedata = wrap_sw_fill2(&ai, rs->s, thequery->s, g_indel_allowed, swp)

      ;  sw_trace2a(&ai, swp, ai.max_ij, SW_TRIMLEFT | SW_TRIMRIGHT)          /* M1 */

      ;  if (id_threshold)
         {  if (id_threshold && ai.aln_identity >= id_threshold)
            do_align(fpo, rs->s, rs->annot, thequery->s, thequery->annot, swp, todo[i])
      ;  }
         else if
         (  ai.aln_identity > best_identity
         || (ai.aln_identity == best_identity && ai.aln_overhang < best_overhang)
         )
         {  best_i = todo[i]
         ;  best_identity = ai.aln_identity
         ;  best_overhang = ai.aln_overhang
      ;  }
         if (thedata)
         free(thedata)
   ;  }

      overhang[0] = best_overhang
   ;  return id_threshold ? -1 : best_i
;  }




                  /* requires 0 terminated list */
int find_dna
(  char* buf
,  int buf_n
,  int* start
)
   {  int best_start = 0, best_len = 0, s = 0, len = 0
   ;  char* a = buf, *z = buf + buf_n
   ;  while (a <= z)
      {  int c = ((unsigned char*) a)[0]
      ;  if (c == 'U')
         c = a[0] = 'T'
      ;  if (BASEMAP(c) < 4)
         {  if (!len)
            s = a - buf
         ;  len++
      ;  }
         else if (len > best_len)  
         {  best_start = s
         ;  best_len = len
         ;  len = 0
      ;  }
         a++
   ;  }
      start[0] = best_start
   ;  return best_len
;  }


void printseq
(  struct req* rs
)
   {  fprintf(stderr, ">%s\n%s\n", rs->annot, rs->s)
;  }


unsigned req_total_kmer_length
(  struct req* rs
,  unsigned k
)
   {  unsigned tl = 0
   ;  while (rs->s)
      {  if (rs->slen >= k)
         tl += rs->slen + 1 -k
      ;  rs++
   ;  }
      return tl
;  }


int kxr_cmp_req
(  const void*  x
,  const void* y
)
   {  const struct kmer_x_req* xx = x
   ;  const struct kmer_x_req* yy = y
   ;  if (xx->req_index < yy->req_index)
      return -1
   ;  if (xx->req_index > yy->req_index)
      return 1
   ;  return 0
;  }


int kxr_cmp_kmr
(  const void*  x
,  const void* y
)
   {  const struct kmer_x_req* xx = x
   ;  const struct kmer_x_req* yy = y
   ;  if (xx->kmer < yy->kmer)
      return -1
   ;  if (xx->kmer > yy->kmer)
      return 1
   ;  return 0
;  }


void make_index
(  struct req* rs
,  struct kmer_index* kidx
,  unsigned k
)
   {  unsigned total_length = req_total_kmer_length(rs, k)
#define KMASK ((1 << 2*k) -1)
   ;  struct req* r = rs
   ;  struct kmer_x_req* kxr
   ;  unsigned xi = 0

   ;  kidx->k = k
   ;  if (!k)
      return

   ;  kidx->kxr_offset = myalloc((KMASK+1) * sizeof kidx->kxr_offset[0])
   ;  kidx->kxr_length = myalloc((KMASK+1) * sizeof kidx->kxr_length[0])
   ;  kxr = kidx->kxr  = myalloc(total_length * sizeof kxr[0])

   ;  argh("swan", "building %d-mer index (use -index 0 to run without index)", (int) k)

   ;  while (r->s)
      {  unsigned i, kmer = 0
      ;  int last_N = -1
      ;  for ( i=0; i < r->slen; i++ )
         {  unsigned base = BASEMAP((unsigned char) r->s[i])
         ;  kmer = (kmer << 2)
         ;  if (base == 4)
            last_N = i
         ;  else if (base < 4)
            kmer |= base
         ;  if (last_N + k <= i)
            {  kxr[xi].kmer = kmer & KMASK
            ;  kxr[xi].req_index = r - rs
            ;  xi++
         ;  }
         }
         r++
   ;  }

      {  int i
      ;  for (i=0; i<=KMASK; i++)         /* all counts to zero */
         kidx->kxr_length[i] = 0

      ;  for (i=0; i<xi; i++)             /* count what's there .. */
         kidx->kxr_length[kxr[i].kmer]++

      ;  kidx->kxr_offset[0] = 0          /* set offsets */
      ;  for (i=1; i<=KMASK; i++)
         kidx->kxr_offset[i] = kidx->kxr_offset[i-1] + kidx->kxr_length[i-1]

      ;  qsort(kxr, xi, sizeof kxr[0], kxr_cmp_kmr)
      ;  for (i=0; i<=KMASK; i++)
         qsort(kxr+kidx->kxr_offset[i], kidx->kxr_length[i], sizeof kxr[0], kxr_cmp_req)
   ;  }
   }



   /* sequences can be arbitrarily long
    * Entire file is read into one chunk of memory,
    * then parsed into its constituent sequences.
   */ 

struct req* read_fasta_file
(  ZFILE fp
,  unsigned*  n_ref
,  const char* substr
,  int minlen
)
   {  unsigned long nbytes
   ;  char* fdata = (char*) read_a_file(fp, &nbytes)
   ;  char* a = fdata
   ;  char* start = NULL
   ;  unsigned nseq, i_rec = 0
   ;  struct req* rs = NULL
   ;  struct req rs0 = { 0 }
   ;  int n_skipped = 0

   ;  if (!fdata)
      exit(1)

   ;  nseq = a[0] == '>'

   ;  if (a[0] != '>')
      {  argh("fasta read", "file does not start with >, worryingly")
      ;  a = strstr(a, "\n>")
      ;  if (!a)
         a = fdata+nbytes
      ;  else
            nseq = 1
         ,  a++
   ;  }
      start = a

   ;  while ((a = strstr(a, "\n>")))
         nseq++
      ,  a++

   ;  fprintf(stderr, "read %d sequences\n", (int) nseq)
   ;  if (!(rs = malloc((nseq+1) * sizeof rs[0])))       /* add sentinel value */
      exit(1)

   ;  a = start
   ;  while (a && a[0] == '>')
      {  char* b = strstr(a, "\n"), *c = b, *x, *dest

      ;  if (!b)                       /* should not happen, really */
            argh("fasta read", "error")
         ,  exit(1)

      ;  rs[i_rec] = rs0
      ;  rs[i_rec].annot = a+1
      ;  c = strstr(b, "\n>")          /* these two are very much order-dependent */
      ;  rs[i_rec].annot[b-a-1] = '\0'

      ;  rs[i_rec].s = b+1
      ;  rs[i_rec].slen = c ? c - b - 1 : (nbytes - 1) - ((b+1) - fdata)
      ;  rs[i_rec].s[rs[i_rec].slen] = '\0'
      ;  dest = rs[i_rec].s

      ;  for (x = rs[i_rec].s; x < rs[i_rec].s + rs[i_rec].slen; x++)
         {  if (isalpha(x[0]))
            x[0] = toupper(x[0])
         ;  if (x[0] == 'U')
            x[0] = 'T'

         ;  if (BASEMAP((unsigned char) x[0]) < 4)
            /* NOTHING */
         ;  else if (isspace(x[0]))
            x[0] = '\0'
         ;  else
            x[0] = 'N'

         ;  if (x[0])
            dest++[0] = x[0]
      ;  }

         dest[0] = '\0'
      ;  rs[i_rec].slen = dest - rs[i_rec].s

      ;  a = c ? c+1 : NULL
      ;  if (substr && !strstr(rs[i_rec].annot, substr))
         n_skipped++
      ;  else if (minlen && rs[i_rec].slen < minlen)
         n_skipped++
      ;  else
         i_rec++
   ;  }
      if (i_rec + n_skipped == nseq)
      rs[i_rec] = rs0
   ;  else
      argh("fasta read", "record count discrepancy (%d/%d) - demons", (int) i_rec, (int) nseq)

   ;  n_ref[0] = i_rec
   ;  return rs
;  }


int main
(  int argc
,  char* argv[]
)
   {  FILE* fpo = stdout
   ;  ZFILE fpref = NULL, fpquery = NULL
   ;  const char* g_fnout = "-"
   ;  const char* g_string_query = NULL
   ;  const char* g_string_ref = NULL

   ;  struct sw_param swp = { 0 }
   ;  unsigned n_truncated = 0
   ;  unsigned index_kmer = 0
   ;  unsigned theid = 0
   ;  unsigned limit = 0

   ;  const char* theref = NULL, *thequery = NULL
   ;  const char* fnref = NULL, *fnquery = NULL

   ;  swp.cost_gapleft      =  3
   ;  swp.cost_gapright     =  3
   ;  swp.cost_subst        =  1
   ;  swp.gain_match        =  4
   ;  swp.left_skip         =  0
   ;  swp.right_limit       =  0
   ;  swp.flags             =  0

   ;  kraken_readline(NULL, NULL, 0, NULL, &n_truncated)            /* resets static buffers */
   ;  themap_init()

/* enter macromagical option world */
      
   ;  arg_switch()

      uniarg("-h")
puts("-o                output file name (STDOUT)");
puts("-r FASTA-file     fasta file for reference");
puts("-q FASTA-file     fasta file for query");
puts("-rs DNA-string    reference string to align (displayed on top)");
puts("-qs DNA-string    query string to align (displayed below)");
puts("-q-len <int>      only consider sequences at least this long");
puts("-r-len <int>      only consider sequences at least this long");
puts("-q-string <string> e.g. hsa, mmu; only matching identifiers are considered");
puts("-r-string <string> e.g. hsa, mmu; only matching identifiers are considered");
puts("-id <int>         display matches with at least <int> identity (0-100)");
puts("-index <int>      k-mer size to build index on (suggest 8 to 12; filters on k-mer match!)");
puts("-swp M/S/G        match/substitution/gap : gain/cost/cost");
puts("-lsrl L/R         reference/left-skip / query/right-limit (adapter specific)");
puts("-trace [olr]+     o: use old code  l: trim left r: trim right");
puts("--noindel         do not consider indels while aligning");
puts("--matrix          dump alignment matrix");
puts("--key-value       output easily parseable line-based key-value output");
puts("--excise          excise the aligned part when printing");
puts("-do <int>         process the top <int> entries from the reference file");
puts("-cell <int>       align from cell <int>");
exit(0);
      endarg()

      optarg("-o")  g_fnout = thearg(); endarg()
      optarg("-r-string") g_string_ref = thearg(); endarg()
      optarg("-q-string") g_string_query = thearg(); endarg()
      optarg("-q-len") g_qlen = atoi(thearg()); endarg()
      optarg("-r-len") g_rlen = atoi(thearg()); endarg()
      optarg("-trace") g_trace = thearg(); endarg()
      optarg("-cell") g_cell = atoi(thearg()); endarg()
      optarg("-swp")
         if
         (  3 != sscanf(thearg(), "%u/%u/%u", &swp.gain_match, &swp.cost_subst, &swp.cost_gapleft)
         )
         die(1, "-swp requires %%u/%%u/%%u format");
         swp.cost_gapright = swp.cost_gapleft;
      endarg()
      optarg("-lsrl")
         if
         (  2 != sscanf(thearg(), "%u/%u", &swp.left_skip, &swp.right_limit)
         )
         die(1, "-lsrl requires %%u/%%u format");
      endarg()
      uniarg("--noindel") g_indel_allowed = 0; endarg()
      uniarg("--matrix") g_dump_matrix = 1; endarg()
      uniarg("--key-value") g_mode_print |= PRINT_KEY_VALUE; endarg()
      uniarg("--excise") g_mode_print |= PRINT_EXCISE; endarg()
      optarg("-do")  limit  = atoi(thearg());  endarg()
      optarg("-qs")  thequery  = thearg();  endarg()
      optarg("-id")  theid     = atoi(thearg()); endarg()
      optarg("-rs")  theref   = thearg();  endarg()
      optarg("-index") index_kmer = atoi(thearg()); if (index_kmer > 12) index_kmer = 12; endarg()
      optarg("-q")  fnquery  = thearg();  endarg()
      optarg("-r")  fnref   = thearg();  endarg()
      failarg()
      arg_done()

/* exit macromagicalitaciousness */

   ;  if (!fpquery && thequery && strlen(thequery) < index_kmer)
      {  argh("swan", "query length is very small; index disabled accordingly")
      ;  index_kmer = 0
   ;  }

      if (!(fpo = myfopen(g_fnout, "w", 0)))
      exit(1)
   ;  if (fnref && !(fpref = myfopen(fnref, "r", 1)))
      exit(1)
   ;  if (fnquery && !(fpquery = myfopen(fnquery, "r", 1)))
      exit(1)

   ;  fprintf(stderr, "Identity set to %d, index to %d\n", (int) theid, (int) index_kmer)

   ;  if (fpref)
      {  unsigned nref  = 0, nquery = 0
      ;  int i, overhang = 0
      ;  struct req* ls_ref = read_fasta_file(fpref, &nref, g_string_ref, g_rlen)

      ;  unsigned* todo = myalloc(nref * sizeof todo[0])
      ;  unsigned n_todo = 0

      ;  struct kmer_index ref_index = { 0 }

      ;  if (!index_kmer)
         {  n_todo = nref
         ;  for (i=0; i<n_todo; i++)
            todo[i] = i
      ;  }

         make_index(ls_ref, &ref_index, index_kmer)   /* works always .. but no index when k == 0 */

      ;  if (thequery)
         {  struct req qreq = { 0 }
         ;  qreq.s = stringle("%s", thequery)
         ;  qreq.annot = stringle("(query)")
         ;  if ((i = search_best_match(fpo, &ref_index, &swp, ls_ref, &qreq, &overhang, todo, n_todo, theid)) >= 0)
            do_align(fpo, ls_ref[i].s, ls_ref[i].annot, thequery, "(query)", &swp, i)
      ;  }
         else
         {  struct req* ls_query = fpquery ? read_fasta_file(fpquery, &nquery, g_string_query, g_qlen) : ls_ref
         ;  struct req* q = ls_query
         ;  g_self = ls_ref == ls_query
         ;  while (q->s)
            {  if ((i = search_best_match(fpo, &ref_index, &swp, ls_ref, q, &overhang, todo, n_todo, theid)) >= 0)      /* M0 */
               do_align(fpo, ls_ref[i].s, ls_ref[i].annot, q->s, q->annot, &swp, i)
            ;  q++
            ;  if (limit && q-ls_query >= limit)
               break
         ;  }
         }
      ;  free(todo)
   ;  }
      else if (theref)
      do_align(fpo, theref, "reference", thequery ? thequery : theref, "query", &swp, 1)
   ;  else
      die(1, "no reference!")

   ;  fclose(fpo)
   ;  return 0
;  }


/* dummy change */
