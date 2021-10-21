
/*
 * (C) Copyright 2011, 2012, 2013, 2014, 2015 European Molecular Biology Laboratory.
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


#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "sw.h"
#include "slib.h"

#define elem(ai,i,j)  ai->data[((i) * ai->nj) + (j)]


/* data is layed out as follows,

                x    0  1  2  .. .. ..  nj-1
                  +----^--^--^--------^---
               0  |  0  0  0  0  ..  0  0
               1  |  0  .  .. .. .. .. ..
               2  |  0  :  .           ..
               :  |  0  :     .        ..
               :  |  0  :        .     ..
               :  |  0  :           .  ..
             ni-1 |  0  .  .. .. .. .. ..

      nj == strlen(right) + 1
      ni == strlen(left) + 1

   In memory this is a concatenation of rows.  The first row and column are set
   to all-zero.  sometimes we use elem(ai, i, j), but other times we use
   ai->data[ij] and compute dependent i and j as

      i = ij / ai->nj
      j = ij - (i * ai->nj)
*/



unsigned sw_fit
(  struct sw_alninfo* ai
,  unsigned ij
)
   {  unsigned i = ij / ai->nj
   ;  unsigned j = ij - (i * ai->nj)
   ;  unsigned ret = ij
   ;  if (j && elem(ai, i, j-1) == elem(ai, i, j))
      ret = ij - 1
   ;  if (i && elem(ai, i-1, j) == elem(ai, i, j))
      ret = ij - ai->nj
;fprintf(stderr, "ij %d ret %d\n", (int) ij, (int) ret)
   ;  return ret
;  }


int sw_fill2
(  struct sw_alninfo* ai
,  SWNUM* data
,  unsigned data_size
,  const char *left           /* fixme, length-encode */
,  const char *right          /* fixme, length-encode */
,  int indel_allowed
,  struct sw_param* swp
)
   {  int i, j
   ;  const int COST_GAPLEFT = swp->cost_gapleft
   ;  const int COST_GAPRIGHT = swp->cost_gapright
   ;  const int COST_SUBST = swp->cost_subst
   ;  const int GAIN_MATCH = swp->gain_match
   ;  const unsigned left_skip = swp->left_skip
   ;  const unsigned right_limit = swp->right_limit
   ;  unsigned flags = swp->flags
   ;  memset(ai, 0, sizeof ai[0])
   ;  ai->data = data
   ;  ai->left = left
   ;  ai->right = right
   ;  ai->ni   = strlen(left) + 1
   ;  ai->nj   = strlen(right) + 1
   ;  unsigned mm = 0, ij = 0
   ;  unsigned e = ai->ni-1, mme = 0, ej = 0     /* at end of left string */
   ;  unsigned lowerleft = flags & SW_ENDTOSTART

   ;  if (indel_allowed)
      swp->flags |= SW_INDELALLOWED

   ;  if ( ai->ni * ai->nj > data_size )
      return 1

   ;  memset(data, 0, ai->nj * ai->ni * sizeof data[0])

   ;  ai->n_cells_used = 0

#define MAX_ABS_DELTA 4
#if 0
#define GAIN_MATCH  2
#define COST_GAP    1
#define COST_SUBST  1
#endif

   ;  for (i=left_skip+1;i<ai->ni;i++)
      {  unsigned maxj = ai->nj
      ;  if (lowerleft && i+1 < maxj)
         maxj = i + 1
      ;  if (right_limit && right_limit+1 < maxj)
         maxj = right_limit + 1
         
      ;  ai->n_cells_used += maxj

      ;  elem(ai, i, 0) = 0

      ;  for (j=1;j<maxj;j++)
         {  unsigned char l = left[i-1]
         ;  unsigned char r = right[j-1]
         ;  unsigned match  = l != r || l == 'N' || r == 'N' ?  0 : 1
         ;  int newval = GAIN_MATCH * match            /* assuming match starts here */
         ;  int iidelta = match ? GAIN_MATCH : -COST_SUBST
         ;  elem(ai,i,j) = 0

         ;  if (indel_allowed)
            {  if (elem(ai, i, j-1) > newval + COST_GAPLEFT)   /* step right; gap in left */
               newval = elem(ai, i, j-1) - COST_GAPLEFT
            ;  if (elem(ai, i-1, j) > newval + COST_GAPRIGHT)  /* step left; gap in right */
               newval = elem(ai, i-1, j) - COST_GAPRIGHT
         ;  }
            if (elem(ai, i-1, j-1) + (MAX_ABS_DELTA + iidelta) > (MAX_ABS_DELTA + newval))
            newval = elem(ai, i-1, j-1) + iidelta
         ;  if (newval >= 0)
            elem(ai, i, j) = newval
         ;  if (elem(ai,i,j) > mm)
               ij = i * ai->nj + j
            ,  mm = elem(ai,i,j)
      ;  }
      }
      for (j=0;j<ai->nj;j++)
      {  if (elem(ai, e, j) > mme)
            ej = e * ai->nj + j
         ,  mme = elem(ai, e, j)
   ;  }
      ai->max_ij = ij   /* sw_fit(ai, ij) */
   ;  ai->max_ej = ej   /* sw_fit(ai, ej) */
   ;  return 0
;  }


void sw_dump
(  struct sw_alninfo* ai
)
   {  int i, j

   ;  for (i=0;i<ai->ni;i++)
      {  fprintf(stdout, "%3d %c: ", i, i ? (int) ai->left[i-1] : '-')
      ;  for (j=0;j<ai->nj;j++)
         fprintf(stdout, "%3d", (int) elem(ai,i,j))
      ;  fprintf(stdout, "  : %d\n", (int) ((i+1) * ai->nj - 1))
   ;  }

   ;  fputs("\n  - -    -", stdout)
   ;  for (i=1;i<ai->nj;i++)
      fprintf(stdout, "%3c", (int) ai->right[i-1])
   ;  fputc('\n', stdout)

   ;  fputs("  - -    -", stdout)
   ;  for (i=1;i<ai->nj;i++)
      fprintf(stdout, "%3d", i)
   ;  fprintf(stdout, "\n   max cell %d value %d\n", ji ai->max_ij, ji ai->data[ai->max_ij])
   ;  fprintf(stdout,   "   max last-row cell %d value %d\n", ji ai->max_ej, ji ai->data[ai->max_ej])
;  }


void sw_trace2a
(  struct sw_alninfo* ai
,  struct sw_param* swp
,  unsigned ij
,  unsigned modes
)
   {  unsigned i = ij / ai->nj
   ;  unsigned j = ij - (i * ai->nj)
   ;  unsigned nj = ai->nj
   ;  SWNUM* data = ai->data
   ;  unsigned n_subst = 0
   ;  unsigned n_match = 0
   ;  unsigned n_insl = 0
   ;  unsigned n_insr = 0
   ;  int indel_allowed = swp->flags & SW_INDELALLOWED

   ;  ai->lft_end = i
   ;  ai->rgt_end = j

   ;  int o_start = SW_ALN_MAXSIZE

   ;  if (ij)
      do
      {  unsigned new11 = ij-nj-1            /* diagonal step          */
      ;  unsigned new10 = ij-1               /* step left; gap in left */
      ;  unsigned new01 = ij-nj              /* step right: gap in right */
      ;  unsigned newij = new11
      ;  unsigned char theleft  = ai->left[i-1]
      ;  unsigned char theright = ai->right[j-1]
      ;  unsigned char status = theleft == 'N' || theright == 'N' ? 'x' : '|'

      ;  if (data[ij] <= 0 || --o_start < 0)
         break

      ;  ai->aln_lft[o_start] = ai->left[i-1]
      ;  ai->aln_rgt[o_start] = ai->right[j-1]

      ;  if (theleft == theright)
         newij = new11
                                             /* prefer match always */
      ;  else if (!indel_allowed || data[new11] > data[ij])     /* substitution */
            status = 'x'
         ,  newij = new11
      ;  else if (indel_allowed && (data[new10] > data[ij]))    /* gap in left  */
            status = 'l'
         ,  ai->aln_lft[o_start] = '-'
         ,  newij = new10
      ;  else if (indel_allowed && data[new01] > data[ij])      /* gap in right  */
            status = 'r'
         ,  ai->aln_rgt[o_start] = '-'
         ,  newij = new01

      ;  ai->aln_aln[o_start] = status

      ;  ij = newij
      ;  i = ij / ai->nj
      ;  j = ij - (i * ai->nj)
;if(0)fprintf(stderr, "\n} %s {\n} %s {\n} %s [%c]{\n", ai->aln_rgt+o_start, ai->aln_aln+o_start, ai->aln_rgt+o_start, (int) ai->aln_rgt[SW_ALN_MAXSIZE])
;if (swp->flags & SW_VERBOSE)
fprintf
(   stderr
,   "i %d j %d status %c val %d\n"
,   (int) i
,   (int) j
,   (int) status
,   (int) data[newij]
)
   ;  }
      while (ij && data[ij] > 0)

   ;  {  int o_right = SW_ALN_MAXSIZE              /* EXCLUSIVE */
      ;  int o_left  = o_start                     /* inclusive */
      ;  int n_left_trim_lr = 0
      ;  int n_left_trim_r = 0
      ;  int n_left_trim_l = 0
      ;  int o

      ;  if (modes & SW_TRIMRIGHT)
         {  int trimscore = 0
         ;  for (o = SW_ALN_MAXSIZE-1; o > o_start; o--)
            {  trimscore += ai->aln_aln[o] == '|' ? 1 : -1
            ;  if (trimscore <= 0 && ai->aln_aln[o] != '|')
               o_right = o
         ;  }
         }

         if (modes & SW_TRIMLEFT)
         {  int trimscore = 0
         ;  for (o = o_start; o < o_right; o++)
            {  trimscore += ai->aln_aln[o] == '|' ? 1 : -1
            ;  if (trimscore <= 0 && ai->aln_aln[o] != '|')
               o_left = o+1
         ;  }
         }

         ai->aln_lft[o_right] = '\0'
      ;  ai->aln_aln[o_right] = '\0'
      ;  ai->aln_rgt[o_right] = '\0'

      ;  for (o=o_start;o<o_left;o++)
         {  switch((unsigned char) ai->aln_aln[o])
            {  case '|' : n_left_trim_lr++; break
            ;  case 'x' : n_left_trim_lr++; break
            ;  case 'l' : n_left_trim_r++; break
            ;  case 'r' : n_left_trim_l++; break
         ;  }
         }

      ;  for (o=o_left;o<o_right;o++)
         {  switch((unsigned char) ai->aln_aln[o])
            {  case '|' : n_match++; break
            ;  case 'x' : n_subst++; break
            ;  case 'l' : n_insl++; ai->aln_aln[o] = '-'; break
            ;  case 'r' : n_insr++; ai->aln_aln[o] = '-'; break
         ;  }
         }

         ai->aln_ofs = o_left
      ;  ai->aln_end = o_right
      ;  ai->n_match = n_match
      ;  ai->n_subst = n_subst
      ;  ai->n_insl  = n_insl
      ;  ai->n_insr  = n_insr

      ;  ai->lft_start =  ( 1 + ij / ai->nj )  + n_left_trim_lr + n_left_trim_l
      ;  ai->rgt_start =  ( 1 + ij - ((ij / ai->nj) * ai->nj) )  + n_left_trim_lr + n_left_trim_r
      ;  ai->lft_end   =  ai->lft_start + ai->n_match + ai->n_subst + ai->n_insr -1
      ;  ai->rgt_end   =  ai->rgt_start + ai->n_match + ai->n_subst + ai->n_insl -1

#define SW_BASE_ID_RIGHT(ai) ((ai)->n_match * 1.0 / ((ai)->nj-1+(ai)->n_insr))
#define SW_BASE_ID_LEFT(ai)  ((ai)->n_match * 1.0 / ((ai)->ni-1+(ai)->n_insl))
      ;  {  double id_right = SW_BASE_ID_RIGHT(ai)
         ;  double id_left  = SW_BASE_ID_LEFT(ai)
         ;  double id_best  = id_right > id_left ? id_right : id_left
         ;  ai->aln_identity = (int) (100.5 * id_best)
         ;  ai->aln_overhang = SW_OVERHANG(ai)
      ;  }
      }

;if(0)fprintf(stderr, "\n! %s {\n} %s {\n} %s [%c]{\n", ai->aln_rgt+o_start, ai->aln_aln+o_start, ai->aln_rgt+o_start, (int) ai->aln_rgt[SW_ALN_MAXSIZE])

;if(0)fprintf(stderr, "match %d gl %d gr %d sub %d ofs %d [%s] [%s]\n"
,(int) ai->n_match
,(int) ai->n_insl
,(int) ai->n_insr
,(int) ai->n_subst
,(int) ai->aln_ofs
, ai->left, ai->right
)
;  }


void sw_trace2
(  struct sw_alninfo* ai
,  struct sw_param* swp
,  int indel_allowed
,  unsigned ij
)
   {  unsigned i = ij / ai->nj
   ;  unsigned j = ij - (i * ai->nj)
   ;  unsigned nj = ai->nj
   ;  SWNUM* data = ai->data
   ;  unsigned n_subst = 0
   ;  unsigned n_match = 0
   ;  unsigned n_insl = 0
   ;  unsigned n_insr = 0

   ;  ai->lft_end = i
   ;  ai->rgt_end = j

   ;  int o = SW_ALN_MAXSIZE

   ;  ai->aln_lft[o] = '\0'
   ;  ai->aln_aln[o] = '\0'
   ;  ai->aln_rgt[o] = '\0'

   ;  if (ij)
      do
      {  unsigned new11 = ij-nj-1            /* diagonal step          */
      ;  unsigned new10 = ij-1               /* step left; gap in left */
      ;  unsigned new01 = ij-nj              /* step right: gap in right */
      ;  unsigned newij = new11
      ;  unsigned char status = 'm'

      ;  if (--o < 0)                        /* fixme: design such that this is unnecessary */
         break

      ;  ai->aln_lft[o] = ai->left[i-1]
      ;  ai->aln_rgt[o] = ai->right[j-1]

      ;  if (data[new11] > data[ij])         /* substitution */
            status = 's'
         ,  newij = new11
      ;  else if (indel_allowed && (data[new10] > data[ij]))    /* gap in left  */
            status = 'l'
         ,  ai->aln_lft[o] = '-'
         ,  newij = new10
      ;  else if (indel_allowed && data[new01] > data[ij])      /* gap in right  */
            status = 'r'
         ,  ai->aln_rgt[o] = '-'
         ,  newij = new01

      ;  ai->aln_aln[o] = status == 'm' ? '|' : status == 's' ? 'x' : '-'

      ;  switch(status)
         {  case 'm' : n_match++; break
         ;  case 's' : n_subst++; break
         ;  case 'l' : n_insl++; break
         ;  case 'r' : n_insr++; break
      ;  }

         ij = newij
      ;  i = ij / ai->nj
      ;  j = ij - (i * ai->nj)
;if(0)fprintf(stderr, "\n} %s {\n} %s {\n} %s [%c]{\n", ai->aln_rgt+o, ai->aln_aln+o, ai->aln_rgt+o, (int) ai->aln_rgt[SW_ALN_MAXSIZE])
;if (swp->flags & SW_VERBOSE)
fprintf
(   stderr
,   "i %d j %d status %c val %d\n"
,   (int) i
,   (int) j
,   (int) status
,   (int) data[newij]
)
   ;  }
      while (ij && data[ij] > 0)
;if(0)fprintf(stderr, "\n! %s {\n} %s {\n} %s [%c]{\n", ai->aln_rgt+o, ai->aln_aln+o, ai->aln_rgt+o, (int) ai->aln_rgt[SW_ALN_MAXSIZE])
   ;  ai->aln_ofs = o
   ;  ai->aln_end = SW_ALN_MAXSIZE
   ;  ai->n_match = n_match
   ;  ai->n_subst = n_subst
   ;  ai->n_insl  = n_insl
   ;  ai->n_insr  = n_insr
   ;  ai->lft_start =  1 + ij / ai->nj
   ;  ai->rgt_start =  1 + ij - (ij / ai->nj * ai->nj)
;  }


void sw_printaln
(  struct sw_alninfo* ai
,  FILE* out
)
   {  unsigned o = ai->aln_ofs
   ;  fprintf(out, "%s :\n%s :\n%s :\n", ai->aln_lft+o, ai->aln_aln+o, ai->aln_rgt+o)
;  }


#if 0
void sw_pp
(  struct sw_alninfo* ai
,  int accept
,  int score
,  FILE* fp
,  int zippit
,  int recno
)
   {  unsigned o = ai->aln_ofs
   ;  char *space = "                                                                                "
   ;  char buf[8192]
   ;  int n
   ;  int leftflush = ai->lft_start > ai->rgt_start
   ;  int shift = (int) ai->lft_start - (int) ai->rgt_start

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
      ,  ai->left + (ai->lft_end) /* mq */
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
      ,  ai->right + (ai->rgt_end) /* mq */
/* 5 */
      ,  (int) ai->rgt_start
      ,  (int) ai->rgt_end
/* 6 */
      ,  (int) accept
      ,  (int) score
/* 7 */
      ,  (int) recno
      )
   ;  if (n > 0)
      fputs(buf, fp)

   ;  fprintf(fp, "%s\n%s\n%s\n", ai->aln_lft+o, ai->aln_aln+o, ai->aln_rgt+o)
;  }
#endif


void sw_trace3
(  struct sw_alninfo* ai
,  struct sw_param* swp
,  unsigned ij
,  unsigned modes
)
   {  return
         0
      ?  sw_trace2(ai, swp, swp->flags & SW_INDELALLOWED , ij)
      :  sw_trace2a(ai, swp, ij, modes)
;  }


#define THE_BUFSIZE 8192

void sw_pp2
(  const char* annot_left
,  const char* annot_right
,  int cell
,  struct sw_alninfo* ai
,  void* fp
,  int zippit
,  int recno
,  int(*cb)(const struct sw_alninfo* ai, char* buf, unsigned bufsize)
)
   {  unsigned o = ai->aln_ofs
   ;  char buf[THE_BUFSIZE]
   ;  char space[512]
   ;  int n
   ;  int leftflush = ai->lft_start > ai->rgt_start
   ;  int shift     = (int) ai->lft_start - (int) ai->rgt_start
   ;  int cell_end  = ai->ni * ai->nj
   ;  int cell_used = cell > ai->nj && cell < cell_end ? cell : cell_end-1
   ;  int score = ai->data[cell_used]
;if(0)fprintf(stderr, "used %d last %d\n", (int) cell_used, cell_end)
   ;  memset(space, ' ', 512)
   ;  space[511] = '\0'

   ;  if (score < 0)
      score = ai->data[ai->max_ij]

   ;  if (shift < 0)
      shift = -shift

   ;  n
   =  snprintf
      (  buf
      ,  THE_BUFSIZE
      ,  " %.*s%.*s%s%s : annot=%s\n %.*s%.*s%s : score=%d id=%d"
       /*   1 1 1 1 1 1          2    2 2 2 2 2          3     3  */
      ,  leftflush ? 0 : shift         /* this much */
      ,  space
      ,  (int) (ai->lft_start-1)       /* this much */
      ,  ai->left
      ,  ai->aln_lft+o
      ,  ai->left + (ai->lft_end)
/* 2 */
      ,  annot_left
      ,  shift
      ,  space
      ,  (int) ((leftflush ? ai->rgt_start : ai->lft_start) - 1)
      ,  space
      ,  ai->aln_aln+o
/* 3 */
      ,  score
      ,  ai->aln_identity
      )

   ;  if (n > 0)
      {
#if WE_USE_ZLIB
         if (zippit) gzwrite(fp, buf, n); else fputs(buf, fp)
#else
         fputs(buf, fp)
#endif
   ;  }

      if (cb)
      {  n = cb(ai, buf, THE_BUFSIZE);
#if WE_USE_ZLIB
         if (zippit) gzwrite(fp, buf, n); else fputs(buf, fp)
#else
         fputs(buf, fp)
#endif
   ;  }

      n
   =  snprintf
      (  buf
      ,  THE_BUFSIZE
      ,  "\n %.*s%.*s%s%s : %s%s\n---- recno=%d\n\n"
       /*     4 4 4 4 4 4    5 5              5 */
      ,  leftflush ? shift : 0
      ,  space
      ,  (int) (ai->rgt_start -1)
      ,  ai->right
      ,  ai->aln_rgt+o
      ,  ai->right + (ai->rgt_end)
/* 5 */
      ,  annot_right ? " annot=" : ""
      ,  annot_right ? annot_right : ""
      ,  (int) recno
      )

   ;  if (n > 0)
      {
#if WE_USE_ZLIB
         if (zippit) gzwrite(fp, buf, n); else fputs(buf, fp)
#else
         fputs(buf, fp)
#endif
   ;  }

   }


void sw_pp_excise
(  const char* annot_left
,  const char* annot_right
,  int cell
,  struct sw_alninfo* ai
,  void* fp
,  int zippit
,  int recno
,  int(*cb)(const struct sw_alninfo* ai, char* buf, unsigned bufsize)
)
   {  unsigned o = ai->aln_ofs
   ;  char buf[THE_BUFSIZE]
   ;  char space[512]
   ;  int n
   ;  int shift     = (int) ai->lft_start - (int) ai->rgt_start
   ;  int cell_end  = ai->ni * ai->nj
   ;  int cell_used = cell > ai->nj && cell < cell_end ? cell : cell_end-1
   ;  int score = ai->data[cell_used]
;if(0)fprintf(stderr, "used %d last %d\n", (int) cell_used, cell_end)
   ;  memset(space, ' ', 512)
   ;  space[511] = '\0'

   ;  if (score < 0)
      score = ai->data[ai->max_ij]

   ;  if (shift < 0)
      shift = -shift

   ;  n
   =  snprintf
      (  buf
      ,  THE_BUFSIZE
      ,  " %s : annot=%s\n %s : score=%d id=%d"
      ,  ai->aln_lft+o
      ,  annot_left
      ,  ai->aln_aln+o
      ,  score
      ,  ai->aln_identity
      )

   ;  if (n > 0)
      {
#if WE_USE_ZLIB
         if (zippit) gzwrite(fp, buf, n); else fputs(buf, fp)
#else
         fputs(buf, fp)
#endif
   ;  }

      if (cb)
      {  n = cb(ai, buf, THE_BUFSIZE);
#if WE_USE_ZLIB
         if (zippit) gzwrite(fp, buf, n); else fputs(buf, fp)
#else
         fputs(buf, fp)
#endif
   ;  }

      n
   =  snprintf
      (  buf
      ,  THE_BUFSIZE
      ,  "\n %s : %s%s\n---- recno=%d\n\n"
      ,  ai->aln_rgt+o
      ,  annot_right ? " annot=" : ""
      ,  annot_right ? annot_right : ""
      ,  (int) recno
      )

   ;  if (n > 0)
      {
#if WE_USE_ZLIB
         if (zippit) gzwrite(fp, buf, n); else fputs(buf, fp)
#else
         fputs(buf, fp)
#endif
   ;  }

   }


void sw_lp                       /* line print */
(  const char* annot_left
,  const char* annot_right
,  int cell
,  struct sw_alninfo* ai
,  void* fp
,  int zippit
,  int recno
,  int(*cb)(const struct sw_alninfo* ai, char* buf, unsigned bufsize)
)
   {  char buf[THE_BUFSIZE]
   ;  int o = ai->aln_ofs
   ;  int n

   ;  n
   =  snprintf
      (  buf
      ,  THE_BUFSIZE
      ,  "aln-aln={%s} aln-r={%s} aln-q={%s} r={%s} q={%s} r-start=%d q-start=%d aln-id=%d r-annot={%s} q-annot={%s}\n"
      ,  ai->aln_aln+o
      ,  ai->aln_lft+o
      ,  ai->aln_rgt+o
      ,  ai->left
      ,  ai->right
      ,  (int) ai->lft_start
      ,  (int) ai->rgt_start
      ,  (int) ai->aln_identity
      ,  annot_left
      ,  annot_right
      )

   ;  if (n > 0)
      {
#if WE_USE_ZLIB
         if (zippit) gzwrite(fp, buf, n); else fputs(buf, fp)
#else
         fputs(buf, fp)
#endif
   ;  }

   }


void sw_pp_sparse
(  struct sw_alninfo* ai
,  void* fp
,  int zippit
)
   {  unsigned o = ai->aln_ofs
   ;  char buf[THE_BUFSIZE]
   ;  char space[512]
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
      ,  THE_BUFSIZE
      ,  " %.*s%.*s%s%s\n %.*s%.*s%s\n %.*s%.*s%s%s\n"

      ,  leftflush ? 0 : shift         /* this much */
      ,  space
      ,  (int) (ai->lft_start-1)       /* this much */
      ,  ai->left
      ,  ai->aln_lft+o
      ,  ai->left + (ai->lft_end)

      ,  shift                         /* this much */
      ,  space
      ,  (int) ((leftflush ? ai->rgt_start : ai->lft_start) - 1)     /* this much */
      ,  space
      ,  ai->aln_aln+o

      ,  leftflush ? shift : 0         /* this much */
      ,  space
      ,  (int) (ai->rgt_start -1)
      ,  ai->right
      ,  ai->aln_rgt+o
      ,  ai->right + (ai->rgt_end)
      )
   ;  if (n > 0)
      {
#if WE_USE_ZLIB
         if (zippit) gzwrite(fp, buf, n)
      ;  else        fputs(buf, fp)
#else
         fputs(buf, fp)
#endif
   ;  }
   }



   /* fixme: n_stretch in callers often computed using SW_ALN_MAXSIZE, but trimming changes this */
   /* returns number of matches if parameters satisfied, 0 otherwise */
int sw_alignment_edit_ok
(  const unsigned char* aln
,  unsigned aln_length
,  unsigned n_stretch
,  unsigned n_edit_max
,  unsigned n_gap_max
)
   {  int i, n_edit = 0, n_gap = 0, n_subst = 0, n_match = 0, ok = 0
   ;  if (!n_stretch)
      return 0
   ;  unsigned hist[SW_ALN_MAXSIZE] = { 0 }
;if(0)fprintf(stderr, "length %d n_stretch %d e %d\n", (int) aln_length ,  (int) n_stretch ,  (int) n_edit_max)
   ;  for (i=0;i<aln_length;i++)
      {  unsigned char a = aln[i]
      ;  unsigned char e = hist[i % n_stretch]
      ;  hist[i % n_stretch] = a
      ;  n_subst += (a == 'x') - (e == 'x')
      ;  n_gap   += (a == '-') - (e == '-')
      ,  n_edit   = n_subst + n_gap
      ;  n_match += a == '|'
      ;  if (  n_edit <= n_edit_max
            && i+1 + n_edit_max >= n_stretch + n_edit
            && n_gap <= n_gap_max
            )
         ok = 1
   ;  }
      return ok ? n_match : 0
;  }




