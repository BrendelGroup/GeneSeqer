/* sahmtd.c
   Sequence Alignment Hidden Markov Tool - DNA version

   Last update:  February 11, 2019. (VB)
*/


/* Corresponding author:                                                      */

/*   Volker Brendel, Department of Biology                                    */
/*   Indiana University, Bloomington, IN 47405                                */
/*   (812) 855-7074, vbrendel@indiana.edu                                     */

/* Past contributing authors:                                                 */
/*   Jonathan Usuka, Department of Chemistry, Stanford University             */

/*******************************************************************************

    Copyright (C) 2012-2019 Volker Brendel.

    This file is part of GeneSeqer.

    GeneSeqer is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    GeneSeqer is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with GeneSeqer.  If not, see <http://www.gnu.org/licenses/>.

*******************************************************************************/



#define GENESEQER
#define LTFLAG	0


/*      INCLUSIONS      */
#define MATHFUNC					/* Gives access to real math functions */
#define MEMOFUNC					/* Gives access to memory allocation */
#define STRGFUNC					/* Gives access to string manipulation */
#include "platform.h"
#include <stdio.h>
#include "html.h"
#include "sahmt.h"
#include "minmax.h"
#include "strapp.h"
#include "RawTextFile.h"	/* Header, raw text file services */

extern char NAUC[];
extern char AAUC[];
extern int codtoaa[4][4][4];
extern int frompA, top, htmlop, ifwdth;
extern int NDEALLOC_GCA;
int NALLOC_PGL = 0, NDEALLOC_PGL = 0;
int NALLOC_ALN = 0, NDEALLOC_ALN = 0;

#ifdef MPITIMING
extern int MPIV_rank, MPIV_nprc;
#endif


int getlns(int Selection, char detsz, char noffset, char *sfname, char *seq);
void complement_seq(char *seq, int numbp, char *seqR);
void prt_seq_segment_to_str(char *str, char *seq, int numbp, int ia, int ib, int rflag, char ABC[], int nflag);


/*      ^^^^^^^^^ FUNCTION DEFINITIONS ^^^^^^^^^ */

int sahmtD(FILE *outfp, char *gdna, char *gdnaR, int sia, int sib,
  int numbp, int rflag, float *pd, float *pa, float *pdR,
  float *paR, char *cdna, struct gcalgnmnt *gca, struct sprmtr sprm, int calln);

static void CompletePathMatrix(int glgth, int clgth, char *gdna, char *cdna, float *pd, float *pa, float pdg);
static void E_n0(int n);
static void I_n0(int n);
static void E_0m(int m);
static void I_0m(int m);
static void E_nm(int n, int m, char *gdna, char *cdna, int glgth, int clgth, float *pd, float *pa, float pdg);
static void I_nm(int n, int m, char *gdna, char *cdna, int glgth, int clgth, float *pd, float *pa, float pdg);


static void OptimalPath(int glgth, int clgth, char *gdna, char *cdna);
static void TracePath(int statev, char *gdna, char *cdna);
static void WritePath(int g, int c, int statev);

static void OutputEnds(int ia);
static int OutputParse(float *pd, float *pa, int ia, int numbp, int rflag, char *gname, int glgth, char *cdna, char *cname, int clgth, int gextnsn, struct gcalgnmnt *gca);
static void OutputAlignment(int numbp, int rflag, struct gcalgnmnt *gca);

static void InitMatch(float ids, float mms, float nns, float dls);
static float GetScore(int g, int c);
static int IndexMax(float *x, int n);
static double weight(int state, int path, float *pd, float *pa,
  int n, int glgth, int m, int clgth, int g, int c, float pdg);
static void OutputBase(int base);

static void allocate_Dspace(int glgth, int clgth);
static void free_Dspace(int glgth, int clgth);

void read_sahmt_prm(FILE *fp, struct sprmtr *sprm);


int insert_gca(struct gcalgnmnt **gcahpp, struct gcalgnmnt *nodep);
int compare_gcas(struct gcalgnmnt *gca1p, struct gcalgnmnt *gca2p);

int rev_sgl(FILE *fp, struct gcalgnmnt *gca, char *gdna, char *gdnaR, int sia, int sib, int numbp, float *pd, float *pa, float *pdR, float *paR, struct sprmtr sprm);
int is_leftend_pgl(struct gcalgnmnt *gca, int rightend, int delta, int fwdni, int revni, float fwdpdpa, float revpdpa);
int maxorflgth(char *seq, int slgth);
int orient_gca(FILE *fp, int npgl, struct gcalgnmnt *lgca, struct gcalgnmnt *rgca, int rflag, char *gdna, char *gdnaR, int sia, int sib, int numbp, float *pd, float *pa, float *pdR, float *paR, struct sprmtr sprm);
int reverse_gca(FILE *fp, struct gcalgnmnt *gca, char *gdna, char *gdnaR, int sia, int sib, int numbp, float *pd, float *pa, float *pdR, float *paR, struct sprmtr sprm);
void gcacpy(struct gcalgnmnt *gca1,struct gcalgnmnt *gca2);


int prt_gca_list(FILE *fp, struct gcalgnmnt *gcap, char *gsgmntn);
void free_gca(struct gcalgnmnt *gca);


void det_pgl(FILE *fp, struct gcalgnmnt **gcaheadpp, char *gname, char *gsgmntn, char *gdna, char *gdnaR, int sia, int sib, int numbp, float *pd, float *pa, float *pdR, float *paR, struct sprmtr sprm, int bflag, int revestallowed);
int build_pgl_list(FILE *fp, struct pgl **pglhpp, struct gcalgnmnt *gcap, char *gdna, char *gdnaR, int sia, int sib, int numbp, float *pd, float *pa, float *pdR, float *paR, struct sprmtr sprm);
int quality_adj_gca(struct gcalgnmnt *gcap, struct sprmtr sprm);
int weak_exon(struct gcalgnmnt *gcap, int i, int o, struct sprmtr sprm);
int consolidate_gca(struct gcalgnmnt **gcahpp, int *ngca);
int overlap_pgl_gca(FILE *fp,struct pgl *pglp, struct gcalgnmnt *gcap, char *gdna, char *gdnaR, int sia, int sib, int numbp, float *pd, float *pa, float *pdR, float *paR, struct sprmtr sprm, int delta, int *gapl);
int is_xwab(int x,int a,int b,int delta);
int xy_contains_ab(int x,int y,int a,int b,int delta);
int gca2pgl(struct pgl **pglhpp, struct pgl *pglp, struct gcalgnmnt *gcap);
int insert_pgl(struct pgl **pglhpp, struct pgl *pglp);
int overlap_pgl_pgl(struct pgl *pgl1, struct pgl *pgl2, int delta);
int merge_pgl_pgl(struct pgl *pgl1, struct pgl *pgl2);

int form_ags_in_pgl(FILE *fp, struct pgl *pglhp);
int sort_gca_in_pgl(struct gcalgnmnt **gcahpp);
int insert_gca_by_link(struct gcalgnmnt **gcahpp, struct gcalgnmnt *nodep);
int order_gcas(struct gcalgnmnt *gca1p, struct gcalgnmnt *gca2p);
int merge_pgl_gca(struct pgl *pgl, int a, struct gcalgnmnt *gca);
int check_pgl_gca(struct pgl *pgl, int a, int p, struct gcalgnmnt *gca, int g);
int is_consistent(struct pgl *pgl, int a, int p, struct gcalgnmnt *gca, int g);
int adj_ags(struct pgl *pgl, int a, int p, struct gcalgnmnt *gca, int g);
int new_ags(struct pgl *pgl, int a, struct gcalgnmnt *gca);

int consolidate_pgl(struct pgl **pglhpp, int *npgl);

int prt_consensus_pgl(FILE *fp, struct pgl *pgl, char *gname, char *gsgmntn, char *gdna, char *gdnaR, int numbp, int npgl, int bflag);
void free_pgl(struct pgl *pgl);

void OutputORF(FILE *fp, char *gdna, char *gdnaR, char *gname, char *gsgmntn, int numbp, int rflag, struct pgl *pgl, int pgln, int a);
void reverse_pgl(struct pgl *pgl, int a, int numbp);
void find_orfs(FILE *fp, char *seq, char *sname, char *sgmntn, int slgth, int minorfl, int numbp, int rflag, struct pgl *pgl, int pgln, int a);
int insert_orf(struct orf **orfhpp, struct orf *nodep);
int consolidate_orfs(struct orf **orfhpp);
void trl_orf(FILE *fp, char *seq, char *sname, char *sgmntn, int slgth, int firstbp, int lastbp, int frame, int numbp, int rflag, int orfn, struct pgl *pgl, int pgln, int a);
void free_orf(struct orf *orf);


int crdfct(int n, int nsgmts, int iab[MAXNEXNS][2], int numbp, int rflag);
extern void fatal_error(char *buf);




/*      ^^^^^^^^^ FUNCTION IMPLEMENTATIONS ^^^^^^^^^ */

static float Log_Pdg,Log_1_Pdg;
static float *Log_Pd,*Log_1_Pd,*Log_Pa,*Log_1_Pa;



int sahmtD(FILE *outfp, char *gdna, char *gdnaR, int sia, int sib,
  int numbp, int rflag, float *pd, float *pa, float *pdR,
  float *paR, char *cdna, struct gcalgnmnt *gca, struct sprmtr sprm, int calln)
{
  int i, ia, ib, glgth, clgth, opflag, giaorig, giborig;
  char gname[257], cname[257];
  char *lgdna;
  float *lpd, *lpa;
  extern float *logV[8];

#ifdef MPITIMING
double MPIV_start, MPIV_finish, MPI_Wtime();
MPIV_start = MPI_Wtime();
#endif
  /* Precalculate log value of pdg: */
  Log_Pdg=log(sprm.pdg);
  Log_1_Pdg=log(1.-sprm.pdg);

ADJUST:
  strcpy(gname, gca->gname);
  ia = gca->gia;
  ib = gca->gib;
  glgth = ib - ia + 1;
  strcpy(cname, gca->cname);
  clgth = gca->clgth;

  if (rflag == 0) {
    lgdna = &(gdna[gca->gia]);
    lpd = &(pd[gca->gia - sia]);
    lpa = &(pa[gca->gia - sia]);
    Log_Pd=logV[0]+gca->gia-sia;
    Log_1_Pd=logV[1]+gca->gia-sia;
    Log_Pa=logV[2]+gca->gia-sia;
    Log_1_Pa=logV[3]+gca->gia-sia;
  }
  else {
    lgdna = &(gdnaR[gca->gia]);
    lpd = &(pdR[gca->gia - sia]);
    lpa = &(paR[gca->gia - sia]);
    Log_Pd=logV[4]+gca->gia-sia;
    Log_1_Pd=logV[5]+gca->gia-sia;
    Log_Pa=logV[6]+gca->gia-sia;
    Log_1_Pa=logV[7]+gca->gia-sia;
  }

  if (glgth > MAXGLGTH) {
    fprintf(outfp, "\nWARNING: %s gDNA length (%d) exceeds allowed maximum (%d).\n", gname, glgth, MAXGLGTH);
    return (0);
  }
  if (clgth > MAXCLGTH) {
    fprintf(outfp, "\nWARNING: %s cDNA length (%d) exceeds allowed maximum (%d).\n", cname, clgth, MAXCLGTH);
    return (0);
  }

  nextFree=0;

  AppendChar('\n');
  for (i = 0; i < 80; ++i)
    AppendChar('*');
  if (cname[strlen(cname)-1] == '+') {
    if (htmlop) {
      CHECKALIGNLEN; nextFree+=sprintf(global_algnmnt+nextFree, "\n<A HREF=\"#HEAD-PGL-%s\">Scroll down to \"Predicted gene locations\"</A>",gca->gsgmntn);
      CHECKALIGNLEN;
#ifdef PLANTGDB
      nextFree+=sprintf(global_algnmnt+nextFree, "\nEST sequence %6d +strand %5d n (File: %s)\n", calln, clgth, getDnaLink(2,cname));
#else
      nextFree+=sprintf(global_algnmnt+nextFree, "\nEST sequence %6d +strand %5d n (File: %s)\n", calln, clgth, getDnaLink(1,cname));
#endif
    }
    else {
      CHECKALIGNLEN; nextFree+=sprintf(global_algnmnt+nextFree, "\nEST sequence %6d +strand %5d n (File: %s)\n", calln, clgth, cname);
    }
  }
  else {
    if (htmlop) {
      CHECKALIGNLEN; nextFree+=sprintf(global_algnmnt+nextFree, "\n<A HREF=\"#HEAD-PGL-%s\">Scroll down to \"Predicted gene locations\"</A>",gca->gsgmntn);
      CHECKALIGNLEN;
#ifdef PLANTGDB
      nextFree+=sprintf(global_algnmnt+nextFree, "\nEST sequence %6d -strand %5d n (File: %s)\n", calln, clgth, getDnaLink(2,cname));
#else
      nextFree+=sprintf(global_algnmnt+nextFree, "\nEST sequence %6d -strand %5d n (File: %s)\n", calln, clgth, getDnaLink(1,cname));
#endif
    }
    else {
      CHECKALIGNLEN; nextFree+=sprintf(global_algnmnt+nextFree, "\nEST sequence %6d -strand %5d n (File: %s)\n", calln, clgth, cname);
    }
  }
  AppendChar('\0');

  prt_seq_segment_to_str(global_algnmnt, cdna, numbp, 0, clgth - 1, 0, NAUC, 1);
  nextFree=strlen(global_algnmnt);

  allocate_Dspace(glgth, clgth);

  JNLGTH = sprm.join_length;
  MIN_EXONLENGTH =  sprm.min_exon_length;
  SEWEIGHT = -100.0;
  MIN_INTRONLENGTH = sprm.min_intron_length;
  SIWEIGHT = -100.0;
  MIN_NBR_ENDMATCHES =  sprm.min_nbr_endmatches;
  if (MIN_NBR_ENDMATCHES > MIN_EXONLENGTH) MIN_NBR_ENDMATCHES = MIN_EXONLENGTH;


  /* ESTABLISH THE SCORING MATRIX: */
  InitMatch(sprm.ids, sprm.mms, sprm.nns, sprm.dls);

  /* FILL OUT THE STATE AND PATH MATRICES: */
  CompletePathMatrix(glgth, clgth, lgdna, cdna, lpd, lpa, sprm.pdg);

  /* BACK TRACE AN OPTIMAL PATH: */
  OptimalPath(glgth, clgth, lgdna, cdna);

  /* PRINT OUT AN OPTIMAL ALIGNMENT: */
  OutputEnds(ia);
  opflag = OutputParse(lpd, lpa, ia, numbp, rflag, gname, glgth, cdna, cname, clgth, 2*sprm.long_intron, gca);
  if (opflag > 0) {
    if (opflag == 1) OutputAlignment(numbp, rflag, gca);
    if (opflag == 2) {
      giaorig = gca->gia;   giborig = gca->gib;
      if (rflag == 1) {
        if (gca->ccds[0][0] > MIN_EXONLENGTH  &&  gca->gcds[0][0] + 2*sprm.long_intron > numbp - gca->gia) {
	  gca->gia-= 2*sprm.long_intron;
	  if (gca->gia < sia) gca->gia = sia;
	}
	if (gca->ccds[gca->exn-1][1] < clgth-MIN_EXONLENGTH  &&  numbp - gca->gib + 2*sprm.long_intron > gca->gcds[gca->exn-1][1]) {
	  gca->gib+= 2*sprm.long_intron;
	  if (gca->gib > sib) gca->gib = sib;
	}
      }
      if (rflag == 0) {
	if (gca->ccds[0][0] > MIN_EXONLENGTH  &&  gca->gia + 2*sprm.long_intron > gca->gcds[0][0]) {
	  gca->gia-= 2*sprm.long_intron;
	  if (gca->gia < sia) gca->gia = sia;
	}
	if (gca->ccds[gca->exn-1][1] < clgth-MIN_EXONLENGTH  &&  gca->gcds[gca->exn-1][1] + 2*sprm.long_intron > gca->gib) {
	  gca->gib+= 2*sprm.long_intron;
	  if (gca->gib > sib) gca->gib = sib;
	}
      }
      if ((gca->gia != giaorig  ||  gca->gib != giborig)  &&  gca->gib - gca->gia < MAXGLGTH) {
        free_Dspace(glgth, clgth);
	goto ADJUST;
      }
      else {
	gca->gia = giaorig;   gca->gib = giborig;   glgth = gca->gib - gca->gia + 1;
	OutputAlignment(numbp, rflag, gca);
      }
    }
  }
  if (gca->score < 0.0) {				/*THE FOLLOWING PREVENTS DISPLAY OF UNSUCCESSFUL ALIGNMENTS */
    free_Dspace(glgth, clgth);
    return (0);
  }
  gca->calln = calln;
  gca->next = NULL;

  free_Dspace(glgth, clgth);

#ifdef MPITIMING
MPIV_finish = MPI_Wtime();
printf("\nProcessor %2d used %9.6f seconds in sahmtD() for EST sequence %6d +strand %5d n (File: %s)",
                MPIV_rank, MPIV_finish - MPIV_start, calln, clgth, cname);
fflush(stdout);
#endif

  return (1);

}							/* end sahmtD() */



static void CompletePathMatrix(int glgth, int clgth, char *gdna, char *cdna, float *pd, float *pa, float pdg)
{
  int m, n;

  scoreD[E][0][0] = scoreD[E][1][0] = 0.0;
  scoreD[I][0][0] = scoreD[I][1][0] = 0.0;

  for (n = 0; n < 2; ++n)
    for (m = 0; m <= clgth; ++m) {
      exonstart[n][m] = 0;
      intronstart[n][m] = 0;
    }

/* FILL OUT THE N = 0 ROW OF THE PATH MATRIX: */
  for (m = 1; m <= clgth; m++) {
    E_0m(m);
    I_0m(m);
  }

/* COMPLETE THE REST OF THE PATH MATRIX FOR n > 0: */
  for (n = 1; n <= glgth; n++) {
    E_n0(n);
    I_n0(n);
    for (m = 1; m <= clgth; m++) {
      E_nm(n, m, gdna, cdna, glgth, clgth, pd, pa, pdg);
      I_nm(n, m, gdna, cdna, glgth, clgth, pd, pa, pdg);
    }
  }

}							/* end CompletePathMatrix() */



static void E_0m(int m)
{

  scoreD[E][0][m] = 0.0;
  pathD[E][0][m] = E_M;

  return;

}							/* end E_0m() */



static void I_0m(int m)
{

/* DISALLOW INTRON STATUS FOR 5' NON-MATCHING cDNA LETTERS: */
  scoreD[I][0][m] = -99999.;
  pathD[I][0][m] = I_M;

  return;

}							/* end I_0m() */



static void E_n0(int n)
{
  pathD[E][n][0] = E_N;
  return;

}							/* end E_n0() */



static void I_n0(int n)
{
  pathD[I][n][0] = I_N;
  return;

}							/* end I_n0() */



static void E_nm(int n, int m, char *gdna, char *cdna, int glgth, int clgth, float *pd, float *pa, float pdg)
{
  float x[6];
  int largest;
  int n_1=MOD2(n - 1),n_=MOD2(n);

  x[0] = scoreD[E][n_1][m - 1] + weight(E, E_NM, pd, pa, n, glgth, m, clgth, gdna[n - 1], cdna[m - 1], pdg);
  x[1] = scoreD[I][n_1][m - 1] + weight(E, I_NM, pd, pa, n, glgth, m, clgth, gdna[n - 1], cdna[m - 1], pdg);
  if (n - intronstart[n_1][m - 1] < MIN_INTRONLENGTH)
    x[1] += SIWEIGHT;

  x[2] = scoreD[E][n_1][m] + weight(E, E_N, pd, pa, n, glgth, m, clgth, gdna[n - 1], DASH, pdg);
  x[3] = scoreD[I][n_1][m] + weight(E, I_N, pd, pa, n, glgth, m, clgth, gdna[n - 1], DASH, pdg);
  if (n - intronstart[n_1][m] < MIN_INTRONLENGTH)
    x[3] += SIWEIGHT;

  x[4] = scoreD[E][n_][m - 1] + weight(E, E_M, pd, pa, n, glgth, m, clgth, DASH, cdna[m - 1], pdg);
  x[5] = scoreD[I][n_][m - 1] + weight(E, I_M, pd, pa, n, glgth, m, clgth, DASH, cdna[m - 1], pdg);
  if (n - intronstart[n_][m - 1] + 1 < MIN_INTRONLENGTH)
    x[5] += SIWEIGHT;

  largest = IndexMax(x, 6);
  scoreD[E][n_][m] = x[largest];
  pathD[E][n][m] = largest;

  if (largest == 1 || largest == 3 || largest == 5) {
    exonstart[n_][m] = n;
  }
  if (largest == 0) {
    exonstart[n_][m] = exonstart[n_1][m - 1];
  }
  if (largest == 2) {
    exonstart[n_][m] = exonstart[n_1][m];
  }
  if (largest == 4) {
    exonstart[n_][m] = exonstart[n_][m - 1];
  }
  return;

}							/* end E_nm() */



static void I_nm(int n, int m, char *gdna, char *cdna, int glgth, int clgth, float *pd, float *pa, float pdg)
{
  float x[2];
  int largest;
  int n_1=MOD2(n - 1),n_=MOD2(n);

  x[0] = scoreD[E][n_1][m] + weight(I, E_N, pd, pa, n, glgth, m, clgth, gdna[n - 1], STAR, pdg);
  if (n > 1  &&  n - exonstart[n_1][m] < MIN_EXONLENGTH)
    x[0] += SEWEIGHT;
  x[1] = scoreD[I][n_1][m] + weight(I, I_N, pd, pa, n, glgth, m, clgth, gdna[n - 1], STAR, pdg);

  largest = (x[0]>=x[1])?0:1;
  scoreD[I][n_][m] = x[largest];
  if (largest == 0) {
    pathD[I][n][m] = E_N;
    intronstart[n_][m] = n;
  }
  if (largest == 1) {
    pathD[I][n][m] = I_N;
    intronstart[n_][m] = intronstart[n_1][m];
  }

  return;

}							/* end I_nm() */



static void OptimalPath(int glgth, int clgth, char *gdna, char *cdna)
{
  int largest;
  float x[2];


  /*FIRST, FIND WHICH OF THE MATRICES HAS THE HIGHEST FINAL VALUE */
  x[0] = scoreD[E][MOD2(glgth)][clgth];
  x[1] = scoreD[I][MOD2(glgth)][clgth];

  largest = (x[0]>=x[1])?0:1;

  optN = glgth;
  optM = clgth;
  optcounter = 0;
  TracePath(largest, gdna, cdna);

}							/* end OptimalPath() */



static void TracePath(int statev, char *gdna, char *cdna)
{
  if (optN == 0 && optM == 0)
    return;

  switch (pathD[statev][optN][optM]) {
  case (E_NM):
    WritePath(gdna[optN - 1], cdna[optM - 1], statev);
    optN--;
    optM--;
    TracePath(E, gdna, cdna);
    return;
  case (I_NM):
    WritePath(gdna[optN - 1], cdna[optM - 1], statev);
    optN--;
    optM--;
    TracePath(I, gdna, cdna);
    return;
  case (E_N):
    if (statev == E) {
      WritePath(gdna[optN - 1], DASH, statev);
    }
    if (statev == I) {
      WritePath(gdna[optN - 1], STAR, statev);
    }
    optN--;
    TracePath(E, gdna, cdna);
    return;
  case (I_N):
    if (statev == E) {
      WritePath(gdna[optN - 1], DASH, statev);
    }
    if (statev == I) {
      WritePath(gdna[optN - 1], STAR, statev);
    }
    optN--;
    TracePath(I, gdna, cdna);
    return;
  case (E_M):
    WritePath(DASH, cdna[optM - 1], statev);
    optM--;
    TracePath(E, gdna, cdna);
    return;
  case (I_M):
    WritePath(DASH, cdna[optM - 1], statev);
    optM--;
    TracePath(I, gdna, cdna);
    return;
  }

}							/* end TracePath() */



static void WritePath(int g, int c, int statev)
{

  optGDNA[optcounter] = g;
  optCDNA[optcounter] = c;

  optState[optcounter] = statev;
  optcounter++;

  return;

}							/* end WritePath() */



static int Obeg, Gbeg, Cbeg, Oend;

static void OutputEnds(int ia)
{
  int i, em = 0, el = 0;

  /* ... alignment ends are trimmed to satisfy the following criteria:
     1) MIN_NBR_ENDMATCHES of non-NN consecutive matches at the ends
     2) no indels in the MIN_EXONLENGTH terminal positions 
  */

  Gbeg = ia;
  Cbeg = 0;

  for (i = 0; i < optcounter; i++) {
    if (optGDNA[optcounter - 1 - i] < 11) {
      ++Gbeg;
      if (optCDNA[optcounter - 1 - i] < 11) {
	++Cbeg;
	if (optCDNA[optcounter - 1 - i] != 10  &&  optCDNA[optcounter - 1 - i] == optGDNA[optcounter - 1 - i]) {
	  ++em;
	}
	else {
	  if (em < MIN_NBR_ENDMATCHES) {
	    em = el = 0;
	    continue;
	  }
	}
	if (++el >= MIN_EXONLENGTH) break;
      }
      else {
	em = el = 0;
      }
    }
    else {
      ++Cbeg;
      em = el = 0;
    }
  }
  Obeg = i - el + 1;
  Gbeg+= (-el+1);
  Cbeg+= (-el+1);

  el = 0;
  for (i = 0; i < optcounter; i++) {
    if (optGDNA[i] < 11) {
      if (optCDNA[i] < 11) {
	if (optCDNA[i] != 10  &&  optCDNA[i] == optGDNA[i]) {
	  ++em;
	}
	else {
	  if (em < MIN_NBR_ENDMATCHES) {
	    em = el = 0;
	    continue;
	  }
	}
	if (++el >= MIN_EXONLENGTH) break;
      }
      else {
	 em = el = 0;
      }
    }
    else {
	em = el = 0;
    }
  }
  Oend = optcounter - 1 - i + el - 1;

  if (Obeg == optcounter) {	/* ... non-matching degenerate alignment */
	--Obeg;
  }
  if (Oend < Obeg) {
	Oend = Obeg;
  }

}							/* end OutputEnds () */



static int OutputParse(float *pd, float *pa, int ia, int numbp, int rflag, char *gname, int glgth, char *cdna, char *cname, int clgth, int gextnsn, struct gcalgnmnt *gca)
{
  int i, j, exn = 0, intrA, WSIZE = 50;
  int currentState, GBegin, GEnd, CBegin, CEnd;
  int mlgth;
  int ppa, mma;
  float escr, oescr, gescr, ogescr, sscr, osscr, cvrge, simscore;

  mlgth = 0;
  escr = oescr = gescr = ogescr = simscore = cvrge = 0.0;
  if (rflag) {
    CHECKALIGNLEN; nextFree+=sprintf(global_algnmnt+nextFree, "\nPredicted gene structure (within gDNA segment %d to %d):\n\n",numbp - gca->gia,numbp - gca->gib);
  }
  else {
    CHECKALIGNLEN; nextFree+=sprintf(global_algnmnt+nextFree, "\nPredicted gene structure (within gDNA segment %d to %d):\n\n",gca->gia + 1,gca->gib + 1);
  }
  strcpy(gca->gname, gname);
  strcpy(gca->cname, cname);

  if (Oend == Obeg) {
    CHECKALIGNLEN; nextFree+=sprintf(global_algnmnt+nextFree, "\nNOMATCH\t%s\t%s\t%5.3f\t%5.3f\t%d\n", gname, cname, gescr, ogescr, mlgth);
    gca->exn = 0;
    gca->score = -99.9;
    return (0);
  }

  currentState = optState[optcounter - 1 - Obeg];
  GBegin = GEnd = Gbeg;
  CBegin = CEnd = Cbeg;
  intrA = Obeg;
  if (optCDNA[optcounter - 1 - Obeg] < STAR) {
    if (optGDNA[optcounter - 1 - Obeg] == optCDNA[optcounter - 1 - Obeg] || optGDNA[optcounter - 1 - Obeg] == DASH)
      escr += GetScore(optGDNA[optcounter - 1 - Obeg], optCDNA[optcounter - 1 - Obeg]);
    oescr += GetScore(optGDNA[optcounter - 1 - Obeg], optGDNA[optcounter - 1 - Obeg]);
  }
  for (i = Obeg + 1; i <= Oend; i++) {
    if (optState[optcounter - 1 - i] != currentState) {
      if (currentState == E) {
	intrA = optcounter - 1 - i;
	if (rflag) {
	  gca->gcds[exn][0] = numbp - GBegin + 1;
	  gca->gcds[exn][1] = numbp - GEnd + 1;
	  CHECKALIGNLEN; nextFree+=sprintf(global_algnmnt+nextFree, " Exon %2d %*d %*d ", ++exn, ifwdth, numbp - GBegin + 1, ifwdth, numbp - GEnd + 1);
	}
	else {
	  gca->gcds[exn][0] = GBegin;
	  gca->gcds[exn][1] = GEnd;
	  CHECKALIGNLEN; nextFree+=sprintf(global_algnmnt+nextFree, " Exon %2d %*d %*d ", ++exn, ifwdth, GBegin, ifwdth, GEnd);
	}
	gca->ccds[exn - 1][0] = CBegin;
	gca->ccds[exn - 1][1] = CEnd;
	CHECKALIGNLEN; nextFree+=sprintf(global_algnmnt+nextFree,
	  "(%4d n);  cDNA %6d %6d (%4d n);", GEnd - GBegin + 1, CBegin, CEnd, CEnd - CBegin + 1);
	if (oescr > 0) {
	  CHECKALIGNLEN; nextFree+=sprintf(global_algnmnt+nextFree, " score: %5.3f\n", escr / oescr);
	  gca->exnscr[exn - 1] = escr / oescr;
	}
	else {
	  CHECKALIGNLEN; nextFree+=sprintf(global_algnmnt+nextFree, " score:     0\n");
	  gca->exnscr[exn - 1] = 0.0;
	}
	if (GEnd - GBegin + 1 >= MINELGTHfS) {
	  gescr += escr;
	  ogescr += oescr;
	}
	mlgth += GEnd - GBegin + 1;
	CBegin = CEnd + 1;
      }
      else {
	if (rflag) {
	  CHECKALIGNLEN; nextFree+=sprintf(global_algnmnt+nextFree, "  Intron %2d %*d %*d ", exn, ifwdth, numbp - GBegin + 1, ifwdth, numbp - GEnd + 1);
	}
	else {
	  CHECKALIGNLEN; nextFree+=sprintf(global_algnmnt+nextFree, "  Intron %2d %*d %*d ", exn, ifwdth, GBegin, ifwdth, GEnd);
	}
	CHECKALIGNLEN; nextFree+=sprintf(global_algnmnt+nextFree, "(%4d n); ", GEnd - GBegin + 1);
	CHECKALIGNLEN; nextFree+=sprintf(global_algnmnt+nextFree, "Pd: %5.3f", pd[GBegin - 1 - ia]);
	gca->itrscr[exn - 1][0] = pd[GBegin - 1 - ia];
	sscr = osscr = 0.0;
	for (j = 1; j <= WSIZE; ++j) {
	  if (intrA + j > optcounter - 1 - Obeg || optState[intrA + j] == I)
	    break;
	  if (optGDNA[intrA + j] == optCDNA[intrA + j] || optGDNA[intrA + j] == DASH)
	    sscr += GetScore(optGDNA[intrA + j], optCDNA[intrA + j]);
	  osscr += GetScore(optGDNA[intrA + j], optGDNA[intrA + j]);
	}
	if (j > (int) (.8 * WSIZE)  &&  osscr > 0.0) {
	  CHECKALIGNLEN; nextFree+=sprintf(global_algnmnt+nextFree, " (s: %4.2f), ", sscr / osscr);
	}
	else {
	  CHECKALIGNLEN; nextFree+=sprintf(global_algnmnt+nextFree, " (s:    0), ");
	}
	CHECKALIGNLEN; nextFree+=sprintf(global_algnmnt+nextFree, "Pa: %5.3f", pa[GEnd - 1 - ia]);
	gca->itrscr[exn - 1][1] = pa[GEnd - 1 - ia];
	sscr = osscr = 0.0;
	for (j = 1; j <= WSIZE; ++j) {
	  if (i + j > Oend + 1 || optState[optcounter - i - j] == I)
	    break;
	  if (optGDNA[optcounter - i - j] == optCDNA[optcounter - i - j] || optGDNA[optcounter - i - j] == DASH)
	    sscr += GetScore(optGDNA[optcounter - i - j], optCDNA[optcounter - i - j]);
	  osscr += GetScore(optGDNA[optcounter - i - j], optGDNA[optcounter - i - j]);
	}
	if (j > (int) (.8 * WSIZE)  &&  osscr > 0.0) {
	  CHECKALIGNLEN; nextFree+=sprintf(global_algnmnt+nextFree, " (s: %4.2f)", sscr / osscr);
	}
	else {
	  CHECKALIGNLEN; nextFree+=sprintf(global_algnmnt+nextFree, " (s:    0)");
	}
	if (GEnd - GBegin + 1 <= MIN_INTRONLENGTH) {
	  CHECKALIGNLEN; nextFree+=sprintf(global_algnmnt+nextFree, " ??");
	}
	AppendChar('\n');;
      }
      GBegin = GEnd + 1;
      escr = oescr = 0.0;
    }
    if (optCDNA[optcounter - 1 - i] < STAR) {
      if (optGDNA[optcounter - 1 - i] == optCDNA[optcounter - 1 - i] || optGDNA[optcounter - 1 - i] == DASH)
	escr += GetScore(optGDNA[optcounter - 1 - i], optCDNA[optcounter - 1 - i]);
      oescr += GetScore(optGDNA[optcounter - 1 - i], optGDNA[optcounter - 1 - i]);
    }
    currentState = optState[optcounter - 1 - i];
    if (optGDNA[optcounter - 1 - i] < 11)
      ++GEnd;
    if (optCDNA[optcounter - 1 - i] < 11)
      ++CEnd;
  }
  if (currentState == E) {
    if (rflag) {
      gca->gcds[exn][0] = numbp - GBegin + 1;
      gca->gcds[exn][1] = numbp - GEnd + 1;
      CHECKALIGNLEN; nextFree+=sprintf(global_algnmnt+nextFree, " Exon %2d %*d %*d ", ++exn, ifwdth, numbp - GBegin + 1, ifwdth, numbp - GEnd + 1);
    }
    else {
      gca->gcds[exn][0] = GBegin;
      gca->gcds[exn][1] = GEnd;
      CHECKALIGNLEN; nextFree+=sprintf(global_algnmnt+nextFree, " Exon %2d %*d %*d ", ++exn, ifwdth, GBegin, ifwdth, GEnd);
    }
    gca->ccds[exn - 1][0] = CBegin;
    gca->ccds[exn - 1][1] = CEnd;
    CHECKALIGNLEN; nextFree+=sprintf(global_algnmnt+nextFree, "(%4d n);  cDNA %6d %6d (%4d n);", GEnd - GBegin + 1, CBegin, CEnd, CEnd - CBegin + 1);
    if (oescr > 0) {
      CHECKALIGNLEN; nextFree+=sprintf(global_algnmnt+nextFree, " score: %5.3f\n", escr / oescr);
      gca->exnscr[exn - 1] = escr / oescr;
    }
    else {
      CHECKALIGNLEN; nextFree+=sprintf(global_algnmnt+nextFree, " score:     0\n");
      gca->exnscr[exn - 1] = 0.0;
    }
    if (GEnd - GBegin + 1 >= MINELGTHfS) {
      gescr += escr;
      ogescr += oescr;
    }
    mlgth += GEnd - GBegin + 1;
  }
  else {
    if (rflag) {
      CHECKALIGNLEN; nextFree+=sprintf(global_algnmnt+nextFree, "  Intron %2d %*d %*d ", exn, ifwdth, numbp - GBegin + 1, ifwdth, numbp - GEnd + 1);
    }
    else {
      CHECKALIGNLEN; nextFree+=sprintf(global_algnmnt+nextFree, "  Intron %2d %*d %*d ", exn, ifwdth, GBegin, ifwdth, GEnd);
    }
    CHECKALIGNLEN; nextFree+=sprintf(global_algnmnt+nextFree, "(%4d n); ", GEnd - GBegin + 1);
    CHECKALIGNLEN; nextFree+=sprintf(global_algnmnt+nextFree, "Pd: %5.3f", pd[GBegin - 1 - ia]);
    gca->itrscr[exn - 1][0] = pd[GBegin - 1 - ia];
    sscr = osscr = 0.0;
    for (j = 1; j <= WSIZE; ++j) {
      if (intrA + j > optcounter - 1 - Obeg || optState[intrA + j] == I)
	break;
      if (optGDNA[intrA + j] == optCDNA[intrA + j] || optGDNA[intrA + j] == DASH)
	sscr += GetScore(optGDNA[intrA + j], optCDNA[intrA + j]);
      osscr += GetScore(optGDNA[intrA + j], optGDNA[intrA + j]);
    }
    if (j > (int) (.8 * WSIZE)) {
      CHECKALIGNLEN; nextFree+=sprintf(global_algnmnt+nextFree, " (s: %4.2f), ", sscr / osscr);
    }
    else {
      CHECKALIGNLEN; nextFree+=sprintf(global_algnmnt+nextFree, " (s:    0), ");
    }
    CHECKALIGNLEN; nextFree+=sprintf(global_algnmnt+nextFree, "Pa: %5.3f", pa[GEnd - 1 - ia]);
    gca->itrscr[exn - 1][1] = pa[GEnd - 1 - ia];
    sscr = osscr = 0.0;
    for (j = 1; j <= WSIZE; ++j) {
      if (i + j > Oend + 1 || optState[optcounter - i - j] == I)
	break;
      if (optGDNA[optcounter - i - j] == optCDNA[optcounter - i - j] || optGDNA[optcounter - i - j] == DASH)
	sscr += GetScore(optGDNA[optcounter - i - j], optCDNA[optcounter - i - j]);
      osscr += GetScore(optGDNA[optcounter - i - j], optGDNA[optcounter - i - j]);
    }
    if (j > (int) (.8 * WSIZE)) {
      CHECKALIGNLEN; nextFree+=sprintf(global_algnmnt+nextFree, " (s: %4.2f)", sscr / osscr);
    }
    else {
      CHECKALIGNLEN; nextFree+=sprintf(global_algnmnt+nextFree, " (s:    0)");
    }
    if (GEnd - GBegin + 1 <= MIN_INTRONLENGTH) {
      CHECKALIGNLEN; nextFree+=sprintf(global_algnmnt+nextFree, " ??");
    }
    AppendChar( '\n');
  }

  gca->ppa[0] = 0;
  gca->ppa[1] = 0;
  ppa = mma = 0;
  for (i = gca->ccds[exn - 1][1] + 1 >= clgth - 50 ? gca->ccds[exn - 1][1] + 1 : clgth - 50; i <= clgth; ++i) {
    if (cdna[i - 1] == 2)
      ++ppa;
    else {
      if (ppa > 0 && mma < 1) {
	++mma;
	continue;
      }
      else {
	if (ppa >= 10)
	  break;
	else {
	  ppa = mma = 0;
	  continue;
	}
      }
    }
  }
  if (ppa >= 10) {
    CHECKALIGNLEN;
    nextFree+=sprintf(global_algnmnt+nextFree, " PPA     ");
    for (j=0;j<2*ifwdth+13;++j) nextFree+=sprintf(global_algnmnt+nextFree, " ");
    nextFree+=sprintf(global_algnmnt+nextFree, "cDNA %6d %6d\n", i - ppa - mma, i - 1);
    gca->ppa[0] = i - ppa - mma;
    gca->ppa[1] = i - 1;
  }
  else {
    ppa = mma = 0;
    for (i = gca->ccds[0][0] - 1 <= 50 ? gca->ccds[0][0] - 1 : 50; i >= 1; --i) {
      if (cdna[i - 1] == 0)
	++ppa;
      else {
	if (ppa > 0 && mma < 1) {
	  ++mma;
	  continue;
	}
	else {
	  if (ppa >= 10)
	    break;
	  else {
	    ppa = mma = 0;
	    continue;
	  }
	}
      }
    }
    if (ppa >= 10) {
      CHECKALIGNLEN;
      nextFree+=sprintf(global_algnmnt+nextFree, " PPA     ");
      for (j=0;j<2*ifwdth+13;++j) nextFree+=sprintf(global_algnmnt+nextFree, " ");
      nextFree+=sprintf(global_algnmnt+nextFree, "cDNA %6d %6d\n", i + ppa + mma, i + 1);
      gca->ppa[0] = i + ppa + mma;
      gca->ppa[1] = i + 1;
    }
  }

  cvrge = (float) mlgth / (float) glgth > (float) mlgth / (float) clgth ?
    (float) mlgth / (float) glgth : (float) mlgth / (float) clgth;
  if (ogescr > 0.0  ||  cvrge >= MINCVRGEfS) {
    if (ogescr > 0.0) simscore = gescr / ogescr;
    CHECKALIGNLEN; nextFree+=sprintf(global_algnmnt+nextFree, "\nMATCH\t%s\t%s\t%5.3f\t%d\t%5.3f", gname, cname, simscore, mlgth, cvrge);
    if ((float) mlgth / (float) glgth > (float) mlgth / (float) (clgth)) {
      CHECKALIGNLEN; nextFree+=sprintf(global_algnmnt+nextFree, "\tG\n");
    }
    else {
      CHECKALIGNLEN; nextFree+=sprintf(global_algnmnt+nextFree, "\tC\n");
    }
    CHECKALIGNLEN; nextFree+=sprintf(global_algnmnt+nextFree, "PGS_%s_%s\t(", gname, cname);
    for (i = 0; i < exn - 1; ++i) {
      CHECKALIGNLEN; nextFree+=sprintf(global_algnmnt+nextFree, "%d  %d,", gca->gcds[i][0], gca->gcds[i][1]);
    }
    CHECKALIGNLEN; nextFree+=sprintf(global_algnmnt+nextFree, "%d  %d)\n", gca->gcds[exn - 1][0], gca->gcds[exn - 1][1]);
  }
  else {
    CHECKALIGNLEN; nextFree+=sprintf(global_algnmnt+nextFree, "\nNOMATCH\t%s\t%s\t%5.3f\t%5.3f\t%d\n", gname, cname, gescr, ogescr, mlgth);
    gca->exn = 0;
    gca->score = -99.9;
    return (0);
  }

  gca->exn = exn;
  gca->score = simscore;

  /* return 2:	The alignment does not exhaust the cDNA sequence but extends almost to the borders
		of the gDNA segment.  In this case, the alignment will be re-calculated in a larger
		gDNA segment (extended by gextnsn at the close borders) to see whether a more
		complete alignment can be obtained that was missed previously.
  */
  if (rflag == 1) {
    if (gca->ccds[0][0] > MIN_EXONLENGTH  &&  gca->gia > 0  &&  gca->gcds[0][0] + gextnsn > numbp - gca->gia) return (2); 
    if (gca->ccds[gca->exn-1][1] < clgth-MIN_EXONLENGTH  &&  gca->gib < numbp-1  &&  numbp - gca->gib + gextnsn > gca->gcds[gca->exn-1][1]) return (2);
  }
  if (rflag == 0) {
    if (gca->ccds[0][0] > MIN_EXONLENGTH  &&  gca->gia > 0  &&  gca->gia + gextnsn > gca->gcds[0][0]) return (2); 
    if (gca->ccds[gca->exn-1][1] < clgth-MIN_EXONLENGTH  &&  gca->gib < numbp-1  &&  gca->gcds[gca->exn-1][1] + gextnsn > gca->gib) return (2);
  }

  return (1);

}							/* end OutputParse() */



static void OutputAlignment(int numbp, int rflag, struct gcalgnmnt *gca)
{
  int i, j;
  int base = 0, cbase, gcount = Gbeg - 1, ccount = Cbeg - 1;

  if (LTFLAG) {
    Obeg = 0;
    Oend = optcounter - 1;
  }

  CHECKALIGNLEN; nextFree+=sprintf(global_algnmnt+nextFree, "\nAlignment (genomic DNA sequence = upper lines):\n\n");

  for (i = Obeg; i <= Oend; i++) {
    OutputBase(optGDNA[optcounter - 1 - i]);
    base++;
    if (optGDNA[optcounter - 1 - i] < 11)
      gcount++;
    if (base % 60 == 0 || i == Oend) {
      if (rflag) {
	CHECKALIGNLEN; nextFree+=sprintf(global_algnmnt+nextFree, "   %*d\n", ifwdth, numbp - gcount + 1);
      }
      else {
	CHECKALIGNLEN; nextFree+=sprintf(global_algnmnt+nextFree, "   %*d\n", ifwdth, gcount);
      }
      cbase = 0;
      for (j = i - base + 1; j <= i; ++j) {
	if (optCDNA[optcounter - 1 - j] == optGDNA[optcounter - 1 - j]) {
	  AppendChar('|');
	}
	else {
	  AppendChar(' ');
	}
	cbase++;
	if (cbase % 10 == 0)
	  AppendChar(' ');
      }
      AppendChar('\n');;
      cbase = 0;
      for (j = i - base + 1; j <= i; ++j) {
	OutputBase(optCDNA[optcounter - 1 - j]);
	if (optCDNA[optcounter - 1 - j] < 11)
	  ccount++;
	cbase++;
	if (cbase % 10 == 0)
	  AppendChar(' ');
      }
      if (j == Oend + 1 && cbase % 10 != 0)
	AppendChar(' ');
      CHECKALIGNLEN; nextFree+=sprintf(global_algnmnt+nextFree, "  %*d\n\n\n", ifwdth, ccount);
      base = 0;
    }
    else if (base % 10 == 0)
      AppendChar(' ');
  }
  AppendChar('\0');
  if ((gca->algnmnt = (char *) calloc(nextFree+1, sizeof(char))) == NULL)
    fatal_error("Error: memory allocation failed. Exit.\n");
  if (LTFLAG) fprintf(stdout,"\nALLOC_ALN1 %ld # %d",(long)(gca->algnmnt),++NALLOC_ALN);
  strcpy(gca->algnmnt, global_algnmnt);

}							/* end OutputAlignment() */



static void InitMatch(float ids, float mms, float nns, float dls)
{
  int n, m;

/*FILL OUT MATRIX FOR MATCHING ACTUAL BASES T, C, A, AND G  */
  for (n = 0; n < 4; n++) {
    for (m = 0; m < 4; m++) {
      if (m == n)
	match[n][m] = ids;
      else
	match[n][m] = mms;
    }
  }

/*FILL OUT MATRIX FOR MATCHING T,C,A,G WITH N */
  for (n = 0; n < 4; n++) {
    for (m = 4; m < 11; m++) {
      match[n][m] = nns;
    }
  }

/*FILL OUT MATRIX FOR MATCHING T,C,A,G WITH DASH */
  for (n = 0; n < 4; n++) {
    match[n][11] = dls;
  }

/*FILL OUT MATRIX FOR MATCHING N WITH T,C,A,G */
  for (n = 4; n < 11; n++) {
    for (m = 0; m < 4; m++) {
      match[n][m] = nns;
    }
  }

/*FILL OUT MATRIX FOR MATCHING N WITH N */
  for (n = 4; n < 11; n++) {
    for (m = 4; m < 11; m++) {
      match[n][m] = nns;
    }
  }

/*FILL OUT MATRIX FOR MATCHING N WITH DASH */
  for (n = 4; n < 11; n++) {
    match[n][11] = dls;
  }

/*FILL OUT MATRIX FOR MATCHING DASH WITH N */
  for (m = 0; m < 11; m++) {
    match[11][m] = dls;
  }

  match[11][11] = 0.0;

}							/* end InitMatch() */



static float GetScore(int g, int c)
{
  float score;

  score = match[g][c];
  return (score);

}							/* end GetScore() */



static int IndexMax(float *x, int n)
{
  int i, imax;
  float maxscore;


  maxscore = x[0];
  imax = 0;
  for (i = 1; i < n; i++) {
    if (x[i] > maxscore) {
      maxscore = x[i];
      imax = i;
    }
  }
  return (imax);
/*
   IndexMax RETURNS OF THE n INDECES FOR THE X ARRAY THE INDEX imax CORRESPONDING
   TO THE ARRAY ELEMENT OF LARGEST VALUE.
 */

}							/* end IndexMax() */



static double weight(int state, int path, float *pd, float *pa, int n,
  int glgth, int m, int clgth, int g, int c, float pdg)
{
  double rval = 0.0;

  if (state == E) {
    /* Transition weights: */

    if (path == E_N) {
      if (n == 1)
	rval += Log_Half;
      else if (m < clgth || n < 80)
	rval += (Log_1_Pdg + Log_1_Pd[n-1]);
    }
    else if (path == E_NM) {
      if (n == 1)
	rval += Log_Half;
      else
	rval += (Log_1_Pdg + Log_1_Pd[n-1]);
    }
    else if (path == E_M) {
      if (n == 1)
	rval += Log_Half;
      else if (n < glgth || m < 80)
	rval += Log_Pdg;
    }
    else if (path == I_N || path == I_NM) {
      if (n == 1)
	rval += Log_Half;
      else
	rval += (Log_Pa[n-2] + Log_1_Pdg);
    }
    else {						/* path==I_M */
      if (n == 1)
	rval += Log_Half;
      else if (n < glgth)
	rval += (Log_Pa[n-1] + Log_Pdg);
    }

    /* Output weights: */

    if (path == E_N || path == I_N) {
      if (m < clgth)
	rval += match[g][c];
    }
    else if (path == E_NM || path == I_NM) {
      rval += match[g][c];
      if ((m < 80 || m > clgth - 80) && g == c)
	rval -= (match[g][c] / 2.0);
    }
    else {						/* path==E_M || path==I_M */
      if (n < glgth)
	rval += match[g][c];
    }
  }

  else {						/* state == I */
    if (n == 1)
      rval += Log_Half;
    else {
      if (path == E_N)
	    rval += (Log_1_Pdg + Log_Pd[n-1]);
      else						/* path==I_N */
	if (m < clgth)   rval +=  Log_1_Pa[n-2];
    }
  }

  return (rval);

}							/* end weight() */



static char numtonu[]="TCAGYRSWKMN-??????????????????.";

static void OutputBase(int base)
{
  if(base<0 || base >30) {
    CHECKALIGNLEN; nextFree+=sprintf(global_algnmnt+nextFree, "\n\n  Unknown symbol in OutputBase function! \n\n");
  }
  else {
    AppendChar(numtonu[base]);
  }

}							/* end OutputBase() */



static void allocate_Dspace(int glgth, int clgth)
{
  int i, j;

  for (i = 0; i < 2; ++i) {
    if ((exonstart[i] = (int *) calloc((clgth + 1) , sizeof(int))) == NULL)
       fatal_error("Error: memory allocation failed. Exit.\n");
    if ((intronstart[i] = (int *) calloc((clgth + 1) , sizeof(int))) == NULL)
       fatal_error("Error: memory allocation failed. Exit.\n");
    if ((pathD[i] = (short int **) calloc((glgth + 1) , sizeof(short int *))) == NULL)
       fatal_error("Error: memory allocation failed. Exit.\n");
    for (j = 0; j < glgth + 1; ++j)
      if ((pathD[i][j] = (short int *) calloc((clgth + 1) , sizeof(short int))) == NULL)
	 fatal_error("Error: memory allocation failed. Exit.\n");
    for (j = 0; j < 2; ++j)
      if ((scoreD[i][j] = (float *) calloc((clgth + 1) , sizeof(float))) == NULL)
	 fatal_error("Error: memory allocation failed. Exit.\n");
  }
  if ((optGDNA = (int *) calloc((3 * glgth + clgth) , sizeof(int))) == NULL)
     fatal_error("Error: memory allocation failed. Exit.\n");
  if ((optCDNA = (int *) calloc((3 * glgth + clgth) , sizeof(int))) == NULL)
     fatal_error("Error: memory allocation failed. Exit.\n");
  if ((optState = (int *) calloc((3 * glgth + clgth) , sizeof(int))) == NULL)
     fatal_error("Error: memory allocation failed. Exit.\n");

}							/* end allocate_Dspace() */



static void free_Dspace(int glgth, int clgth)
{
  int i, j;

  for (i = 0; i < 2; ++i) {
    free(exonstart[i]);
    free(intronstart[i]);
    for (j = 0; j < glgth + 1; ++j)
      free(pathD[i][j]);
    free(pathD[i]);
    for (j = 0; j < 2; ++j)
      free(scoreD[i][j]);
  }
  free(optGDNA);
  free(optCDNA);
  free(optState);

}							/* end free_Dspace() */



void read_sahmt_prm(FILE *fp, struct sprmtr *sprm)
{
  char buf[257];

  fscanf(fp, "%s", buf);
  fscanf(fp, "%f", &(sprm->pdg));
  fscanf(fp, "%s", buf);
  fscanf(fp, "%f", &(sprm->ids));
  fscanf(fp, "%s", buf);
  fscanf(fp, "%f", &(sprm->mms));
  fscanf(fp, "%s", buf);
  fscanf(fp, "%f", &(sprm->nns));
  fscanf(fp, "%s", buf);
  fscanf(fp, "%f", &(sprm->dls));
  fscanf(fp, "%s", buf);
  fscanf(fp, "%d", &(sprm->unfrm));
  fscanf(fp, "%s", buf);
  fscanf(fp, "%d", &(sprm->min_intron_length));
  fscanf(fp, "%s", buf);
  fscanf(fp, "%d", &(sprm->min_exon_length));
  fscanf(fp, "%s", buf);
  fscanf(fp, "%d", &(sprm->min_nbr_endmatches));
  fscanf(fp, "%s", buf);
  fscanf(fp, "%d", &(sprm->tiny_exon));
  fscanf(fp, "%s", buf);
  fscanf(fp, "%d", &(sprm->short_exon));
  fscanf(fp, "%s", buf);
  fscanf(fp, "%d", &(sprm->long_intron));
  fscanf(fp, "%s", buf);
  fscanf(fp, "%f", &(sprm->poor_exon_score));
  fscanf(fp, "%s", buf);
  fscanf(fp, "%f", &(sprm->poor_donor_score));
  fscanf(fp, "%s", buf);
  fscanf(fp, "%f", &(sprm->poor_acptr_score));
  fscanf(fp, "%s", buf);
  fscanf(fp, "%d", &(sprm->join_length));

}							/* end read_sahmt_prm() */



int insert_gca(gcahpp, nodep)
  struct gcalgnmnt **gcahpp, *nodep;
{
int c;
struct gcalgnmnt *tmp;

  if (*gcahpp == NULL) {
	*gcahpp = nodep;
	nodep->next = NULL;
  }
  else {
	while ((*gcahpp)->next != NULL)
	 {c = compare_gcas(*gcahpp,nodep);
	  if (c == 0) return(0);
	  if (c == 1) {
	    nodep->next = (*gcahpp)->next;
	    tmp = *gcahpp;
	    *gcahpp= nodep;
	    free_gca(tmp);
	    return(1);
	  }
	  if (c == 2)
            gcahpp = &((*gcahpp)->next);
	  else
	    break;
         }
        if ((*gcahpp)->next == NULL) {
        	c = compare_gcas(*gcahpp,nodep);
		if (c == 0) return(0);
		if (c == 1) {
		  nodep->next = NULL;
		  tmp = *gcahpp;
		  *gcahpp= nodep;
		  free_gca(tmp);
		  return(1);
		}
		if (c == 2) {(*gcahpp)->next = nodep; nodep->next = NULL;}
		else
		  {nodep->next = *gcahpp; *gcahpp = nodep;}
          }
        else
          { nodep->next = *gcahpp; *gcahpp = nodep; }
  }
  return (1);

}							/* end insert_gca() */



int compare_gcas(gca1p, gca2p)
  struct gcalgnmnt *gca1p, *gca2p;
{
	/* return 0: discard gca2p           */
	/* return 1: replace gca1p by gca2p  */
	/* return 2: gca1p "<" gca2p         */
	/* return 3: gca1p ">" gca2p         */

int beg1 = MIN(gca1p->gcds[0][0],gca1p->gcds[gca1p->exn - 1][1]);
int end1 = MAX(gca1p->gcds[0][0],gca1p->gcds[gca1p->exn - 1][1]);
int beg2 = MIN(gca2p->gcds[0][0],gca2p->gcds[gca2p->exn - 1][1]);
int end2 = MAX(gca2p->gcds[0][0],gca2p->gcds[gca2p->exn - 1][1]);

  if ( strcmp(gca1p->gname, gca2p->gname) == 0   &&
       strcmp(gca1p->cname, gca2p->cname) == 0     ) {
    if ( beg1 <= beg2  &&  end1 >= end2 )  return(0);
    if ( beg2 <= beg1  &&  end2 >= end1 )  return(1);
  }

  if (is_xwab(end2,beg1,end1,JNLGTH)             ||
      is_xwab(beg2,beg1,end1,JNLGTH)             ||
      xy_contains_ab(beg2,end2,beg1,end1,JNLGTH)   ) {
	/* ... overlap -> multi-exon gca's before singlets */
    if (gca2p->exn == 1  &&  gca1p->exn > 1)  return(2);
    if (gca1p->exn == 1  &&  gca2p->exn > 1)  return(3);
  }
  if ( (beg1 < beg2) || (beg1 == beg2  &&  end1 >= end2) )  return(2);

  return(3);

}							/* end compare_gcas() */


void free_gca(struct gcalgnmnt *gca)
{
  if (gca->algnmnt != NULL) {
    if (LTFLAG) fprintf(stdout,"\nDEALLOC_ALN1 %ld # %d",(long)(gca->algnmnt),++NDEALLOC_ALN);
    free((char *)gca->algnmnt);
  }
  if (LTFLAG) fprintf(stdout,"\nDEALLOC_GCA %ld # %d",(long)gca,++NDEALLOC_GCA);
  free((struct gcalgnmnt *) gca);
}							/* end free_gca() */



int rev_sgl(FILE *fp, struct gcalgnmnt *gca, char *gdna, char *gdnaR, int sia, int sib, int numbp, float *pd, float *pa, float *pdR, float *paR, struct sprmtr sprm)
{
int i, npgl = 0, r = 0;
struct gcalgnmnt *leftgca = NULL;
int leftend = 0, rightend = 0;
int fwdni = 0, fwdppa = 0, revni = 0, revppa = 0, fwdmol = 0, revmol = 0;
float fwdpdpa = 0.0, revpdpa = 0.0;

  if (gca == NULL) return(0);

  if (LTFLAG) fprintf(stdout,"\nrev_sgl:\n");
  if (gca->exn > 0) {		/* ... safeguard against severely quality-adjusted gca ! */
    leftgca = gca;
    leftend  = MIN(gca->gcds[0][0],gca->gcds[gca->exn-1][1]);
    rightend = MAX(gca->gcds[0][0],gca->gcds[gca->exn-1][1]);
    if (gca->exn > 1) {
      if (gca->gcds[0][0] < gca->gcds[gca->exn-1][1]) {
	for (i = 0; i < gca->exn - 1; ++i) {
	  fwdpdpa += gca->itrscr[i][0];
	  fwdpdpa += gca->itrscr[i][1];
	  ++fwdni;
	}
      }
      else {
	for (i = 0; i < gca->exn - 1; ++i) {
	  revpdpa += gca->itrscr[i][0];
	  revpdpa += gca->itrscr[i][1];
	  ++revni;
	}
      }
    }
    if (gca->ppa[0] < gca->ppa[1]) {
      if (gca->gcds[0][0] < gca->gcds[gca->exn-1][1]) ++fwdppa;
      else                                            ++revppa;
    }
  }
  if (LTFLAG) fprintf(stdout,"\nDETPGL leftend %d rightend %d fwdni %d revni %d gca %s",leftend,rightend,fwdni,revni,gca->cname);

  gca = gca->next;
  while (gca != NULL) {
    if (gca->exn < 1) {		/* ... safeguard against severely quality-adjusted gca ! */
      gca = gca->next;
      continue;
    }
    if (is_leftend_pgl(gca,rightend,JNLGTH,fwdni,revni,fwdpdpa,revpdpa)) {
      ++npgl;
      if (LTFLAG) fprintf(stdout,"\nDETPGL npgl %d fwdni %d revni %d",npgl,fwdni,revni);
      if (fwdni+revni > 0) {
	if (fwdpdpa >= revpdpa) {
	  r+= orient_gca(fp,npgl,leftgca,gca,0,gdna,gdnaR,sia,sib,numbp,pd,pa,pdR,paR,sprm);
	  if (LTFLAG) fprintf(stdout,"\nPGL %3d (+ strand): %d to %d (f %d r %d\tfp %d rp %d)\n",npgl,leftend,rightend,fwdni,revni,fwdppa,revppa);
	}
        else {
	  r+= orient_gca(fp,npgl,leftgca,gca,1,gdna,gdnaR,sia,sib,numbp,pd,pa,pdR,paR,sprm);
	  if (LTFLAG) fprintf(stdout,"\nPGL %3d (- strand): %d to %d (f %d r %d\tfp %d rp %d)\n",npgl,rightend,leftend,fwdni,revni,fwdppa,revppa);
	}
      }
      else if (fwdppa+revppa > 0) {
	if (fwdppa >= revppa) {
	  r+= orient_gca(fp,npgl,leftgca,gca,0,gdna,gdnaR,sia,sib,numbp,pd,pa,pdR,paR,sprm);
	  if (LTFLAG) fprintf(stdout,"\nPGL %3d (+ strand): %d to %d (f %d r %d\tfp %d rp %d)\n",npgl,leftend,rightend,fwdni,revni,fwdppa,revppa);
	}
        else {
	  r+= orient_gca(fp,npgl,leftgca,gca,1,gdna,gdnaR,sia,sib,numbp,pd,pa,pdR,paR,sprm);
	  if (LTFLAG) fprintf(stdout,"\nPGL %3d (- strand): %d to %d (f %d r %d\tfp %d rp %d)\n",npgl,rightend,leftend,fwdni,revni,fwdppa,revppa);
	}
      }
      else {
	fwdmol = maxorflgth(gdna+leftend-1,rightend-leftend+1);
	revmol = maxorflgth(gdnaR+numbp-rightend,rightend-leftend+1);
	if (fwdmol >= revmol) {
	  r+= orient_gca(fp,npgl,leftgca,gca,0,gdna,gdnaR,sia,sib,numbp,pd,pa,pdR,paR,sprm);
	  if (LTFLAG) fprintf(stdout,"\nPGL %3d (+ strand): %d to %d (f %d r %d\tfp %d rp %d\tfs %d rs %d)\n",npgl,leftend,rightend,fwdni,revni,fwdppa,revppa,fwdmol,revmol);
	}
        else {
	  r+= orient_gca(fp,npgl,leftgca,gca,1,gdna,gdnaR,sia,sib,numbp,pd,pa,pdR,paR,sprm);
	  if (LTFLAG) fprintf(stdout,"\nPGL %3d (- strand): %d to %d (f %d r %d\tfp %d rp %d\tfs %d rs %d)\n",npgl,rightend,leftend,fwdni,revni,fwdppa,revppa,fwdmol,revmol);
	}
      }
      leftgca = gca;
      leftend  = MIN(gca->gcds[0][0],gca->gcds[gca->exn-1][1]);
      rightend = MAX(gca->gcds[0][0],gca->gcds[gca->exn-1][1]);
      fwdpdpa = revpdpa = 0.0;
      fwdni = revni = 0;
      fwdppa = revppa = 0;
      fwdmol = revmol = 0;
      if (gca->exn > 1) {
	if (gca->gcds[0][0] < gca->gcds[gca->exn-1][1]) {
	  for (i = 0; i < gca->exn - 1; ++i) {
	    fwdpdpa += gca->itrscr[i][0];
	    fwdpdpa += gca->itrscr[i][1];
	    ++fwdni;
	  }
	}
	else {
	  for (i = 0; i < gca->exn - 1; ++i) {
	    revpdpa += gca->itrscr[i][0];
	    revpdpa += gca->itrscr[i][1];
	    ++revni;
	  }
	}
      }
      if (gca->ppa[0] < gca->ppa[1]) {
	if (gca->gcds[0][0] < gca->gcds[gca->exn-1][1]) ++fwdppa;
	else                                            ++revppa;
      }
      if (LTFLAG) fprintf(stdout,"\nDETPGL leftend %d rightend %d fwdni %d revni %d gca %s",leftend,rightend,fwdni,revni,gca->cname);
    }
    else {
      rightend = MAX(rightend,MAX(gca->gcds[0][0],gca->gcds[gca->exn-1][1]));
      if (gca->exn > 1) {
	if (gca->gcds[0][0] < gca->gcds[gca->exn-1][1]) {
	  for (i = 0; i < gca->exn - 1; ++i) {
	    fwdpdpa += gca->itrscr[i][0];
	    fwdpdpa += gca->itrscr[i][1];
	    ++fwdni;
	  }
	}
	else {
	  for (i = 0; i < gca->exn - 1; ++i) {
	    revpdpa += gca->itrscr[i][0];
	    revpdpa += gca->itrscr[i][1];
	    ++revni;
	  }
	}
      }
      if (gca->ppa[0] < gca->ppa[1]) {
	if (gca->gcds[0][0] < gca->gcds[gca->exn-1][1]) ++fwdppa;
	else                                            ++revppa;
      }
      if (LTFLAG) fprintf(stdout,"\nDETPGL leftend %d rightend %d fwdni %d revni %d gca %s",leftend,rightend,fwdni,revni,gca->cname);
    }
    gca = gca->next;
  }
  ++npgl;
  if (fwdni+revni > 0) {
    if (fwdpdpa >= revpdpa) {
      r+= orient_gca(fp,npgl,leftgca,gca,0,gdna,gdnaR,sia,sib,numbp,pd,pa,pdR,paR,sprm);
      if (LTFLAG) fprintf(stdout,"\nPGL %3d (+ strand): %d to %d (f %d r %d\tfp %d rp %d)\n",npgl,leftend,rightend,fwdni,revni,fwdppa,revppa);
    }
    else {
      r+= orient_gca(fp,npgl,leftgca,gca,1,gdna,gdnaR,sia,sib,numbp,pd,pa,pdR,paR,sprm);
      if (LTFLAG) fprintf(stdout,"\nPGL %3d (- strand): %d to %d (f %d r %d\tfp %d rp %d)\n",npgl,rightend,leftend,fwdni,revni,fwdppa,revppa);
    }
  }
  else if (fwdppa+revppa > 0) {
    if (fwdppa >= revppa) {
      r+= orient_gca(fp,npgl,leftgca,gca,0,gdna,gdnaR,sia,sib,numbp,pd,pa,pdR,paR,sprm);
      if (LTFLAG) fprintf(stdout,"\nPGL %3d (+ strand): %d to %d (f %d r %d\tfp %d rp %d)\n",npgl,leftend,rightend,fwdni,revni,fwdppa,revppa);
    }
    else {
      r+= orient_gca(fp,npgl,leftgca,gca,1,gdna,gdnaR,sia,sib,numbp,pd,pa,pdR,paR,sprm);
      if (LTFLAG) fprintf(stdout,"\nPGL %3d (- strand): %d to %d (f %d r %d\tfp %d rp %d)\n",npgl,rightend,leftend,fwdni,revni,fwdppa,revppa);
    }
  }
  else {
    fwdmol = maxorflgth(gdna+leftend-1,rightend-leftend+1);
    revmol = maxorflgth(gdnaR+numbp-rightend,rightend-leftend+1);
    if (fwdmol >= revmol) {
      r+= orient_gca(fp,npgl,leftgca,gca,0,gdna,gdnaR,sia,sib,numbp,pd,pa,pdR,paR,sprm);
      if (LTFLAG) fprintf(stdout,"\nPGL %3d (+ strand): %d to %d (f %d r %d\tfp %d rp %d\tfs %d rs %d)\n",npgl,leftend,rightend,fwdni,revni,fwdppa,revppa,fwdmol,revmol);
    }
    else {
      r+= orient_gca(fp,npgl,leftgca,gca,1,gdna,gdnaR,sia,sib,numbp,pd,pa,pdR,paR,sprm);
      if (LTFLAG) fprintf(stdout,"\nPGL %3d (- strand): %d to %d (f %d r %d\tfp %d rp %d\tfs %d rs %d)\n",npgl,rightend,leftend,fwdni,revni,fwdppa,revppa,fwdmol,revmol);
    }
  }
  if (LTFLAG) fprintf(stdout,"\nend rev_sgl; r=%d\n",r);

return(npgl);

}							/* end rev_sgl() */



int is_leftend_pgl(struct gcalgnmnt *gca, int rightend, int delta, int fwdni, int revni, float fwdpdpa, float revpdpa)
{
int i;
float pdpa = 0.0, MINAVGPDPA = 0.35;

  if (MIN(gca->gcds[0][0],gca->gcds[gca->exn-1][1]) > rightend + delta) {
    return(1);
  }
  else {
    if (gca->exn == 1) return(0);
    if (gca->gcds[0][0] < gca->gcds[gca->exn-1][1]) {		/* ... forward */
      if (revni == 0) return(0);				/* ... no conflict */
      if (fwdpdpa > revpdpa) return(0);				/* ... consistent, forward */
      for (i = 0; i < gca->exn - 1; ++i) {
	pdpa += gca->itrscr[i][0];
	pdpa += gca->itrscr[i][1];
      }
      if (revpdpa/(float)(2*revni) > MINAVGPDPA     &&
	  pdpa/(float)(2*(gca->exn-1)) > MINAVGPDPA   ) return(1);	/* ... conflict */
    }
    else {							/* ... reverse */
      if (fwdni == 0) return (0);				/* ... no conflict */
      if (fwdpdpa < revpdpa) return(0);				/* ... consistent, reverse */
      for (i = 0; i < gca->exn - 1; ++i) {
	pdpa += gca->itrscr[i][0];
	pdpa += gca->itrscr[i][1];
      }
      if (LTFLAG) fprintf(stdout,"\nDETPGL fwdpdpa %f av %f MIN %f gca %f\n",fwdpdpa,fwdpdpa/(float)fwdni,MINAVGPDPA,pdpa/(float)(2*(gca->exn-1)));
      if (fwdpdpa/(float)(2*fwdni) > MINAVGPDPA     &&
	  pdpa/(float)(2*(gca->exn-1)) > MINAVGPDPA   ) return(1);	/* ... conflict */
    }
  }
  return(0);

}							/* end is_leftend_pgl() */



int maxorflgth(char *seq, int slgth)
{
int i, j, aa;
int last_stp[3];
int frame, l_orf, maxorf[3], molgth = 0;

  for (i = 0; i <= 2; ++i) {
    last_stp[i] = -1;
    maxorf[i] = 0;
  }
  for (i = 0; i <= slgth - 3; ++i) {
    if (seq[i] > 3 || seq[i + 1] > 3 || seq[i + 2] > 3)
      aa = 22;
    else
      aa = CODON2NUMBER(i+1);
    frame = (i + 1) % 3;
    if (i <= last_stp[frame])
      continue;
    j = i;
    while (aa != 23) {
      j += 3;
      if (j > slgth - 3)
	break;
      if (seq[j] > 3 || seq[j + 1] > 3 || seq[j + 2] > 3)
	aa = 22;
      else
	aa = CODON2NUMBER(j+1);
    }
    last_stp[frame] = j;
    l_orf = (j - i) / 3;
    if (l_orf > maxorf[frame]) maxorf[frame] = l_orf;
    if (maxorf[frame] > molgth) molgth = maxorf[frame];
  }
return (molgth);

}							/* end maxorflgth() */



int orient_gca(FILE *fp, int npgl, struct gcalgnmnt *lgca, struct gcalgnmnt *rgca, int rflag, char *gdna, char *gdnaR, int sia, int sib, int numbp, float *pd, float *pa, float *pdR, float *paR, struct sprmtr sprm)
{
struct gcalgnmnt *gca = lgca;
int r = 0;

  if (LTFLAG) fprintf(stdout,"\nbeg orient_gca for PGL %d",npgl);
  while (gca != rgca) {
    if ((rflag == 0  &&  gca->gcds[0][0] > gca->gcds[gca->exn-1][1])  ||
	(rflag == 1  &&  gca->gcds[0][0] < gca->gcds[gca->exn-1][1])    ) {
      if (LTFLAG) fprintf(stdout,"\nREVERSE %s (%d)\n",gca->cname,gca->calln);
      if (reverse_gca(fp,gca,gdna,gdnaR,sia,sib,numbp,pd,pa,pdR,paR,sprm) == 0) {
	gca->exn = 0;
      }
      if (LTFLAG && gca->exn>0) fprintf(stdout,"\nNOW%s\n",gca->algnmnt);
      if (quality_adj_gca(gca,sprm) == 0)
	if (LTFLAG) fprintf(stdout,"\nLOW QUALITY: REVERSED GCA %s (%d) WILL BE DISCARDED\n",gca->cname,gca->calln);
      ++r;
    }
    gca = gca->next;
  }
  if (LTFLAG) fprintf(stdout,"\nend orient_gca for PGL %d; r=%d",npgl,r);

return(r);

}							/* end orient_gca() */



int reverse_gca(FILE *fp, struct gcalgnmnt *gca, char *gdna, char *gdnaR, int sia, int sib, int numbp, float *pd, float *pa, float *pdR, float *paR, struct sprmtr sprm)
{
  int DBEST = 7;
  char estfn[257], *cdna=NULL;
  int estlg, revest;
  char sfname[257], gname[257], gnameR[257], cname[257], cnameR[257];
  char *cdnaR;
  struct gcalgnmnt gcaR;

  if (openRawTextFile(DBEST, gca->estdbn, 1) != 0) {
    fprintf(stderr, "File %s cannot be opened.\n", gca->estdbn);
    perror(gca->estdbn);
    exit(-1);
  }
  setPositionRawTextFile(DBEST,gca->offset);
  if ((estlg = getlns(DBEST, 1, 1, sfname, cdna)) == -1) {
    closeRawTextFile(DBEST);
    return (0);
  }
  if ((cdna = (char *) calloc(estlg , sizeof(char))) == NULL) {
    fatal_error("Error: memory allocation failed. Exit.\n");
  }
  estlg = getlns(DBEST, 0, 1, estfn, cdna);
  closeRawTextFile(DBEST);

  strcpy(gname, gca->gname);
  strcpy(gnameR, gname);
  if (gname[strlen(gname) - 1] == '+')
    gnameR[strlen(gnameR) - 1] = '-';
  else
    gnameR[strlen(gnameR) - 1] = '+';

  strcpy(cname, gca->cname);
  strcpy(cnameR, cname);
  if (cname[strlen(cname) - 1] == '+') {
    revest = 1;
    cnameR[strlen(cnameR) - 1] = '-';
  }
  else {
    revest = 0;
    cnameR[strlen(cnameR) - 1] = '+';
  }

  strcpy(gcaR.gname, gnameR);
  strcpy(gcaR.gsgmntn, gca->gsgmntn);
  strcpy(gcaR.cname, cnameR);
  strcpy(gcaR.estdbn, gca->estdbn);
  gcaR.calln = gca->calln;
  gcaR.offset = gca->offset;
  gcaR.gia = numbp - 1 - gca->gib;
  gcaR.gib = numbp - 1 - gca->gia;
  gcaR.clgth = gca->clgth;

  if (gca->gcds[0][0] < gca->gcds[gca->exn-1][1]) {
    if (revest) {
      if ((cdnaR = (char *) calloc(estlg , sizeof(char))) == NULL)
        fatal_error("Error: memory allocation failed. Exit.\n");
      complement_seq(cdna, estlg, cdnaR);
      if (sahmtD(fp, gdna, gdnaR, numbp - 1 - sib, numbp - 1 - sia, numbp, 1, pd, pa, pdR, paR, cdnaR, &gcaR, sprm, gca->calln) == 0) {
	free(cdna);
	free(cdnaR);
	return (0);
      }
      free(cdnaR);
    }
    else {
      if (sahmtD(fp, gdna, gdnaR, numbp - 1 - sib, numbp - 1 - sia, numbp, 1, pd, pa, pdR, paR, cdna, &gcaR, sprm, gca->calln) == 0) {
	free(cdna);
	return (0);
      }
    }
  }
  else {
    if (revest) {
      if ((cdnaR = (char *) calloc(estlg , sizeof(char))) == NULL)
        fatal_error("Error: memory allocation failed. Exit.\n");
      complement_seq(cdna, estlg, cdnaR);
      if (sahmtD(fp, gdna, gdnaR, sia, sib, numbp, 0, pd, pa, pdR, paR, cdnaR, &gcaR, sprm, gca->calln) == 0) {
	free(cdna);
	free(cdnaR);
	return (0);
      }
      free(cdnaR);
    }
    else {
      if (sahmtD(fp, gdna, gdnaR, sia, sib, numbp, 0, pd, pa, pdR, paR, cdna, &gcaR, sprm, gca->calln) == 0) {
	free(cdna);
	return (0);
      }
    }
  }

  gcacpy(gca,&gcaR);
  free(cdna);
  return (1);
}							/* end reverse_gca() */



void gcacpy(struct gcalgnmnt *gca1,struct gcalgnmnt *gca2)
{
int i,j;

  strcpy(gca1->gname,gca2->gname);
  strcpy(gca1->gsgmntn,gca2->gsgmntn);
  strcpy(gca1->cname,gca2->cname);
  strcpy(gca1->estdbn,gca2->estdbn);
  gca1->calln = gca2->calln;
  gca1->offset = gca2->offset;
  gca1->gia = gca2->gia;
  gca1->gib = gca2->gib;
  gca1->clgth = gca2->clgth;
  gca1->exn = gca2->exn;
  gca1->ppa[0] = gca2->ppa[0];
  gca1->ppa[1] = gca2->ppa[1];
  for (i=0;i<gca1->exn;++i) {
    gca1->exnscr[i] = gca2->exnscr[i];
    for (j=0;j<2;++j) {
      gca1->gcds[i][j] = gca2->gcds[i][j];
      gca1->ccds[i][j] = gca2->ccds[i][j];
      gca1->itrscr[i][j] = gca2->itrscr[i][j];
    }
  }
  gca1->score = gca2->score;
  if (gca1->algnmnt!=NULL) {
    if (LTFLAG) fprintf(stdout,"\nDEALLOC_ALN2 %ld # %d",(long)(gca1->algnmnt),++NDEALLOC_ALN);
    free((char *)gca1->algnmnt);
  }
  if ((gca1->algnmnt = (char *) calloc(strlen(gca2->algnmnt)+1, sizeof(char))) == NULL)
    fatal_error("Error: memory allocation failed. Exit.\n");
  if (LTFLAG) fprintf(stdout,"\nALLOC_ALN2 %ld # %d",(long)(gca1->algnmnt),++NALLOC_ALN);
  strcpy(gca1->algnmnt,gca2->algnmnt);
  if (gca2->algnmnt!=NULL) {
    if (LTFLAG) fprintf(stdout,"\nDEALLOC_ALN3 %ld # %d",(long)(gca2->algnmnt),++NDEALLOC_ALN);
    free((char *)gca2->algnmnt);
  }

}							/* end gcacpy() */



int prt_gca_list(FILE *fp, struct gcalgnmnt *gca, char *gsgmntn)
{
  int i, ngca = 0;
  char anchorName[257];

  while (gca != NULL) {
    if (htmlop) {
      sprintf(anchorName,"PGS%d@%s",gca->calln,gsgmntn);
      ADD_NAME(fp,anchorName);
    }
    fprintf(fp, "%s", gca->algnmnt);
    fprintf(fp, "hqPGS_%s_%s\t(", gca->gname, gca->cname);
    for (i = 0; i < gca->exn - 1; ++i) {
      fprintf(fp, "%d  %d,", gca->gcds[i][0], gca->gcds[i][1]);
    }
    fprintf(fp, "%d  %d)\n\n", gca->gcds[gca->exn - 1][0], gca->gcds[gca->exn - 1][1]);
    ++ngca;
    gca = gca->next;
  }

#ifdef HTMLWS
  if (htmlop) {
    fprintf(imageDataFh,"@%s %d ",gsgmntn,ngca);
  }
#endif

  return (ngca);

}							/* end prt_gca_list() */



void det_pgl(FILE *fp, struct gcalgnmnt **gcaheadpp, char *gname, char *gsgmntn, char *gdna, char *gdnaR, int sia, int sib, int numbp, float *pd, float *pa, float *pdR, float *paR, struct sprmtr sprm, int bflag, int revestallowed)
{
  int i;
  struct gcalgnmnt *gca;
  struct pgl *pglheadp = NULL;
  int npgl = 0;
  int ngca = 0;
  int rngca;

  gca= *gcaheadpp;
  while (gca != NULL) {
    ++ngca;
    quality_adj_gca(gca,sprm);
    gca= gca->next;
  }
  consolidate_gca(&(*gcaheadpp), &rngca);
  if (LTFLAG  &&  rngca < ngca) fprintf(stdout,"\nNUMBER OF GCAs DISCARDED BEFORE rev_sgl: %d-%d=%d\n",ngca,rngca,ngca-rngca);

  if (bflag == 1  &&  revestallowed == 1) {
    rev_sgl(fp,*gcaheadpp,gdna,gdnaR,sia,sib,numbp,pd,pa,pdR,paR,sprm);

    ngca= rngca;
    consolidate_gca(&(*gcaheadpp), &rngca);
    if (LTFLAG  &&  rngca < ngca) fprintf(stdout,"\nNUMBER OF GCAs DISCARDED AFTER rev_sgl: %d-%d=%d\n",ngca,rngca,ngca-rngca);
  }

  gca= *gcaheadpp;
  while (gca != NULL) {
    for (i = 0; i < MAXNAGS; ++i)
      gca->ags[i] = 0;
    if (build_pgl_list(fp,&pglheadp, gca, gdna,gdnaR,sia,sib,numbp,pd,pa,pdR,paR,sprm))
      ++npgl;
    if (LTFLAG) {
      fprintf(stdout,"\nCURRENT LIST OF PREDICTED GENE STRUCTURES:");
      prt_consensus_pgl(stdout, pglheadp, gname, gsgmntn, gdna, gdnaR, numbp, npgl, bflag);
    }
    gca = gca->next;
  }
  consolidate_gca(&(*gcaheadpp), &rngca);

  form_ags_in_pgl(fp,pglheadp);

  while (consolidate_pgl(&pglheadp, &npgl));
  ngca = prt_gca_list(fp, *gcaheadpp, gsgmntn);
  if (htmlop)
    fprintf(fp, "<A NAME=\"HEAD-PGL-%s\"></A>\n",gsgmntn);
  fprintf(fp,"\n\nTotal number of EST alignments reported: %d\n",ngca);
  fprintf(fp, "\n\n\n________________________________________\
________________________________________\n\n");
  fflush(fp);
  prt_consensus_pgl(fp, pglheadp, gname, gsgmntn, gdna, gdnaR, numbp, npgl, bflag);
  free_pgl(pglheadp);

}							/* end det_pgl() */



int quality_adj_gca(struct gcalgnmnt *gcap, struct sprmtr sprm)
{
int i, ia, ib, j;

if (LTFLAG) fprintf(stdout,"\nQADJ gca %s %d %d (%d exons)",gcap->cname,gcap->gcds[0][0],gcap->gcds[gcap->exn-1][1],gcap->exn);

if (gcap->exn < 1) return(0);
ia = 0;
ib = gcap->exn;

if (gcap->exn == 1) {
  if (weak_exon(gcap,0,-1,sprm)) ++ia;
  }
if (gcap->exn > 1) {
  if (weak_exon(gcap,0,0,sprm)) ++ia;
  if (weak_exon(gcap,gcap->exn-1,1,sprm)) --ib;
  }
if (ia > 0  ||  ib < gcap->exn) {
  for (i = ia, j= 0; i < ib; ++i, ++j) {
    gcap->gcds[j][0] = gcap->gcds[i][0];
    gcap->gcds[j][1] = gcap->gcds[i][1];
    gcap->exnscr[j] = gcap->exnscr[i];
    if (i < gcap->exn - 1) {
      gcap->itrscr[j][0] = gcap->itrscr[i][0];
      gcap->itrscr[j][1] = gcap->itrscr[i][1];
    }
  }
 }

if (LTFLAG) fprintf(stdout,"\tia= %d, ib= %d, ib-ia= %d\n",ia,ib,ib-ia);
if (gcap->exn == ib-ia) return(1);

gcap->exn = ib-ia;
if (gcap->exn > 0) return(quality_adj_gca(gcap,sprm));
else               return(0);

}							/* end quality_adj_gca() */



int weak_exon(struct gcalgnmnt *gcap, int i, int o, struct sprmtr sprm)
{

if (o == -1) {
  if (gcap->exnscr[i] <= sprm.poor_exon_score) {
    if (abs(gcap->gcds[i][1]-gcap->gcds[i][0])+1 <= sprm.short_exon) return(1);
    return(0);
  }
  else {
    if (abs(gcap->gcds[i][1]-gcap->gcds[i][0])+1 <= sprm.tiny_exon) return(1);
    return(0);
  }
 }
if (o == 0) {
  if (gcap->exnscr[i] <= sprm.poor_exon_score) {
    if (abs(gcap->gcds[i][1]-gcap->gcds[i][0])+1 <= sprm.tiny_exon) return(1);
    if (gcap->itrscr[i][0] <= sprm.poor_donor_score) return(1);
    if (abs(gcap->gcds[i+1][0]-gcap->gcds[i][1])-1 >= sprm.long_intron) return(1);
    if (i < gcap->exn-2  &&  weak_exon(gcap,i+1,0,sprm)) return(1);
    return(0);
  }
  else {
    if (abs(gcap->gcds[i][1]-gcap->gcds[i][0])+1 > sprm.tiny_exon) return(0);
    if (gcap->itrscr[i][0] <= sprm.poor_donor_score) return(1);
    if (abs(gcap->gcds[i+1][0]-gcap->gcds[i][1])-1 >= sprm.long_intron) return(1);
    if (i < gcap->exn-2  &&  weak_exon(gcap,i+1,0,sprm)) return(1);
    return(0);
  }
 }
if (o == 1) {
  if (gcap->exnscr[i] <= sprm.poor_exon_score) {
    if (abs(gcap->gcds[i][1]-gcap->gcds[i][0])+1 <= sprm.tiny_exon) return(1);
    if (gcap->itrscr[i-1][1] <= sprm.poor_acptr_score) return(1);
    if (abs(gcap->gcds[i][0]-gcap->gcds[i-1][1])-1 >= sprm.long_intron) return(1);
    if (i > 1  &&  weak_exon(gcap,i-1,1,sprm)) return(1);
    return(0);
  }
  else {
    if (abs(gcap->gcds[i][1]-gcap->gcds[i][0])+1 > sprm.tiny_exon) return(0);
    if (gcap->itrscr[i-1][1] <= sprm.poor_acptr_score) return(1);
    if (abs(gcap->gcds[i][0]-gcap->gcds[i-1][1])-1 >= sprm.long_intron) return(1);
    if (i > 1  &&  weak_exon(gcap,i-1,1,sprm)) return(1);
    return(0);
  }
 }
 return (0);
}							/* end weak_exon() */



int consolidate_gca(struct gcalgnmnt **gcahpp, int *ngca)
{
  struct gcalgnmnt *ip = *gcahpp, *tmp, *newheadp = NULL;
  int rval = 0;

  *ngca = 0;
  while (ip != NULL) {
    if (ip->exn < 1) {			/* ... discard low quality gcas */
      if (LTFLAG) fprintf(stdout,"\nLOW QUALITY: GCA %s (%d %ld) DISCARDED\n",ip->cname,ip->calln,(long)ip);
      tmp = ip->next;
      free_gca(ip);
      ip = tmp;
      continue;
    }
    tmp = ip->next;
    ip->next = NULL;
    if (insert_gca(&newheadp, ip) == 0) {
      if (LTFLAG) fprintf(stdout,"\nNO-INSERT %s %d %ld (%d exons; %d %d)",ip->cname,ip->calln,(long)ip,ip->exn,ip->gcds[0][0],ip->gcds[ip->exn-1][1]);
      free_gca(ip);
      rval = 1;
    }
    else
      ++(*ngca);
    ip = tmp;
  }
  *gcahpp = newheadp;
  return (rval);

}							/* end consolidate_gca() */



int build_pgl_list(FILE *fp, struct pgl **pglhpp, struct gcalgnmnt *gca, char *gdna, char *gdnaR, int sia, int sib, int numbp, float *pd, float *pa, float *pdR, float *paR, struct sprmtr sprm)
{
  int o, gapl1, gapl2;
  struct pgl *pgl;

  if (*pglhpp == NULL) {
    gca2pgl(pglhpp, NULL, gca);
    return (1);
  }
  else {
    pgl = *pglhpp;
    while (pgl != NULL) {
      o = overlap_pgl_gca(fp,pgl,gca,gdna,gdnaR,sia,sib,numbp,pd,pa,pdR,paR,sprm,JNLGTH,&gapl1);
      if (o == -1) return (0);
      if (o ==  1) {
        if (pgl->next != NULL) {	/* ... always append rightmost PGL */
          o = overlap_pgl_gca(fp,pgl->next,gca,gdna,gdnaR,sia,sib,numbp,pd,pa,pdR,paR,sprm,0,&gapl2);
	  if (o == -1) return (0);
	  if (o == 1  &&  gapl2<gapl1) {
	    pgl = pgl->next;
	    continue;
	  }
	}
	gca2pgl(pglhpp, pgl, gca);
	return (0);
      }
      pgl = pgl->next;
    }
    gca2pgl(pglhpp, NULL, gca);
    return (1);
  }

}							/* end build_pgl_list() */



int overlap_pgl_gca(FILE *fp, struct pgl *pglp, struct gcalgnmnt *gcap, char *gdna, char *gdnaR, int sia, int sib, int numbp, float *pd, float *pa, float *pdR, float *paR, struct sprmtr sprm, int delta, int *gapl)
{
  int i;

  if (LTFLAG) fprintf(stdout,"\nOVERLAP? gca %s %d %d versus pgl %d %d delta %d",gcap->cname,gcap->gcds[0][0],gcap->gcds[gcap->exn - 1][1],pglp->gcrd[0],pglp->gcrd[1],delta);

  if (gcap->gcds[0][0] < gcap->gcds[gcap->exn-1][1]) {
    if (pglp->rflag == 0) {
      if (LTFLAG) fprintf(stdout,"\nOVERLAP11? gca %s %d %d versus pgl %d %d",gcap->cname,gcap->gcds[0][0],gcap->gcds[gcap->exn - 1][1],pglp->gcrd[0],pglp->gcrd[1]);
      if (is_xwab(gcap->gcds[gcap->exn - 1][1],pglp->gcrd[0],pglp->gcrd[1],delta)                         ||
          is_xwab(gcap->gcds[0][0],pglp->gcrd[0],pglp->gcrd[1],delta)                                     ||
	  xy_contains_ab(gcap->gcds[0][0],gcap->gcds[gcap->exn - 1][1],pglp->gcrd[0],pglp->gcrd[1],delta)   )
						/* ... same strand gca overlaps pgl */
	{
	  *gapl = MAX(gcap->gcds[0][0] - pglp->gcrd[1], pglp->gcrd[0] - gcap->gcds[gcap->exn - 1][1]);
	  return(1);
	}
      else
	return(0);
    }
    else {
      if (gcap->exn > 1) return(0);
      for (i=0;i<pglp->nags;++i) {
	if (pglp->exn[i] > 1) break;
      }
      if (i==pglp->nags) return (0);
						/* ... singlet versus opposite strand multi-exon pgl case */
      if (LTFLAG) fprintf(stdout,"\nOVERLAP12? gca %s %d %d versus pgl %d %d",gcap->cname,gcap->gcds[0][0],gcap->gcds[0][1],pglp->gcrd[0],pglp->gcrd[1]);
      if (is_xwab(gcap->gcds[gcap->exn - 1][1],pglp->gcrd[1],pglp->gcrd[0],delta)                         ||
          is_xwab(gcap->gcds[0][0],pglp->gcrd[1],pglp->gcrd[0],delta)                                     ||
	  xy_contains_ab(gcap->gcds[0][0],gcap->gcds[gcap->exn - 1][1],pglp->gcrd[1],pglp->gcrd[0],delta)   ) {
	if (LTFLAG) fprintf(stdout,"\nREVERSE1 GCA %s (%d) relative to PGS %s AGS %d (%d exons)\n",gcap->cname,gcap->calln,(pglp->gca)->cname,i,pglp->exn[i]);
	if (reverse_gca(fp,gcap,gdna,gdnaR,sia,sib,numbp,pd,pa,pdR,paR,sprm) == 0) {
	  gcap->exn = 0;
	}
	if (LTFLAG && gcap->exn>0) fprintf(stdout,"\nNOW%s\n",gcap->algnmnt);
	if (quality_adj_gca(gcap,sprm) == 0) {
	  if (LTFLAG) fprintf(stdout,"\nLOW QUALITY: GCA %s (%d) DISCARDED\n",gcap->cname,gcap->calln);
	  return (-1);
	}
	*gapl = MAX(gcap->gcds[gcap->exn - 1][1] - pglp->gcrd[0], pglp->gcrd[1] - gcap->gcds[0][0]);
	return (1);
      }
      else
	return (0);
    }
  }
  else {
    if (pglp->rflag == 1) {
      if (LTFLAG) fprintf(stdout,"\nOVERLAP21? gca %s %d %d versus pgl %d %d",gcap->cname,gcap->gcds[0][0],gcap->gcds[gcap->exn - 1][1],pglp->gcrd[0],pglp->gcrd[1]);
      if (is_xwab(gcap->gcds[0][0],pglp->gcrd[1],pglp->gcrd[0],delta)                                     ||
          is_xwab(gcap->gcds[gcap->exn - 1][1],pglp->gcrd[1],pglp->gcrd[0],delta)                         ||
	  xy_contains_ab(gcap->gcds[gcap->exn - 1][0],gcap->gcds[0][0],pglp->gcrd[1],pglp->gcrd[0],delta)   )
						/* ... same strand gca overlaps pgl */
	{
	  *gapl = MAX(gcap->gcds[gcap->exn - 1][1] - pglp->gcrd[0], pglp->gcrd[1] - gcap->gcds[0][0]);
	  return(1);
	}
      else
	return(0);
    }
    else {
      if (gcap->exn > 1) return(0);
      for (i=0;i<pglp->nags;++i) {
	if (pglp->exn[i] > 1) break;
      }
      if (i==pglp->nags) return (0);
						/* ... singlet versus opposite strand multi-exon pgl case */
      if (LTFLAG) fprintf(stdout,"\nOVERLAP22? gca %s %d %d versus pgl %d %d",gcap->cname,gcap->gcds[0][0],gcap->gcds[gcap->exn - 1][1],pglp->gcrd[0],pglp->gcrd[1]);
      if (is_xwab(gcap->gcds[gcap->exn - 1][1],pglp->gcrd[0],pglp->gcrd[1],delta)                         ||
          is_xwab(gcap->gcds[0][0],pglp->gcrd[0],pglp->gcrd[1],delta)                                     ||
	  xy_contains_ab(gcap->gcds[gcap->exn - 1][0],gcap->gcds[0][0],pglp->gcrd[0],pglp->gcrd[1],delta)   ) {
	if (LTFLAG) fprintf(stdout,"\nREVERSE2 GCA %s (%d) relative to PGS %s AGS %d (%d exons)\n",gcap->cname,gcap->calln,(pglp->gca)->cname,i,pglp->exn[i]);
	if (reverse_gca(fp,gcap,gdna,gdnaR,sia,sib,numbp,pd,pa,pdR,paR,sprm) == 0) {
	  gcap->exn = 0;
	}
	if (LTFLAG && gcap->exn>0) fprintf(stdout,"\nNOW%s\n",gcap->algnmnt);
	if (quality_adj_gca(gcap,sprm) == 0) {
	  if (LTFLAG) fprintf(stdout,"\nLOW QUALITY: GCA %s (%d) DISCARDED\n",gcap->cname,gcap->calln);
	  return (-1);
	}
	*gapl = MAX(gcap->gcds[0][0] - pglp->gcrd[1], pglp->gcrd[0] - gcap->gcds[gcap->exn - 1][1]);
	return (1);
      }
      else
	return (0);
    }
  }

}							/* end overlap_pgl_gca() */



int is_xwab(int x,int a,int b,int delta)
{

if (a-delta <= x  &&  x <= b+delta)  return(1);
else  return(0);

}							/* end is_xwab() */



int xy_contains_ab(int x,int y,int a,int b,int delta)
{

if (x <= a-delta  &&  y >= b+delta)  return(1);
else  return(0);

}							/* end xy_contains_ab() */


int gca2pgl(pglhpp, pglp, gcap)
  struct pgl **pglhpp, *pglp;
  struct gcalgnmnt *gcap;
{
  struct pgl *pgl;
  struct gcalgnmnt *gca;

  if (pglp == NULL) {
    if ((pgl = (struct pgl *) calloc(1 , sizeof(struct pgl))) == NULL)
       fatal_error("Error: memory allocation failed. Exit.\n");
    if (LTFLAG) fprintf(stdout,"\nALLOC_PGL %ld # %d",(long)pgl,++NALLOC_PGL);
    gcap->link = NULL;
    if (LTFLAG) fprintf(stdout,"\nCREATING new pgl from gca %s %d %d\n",gcap->cname,gcap->gcds[0][0],gcap->gcds[gcap->exn-1][1]);
    pgl->gca = gcap;
    pgl->gcrd[0] = gcap->gcds[0][0];
    pgl->gcrd[1] = gcap->gcds[gcap->exn - 1][1];
    if (pgl->gcrd[0] < pgl->gcrd[1])
      pgl->rflag = 0;
    else
      pgl->rflag = 1;
    pgl->nags = 0;
    if (insert_pgl(pglhpp, pgl) == 0) {
      if (LTFLAG) fprintf(stdout,"\nDEALLOC_PGL %ld # %d NO-INSERT IN gca2pgl",(long)pgl,++NDEALLOC_PGL);
      free(pgl);
    }
  }
  else {
    gca = pglp->gca;
    while (gca->link != NULL) {
      if (compare_gcas(gca,gcap) <= 1) {
	gcap->exn = 0;
	return (0);
      }
      gca = gca->link;
    }
    if (compare_gcas(gca,gcap) <= 1) {
	gcap->exn = 0;
	return (0);
    }
    gca->link = gcap;
    gcap->link = NULL;
    if (pglp->rflag == 0) {
      if (gcap->gcds[0][0] < pglp->gcrd[0])
	pglp->gcrd[0] = gcap->gcds[0][0];
      if (gcap->gcds[gcap->exn - 1][1] > pglp->gcrd[1])
	pglp->gcrd[1] = gcap->gcds[gcap->exn - 1][1];
    }
    else {
      if (gcap->gcds[0][0] > pglp->gcrd[0])
	pglp->gcrd[0] = gcap->gcds[0][0];
      if (gcap->gcds[gcap->exn - 1][1] < pglp->gcrd[1])
	pglp->gcrd[1] = gcap->gcds[gcap->exn - 1][1];
    }
  }
  return (1);
}							/* end gca2pgl() */



int insert_pgl(pglhpp, pgl)
  struct pgl **pglhpp, *pgl;
{
  int iflag;

  if (*pglhpp == NULL) {
    *pglhpp = pgl;
    pgl->next = NULL;
  }
  else {
    if (overlap_pgl_pgl(*pglhpp, pgl, JNLGTH)) {
      merge_pgl_pgl(*pglhpp, pgl);
      return (0);
    }
    iflag = 1;
    if ((*pglhpp)->rflag == 0) {
      if (pgl->rflag == 0) {
	if (pgl->gcrd[0] > (*pglhpp)->gcrd[1])
	  iflag = 0;
      }
      else {
	if (pgl->gcrd[1] > (*pglhpp)->gcrd[0])
	  iflag = 0;
      }
    }
    else {
      if (pgl->rflag == 0) {
	if (pgl->gcrd[0] > (*pglhpp)->gcrd[1])
	  iflag = 0;
      }
      else {
	if (pgl->gcrd[1] > (*pglhpp)->gcrd[0])
	  iflag = 0;
      }
    }
    while ((*pglhpp)->next != NULL && iflag == 0) {
      pglhpp = &((*pglhpp)->next);
      if (overlap_pgl_pgl(*pglhpp, pgl, JNLGTH)) {
	merge_pgl_pgl(*pglhpp, pgl);
	return (0);
      }
      iflag = 1;
      if ((*pglhpp)->rflag == 0) {
	if (pgl->rflag == 0) {
	  if (pgl->gcrd[0] > (*pglhpp)->gcrd[1])
	    iflag = 0;
	}
	else {
	  if (pgl->gcrd[1] > (*pglhpp)->gcrd[0])
	    iflag = 0;
	}
      }
      else {
	if (pgl->rflag == 0) {
	  if (pgl->gcrd[0] > (*pglhpp)->gcrd[1])
	    iflag = 0;
	}
	else {
	  if (pgl->gcrd[1] > (*pglhpp)->gcrd[0])
	    iflag = 0;
	}
      }
    }
    if ((*pglhpp)->next == NULL) {
      if (iflag == 0) {
	(*pglhpp)->next = pgl;
	pgl->next = NULL;
      }
      else {
	pgl->next = *pglhpp;
	*pglhpp = pgl;
      }
    }
    else {
      pgl->next = *pglhpp;
      *pglhpp = pgl;
    }
  }
  return (1);

}							/* end insert_pgl() */



int overlap_pgl_pgl(struct pgl *pgl1, struct pgl *pgl2, int delta)
{

  if (pgl1->rflag == 0) {
    if (pgl2->rflag == 1)
      return (0);
    if (is_xwab(pgl2->gcrd[1],pgl1->gcrd[0],pgl1->gcrd[1],delta)                      ||
        is_xwab(pgl2->gcrd[0],pgl1->gcrd[0],pgl1->gcrd[1],delta)                      ||
        xy_contains_ab(pgl2->gcrd[0],pgl2->gcrd[1],pgl1->gcrd[0],pgl1->gcrd[1],delta)   )
      return (1);
    else
      return (0);
  }
  else {
    if (pgl2->rflag == 0)
      return (0);
    if (is_xwab(pgl2->gcrd[0],pgl1->gcrd[1],pgl1->gcrd[0],delta)                      ||
        is_xwab(pgl2->gcrd[1],pgl1->gcrd[1],pgl1->gcrd[0],delta)                      ||
        xy_contains_ab(pgl2->gcrd[1],pgl2->gcrd[0],pgl1->gcrd[1],pgl1->gcrd[0],delta)   )
      return (1);
    else
      return (0);
  }

}							/* end overlap_pgl_pgl() */



int merge_pgl_pgl(pgl1, pgl2)
  struct pgl *pgl1, *pgl2;
{
  int i;
  struct gcalgnmnt *gca, *gcah;

  if (pgl1->rflag == 0) {
    if (pgl2->gcrd[0] < pgl1->gcrd[0])
      pgl1->gcrd[0] = pgl2->gcrd[0];
    if (pgl2->gcrd[1] > pgl1->gcrd[1])
      pgl1->gcrd[1] = pgl2->gcrd[1];
  }
  else {
    if (pgl2->gcrd[0] > pgl1->gcrd[0])
      pgl1->gcrd[0] = pgl2->gcrd[0];
    if (pgl2->gcrd[1] < pgl1->gcrd[1])
      pgl1->gcrd[1] = pgl2->gcrd[1];
  }

  gca = pgl2->gca;
  gcah = pgl1->gca;
  while (gcah->link != NULL)
    gcah = gcah->link;
  gcah->link = gca;
  while (gca != NULL) {
    for (i = 0; i < MAXNAGS; ++i)
      gca->ags[i] = 0;
    for (i = 0; i < pgl1->nags; ++i) {
      if (merge_pgl_gca(pgl1, i, gca) == 1)
	break;
    }
    if (i == pgl1->nags) {
      if (i < MAXNAGS - 1)
        new_ags(pgl1, i, gca);
      else {
	fprintf(stdout,"\nPROBLEM:  TOO MANY AGS\t");
	fprintf(stdout,"GCA %s (%d) NOT COPIED INTO PGL %d %d.\n", gca->cname, gca->calln, pgl1->gcrd[0], pgl1->gcrd[1]);
      }
    }
    gca = gca->link;
  }
  return 0;
}							/* end merge_pgl_pgl() */



int form_ags_in_pgl(FILE *fp, struct pgl *pglhp)
{
  int i;
  struct pgl *pglp = pglhp;
  struct gcalgnmnt *gcap;

  while (pglp != NULL) {
    gcap = pglp->gca;
    sort_gca_in_pgl(&gcap);
    pglp->gca = gcap;
    while (gcap != NULL) {
      if (pglp->nags == 0) {
	new_ags(pglp, 0, gcap);
	gcap = gcap->link;
	continue;
      }
      for (i = 0; i < pglp->nags; ++i) {
        if (LTFLAG) fprintf(stdout,"\nCHECKING pgl %d %d ags %d versus gca %s %d %d\n",pglp->gcrd[0],pglp->gcrd[1],i,gcap->cname,gcap->gcds[0][0],gcap->gcds[gcap->exn-1][1]);
        if (merge_pgl_gca(pglp, i, gcap) == 1)
	  break;
      }
      if (i >= MAXNAGS) {
        fprintf(stdout,"\nPROBLEM:  TOO MANY AGS\t");
        fprintf(stdout,"GCA %s (%d) NOT MERGED INTO PGL %d %d.\n", gcap->cname, gcap->calln, pglp->gcrd[0], pglp->gcrd[1]);
      }
      else if (i == pglp->nags) {
        new_ags(pglp, i, gcap);
      }
      gcap = gcap->link;
    }
  pglp = pglp->next;
  }
  return (1);
}							/* end form_ags_in_pgl() */



int sort_gca_in_pgl(struct gcalgnmnt **gcahpp)
{
  struct gcalgnmnt *ip = *gcahpp, *tmp, *newheadp = NULL;
  int rval = 0;

  while (ip != NULL) {
    tmp = ip->link;
    ip->link = NULL;
    if (insert_gca_by_link(&newheadp, ip) == 0) {
      free(ip);
      rval = 1;
    }
    ip = tmp;
  }
  *gcahpp = newheadp;
  return (rval);

}							/* end sort_gca_in_pgl() */



int insert_gca_by_link(struct gcalgnmnt **gcahpp, struct gcalgnmnt *nodep)
{
int c;

  if (*gcahpp == NULL) {
	*gcahpp = nodep;
	nodep->link = NULL;
  }
  else {
	while ((*gcahpp)->link != NULL) {
		c = order_gcas(*gcahpp,nodep);
		if (c == 1)
		  gcahpp = &((*gcahpp)->link);
		else
		  break;
        }
        if ((*gcahpp)->link == NULL) {
        	c = order_gcas(*gcahpp,nodep);
		if (c == 1) {(*gcahpp)->link = nodep; nodep->link = NULL;}
		else
		  {nodep->link = *gcahpp; *gcahpp = nodep;}
        }
	else {
		nodep->link = *gcahpp; *gcahpp = nodep;
	}
  }
  return (1);

}							/* end insert_gca_by_link() */



int order_gcas(gca1p, gca2p)
  struct gcalgnmnt *gca1p, *gca2p;
{
	/* return 1: gca1p "<" gca2p         */
	/* return 2: gca1p ">" gca2p         */

int beg1 = MIN(gca1p->gcds[0][0],gca1p->gcds[gca1p->exn - 1][1]);
int end1 = MAX(gca1p->gcds[0][0],gca1p->gcds[gca1p->exn - 1][1]);
int beg2 = MIN(gca2p->gcds[0][0],gca2p->gcds[gca2p->exn - 1][1]);
int end2 = MAX(gca2p->gcds[0][0],gca2p->gcds[gca2p->exn - 1][1]);

  if ( (beg1 < beg2) || (beg1 == beg2  &&  end1 >= end2) )  return(1);
  return(2);

}							/* end order_gcas() */



int merge_pgl_gca(pgl, a, gca)
  struct pgl *pgl;
  int a;
  struct gcalgnmnt *gca;
{
  int i, j = 0, ovrlp = 0;

  for (i = 0; i < pgl->exn[a]; ++i) {
    for (j = 0; j < gca->exn; ++j) {
      if (pgl->rflag == 0) {
	if ((pgl->gcds[a][i][0] >= gca->gcds[j][0] &&
	     pgl->gcds[a][i][0] <= gca->gcds[j][1]   ) ||
	    (gca->gcds[j][0] >= pgl->gcds[a][i][0] && gca->gcds[j][0] <= pgl->gcds[a][i][1])) {
	  ovrlp = 1;
	  break;
	}
      }
      else {
	if ((pgl->gcds[a][i][0] <= gca->gcds[j][0] &&
	     pgl->gcds[a][i][0] >= gca->gcds[j][1]   ) ||
	    (gca->gcds[j][0] <= pgl->gcds[a][i][0] && gca->gcds[j][0] >= pgl->gcds[a][i][1])) {
	  ovrlp = 1;
	  break;
	}
      }
    }
    if (ovrlp)
      break;
  }

  if (!ovrlp) {
    if (LTFLAG)
      fprintf(stdout,"NO OVERLAP AGS-%d %d %d gca %d %d\n",
	a + 1, pgl->gcds[a][0][0], pgl->gcds[a][pgl->exn[a] - 1][1], gca->gcds[0][0], gca->gcds[gca->exn - 1][1]);
    return (0);
  }
  else {
    if (i == 0) {
      if (LTFLAG)
	fprintf(stdout,"OVERLAP AGS-%d exon %d (%d %d) with GCA exon %d (%d %d)\n",
	  a + 1, i + 1, pgl->gcds[a][i][0], pgl->gcds[a][i][1], j + 1, gca->gcds[j][0], gca->gcds[j][1]);
      return (check_pgl_gca(pgl, a, 0, gca, j));
    }
    else if (j == 0) {
      if (LTFLAG)
	fprintf(stdout,"OVERLAP AGS-%d exon %d (%d %d) with GCA exon %d (%d %d)\n",
	  a + 1, i + 1, pgl->gcds[a][i][0], pgl->gcds[a][i][1], j + 1, gca->gcds[j][0], gca->gcds[j][1]);
      return (check_pgl_gca(pgl, a, i, gca, 0));
    }
    else {						/* i>0 && j>0 */
      if (LTFLAG)
	fprintf(stdout,"INCONSISTENT: AGS-%d exon %d (%d %d) with GCA exon %d (%d %d)\n",
	  a + 1, i + 1, pgl->gcds[a][i][0], pgl->gcds[a][i][1], j + 1, gca->gcds[j][0], gca->gcds[j][1]);
      return (0);
    }
  }

}							/* end merge_pgl_gca() */



int check_pgl_gca(pgl, a, p, gca, g)
  struct pgl *pgl;
  int a, p, g;
  struct gcalgnmnt *gca;
{
  int i, j;


  for (i = p + 1, j = g + 1; i < pgl->exn[a] - 1 && j < gca->exn - 1; ++i, ++j) {
    if (pgl->gcds[a][i][0] != gca->gcds[j][0] || pgl->gcds[a][i][1] != gca->gcds[j][1]) {
      if (LTFLAG)
	fprintf(stdout,"internal mm AGS-%d exon-%d %d %d vs gca exon-%d %d %d\n",
	  a + 1, i, pgl->gcds[a][i][0], pgl->gcds[a][i][1], j, gca->gcds[j][0], gca->gcds[j][1]);
      return (0);
    }
  }
  if (LTFLAG)
    fprintf(stdout,"pgl last (%d) = %d %d; gca last (%d) = %d %d\n",
      pgl->exn[a],
      pgl->gcds[a][pgl->exn[a] - 1][0], pgl->gcds[a][pgl->exn[a] - 1][1],
      gca->exn, gca->gcds[gca->exn - 1][0], gca->gcds[gca->exn - 1][1]);
  if (LTFLAG)
    fprintf(stdout,"beg to check: pgl-%d %d %d vs gca-%d %d %d\n",
      p + 1, pgl->gcds[a][p][0], pgl->gcds[a][p][1], g + 1, gca->gcds[g][0], gca->gcds[g][1]);
  if (i < pgl->exn[a] && j < gca->exn) {
    if (LTFLAG)
      fprintf(stdout,"end to check: pgl-%d %d %d vs gca-%d %d %d\n",
	i + 1, pgl->gcds[a][i][0], pgl->gcds[a][i][1], j + 1, gca->gcds[j][0], gca->gcds[j][1]);
    if (is_consistent(pgl, a, p, gca, g) && is_consistent(pgl, a, i, gca, j)) {
      adj_ags(pgl, a, p, gca, g);
      return (1);
    }
    else
      return (0);
  }
  else {
    if (is_consistent(pgl, a, p, gca, g)) {
      adj_ags(pgl, a, p, gca, g);
      return (1);
    }
    else
      return (0);
  }

}							/* end check_pgl_gca() */



int is_consistent(pgl, a, p, gca, g)
  struct pgl *pgl;
  int a, p, g;
  struct gcalgnmnt *gca;
{
int TOLERANCE = 25;

  if (p == 0) {
    if (p == pgl->exn[a] - 1) {
      if (g == 0) {
	if (g == gca->exn - 1) {
	  if (LTFLAG)
	    fprintf(stdout,"case  1: %d (%d %d) vs %d (%d %d)\n",
	      p + 1, pgl->gcds[a][p][0], pgl->gcds[a][p][1], g + 1, gca->gcds[g][0], gca->gcds[g][1]);
	}
	else {
	  if (LTFLAG)
	    fprintf(stdout,"case  2: %d (%d %d) vs %d (%d %d)\n",
	      p + 1, pgl->gcds[a][p][0], pgl->gcds[a][p][1], g + 1, gca->gcds[g][0], gca->gcds[g][1]);
	  if (pgl->rflag == 0) {
	    if (gca->gcds[g][1] < pgl->gcds[a][p][1] - TOLERANCE)
	      return (0);
	  }
	  else {
	    if (gca->gcds[g][1] > pgl->gcds[a][p][1] + TOLERANCE)
	      return (0);
	  }
	}
      }
      else {
	if (g == gca->exn - 1) {
	  if (LTFLAG)
	    fprintf(stdout,"case  3: %d (%d %d) vs %d (%d %d)\n",
	      p + 1, pgl->gcds[a][p][0], pgl->gcds[a][p][1], g + 1, gca->gcds[g][0], gca->gcds[g][1]);
	  if (pgl->rflag == 0) {
	    if (gca->gcds[g][0] > pgl->gcds[a][p][0] + TOLERANCE)
	      return (0);
	  }
	  else {
	    if (gca->gcds[g][0] < pgl->gcds[a][p][0] - TOLERANCE)
	      return (0);
	  }
	}
	else {
	  if (LTFLAG)
	    fprintf(stdout,"case  4: %d (%d %d) vs %d (%d %d)\n",
	      p + 1, pgl->gcds[a][p][0], pgl->gcds[a][p][1], g + 1, gca->gcds[g][0], gca->gcds[g][1]);
	  if (pgl->rflag == 0) {
	    if (gca->gcds[g][0] > pgl->gcds[a][p][0] + TOLERANCE || gca->gcds[g][1] < pgl->gcds[a][p][1] - TOLERANCE)
	      return (0);
	  }
	  else {
	    if (gca->gcds[g][0] < pgl->gcds[a][p][0] - TOLERANCE || gca->gcds[g][1] > pgl->gcds[a][p][1] + TOLERANCE)
	      return (0);
	  }
	}
      }
    }
    else {
      if (g == 0) {
	if (g == gca->exn - 1) {
	  if (LTFLAG)
	    fprintf(stdout,"case  5: %d (%d %d) vs %d (%d %d)\n",
	      p + 1, pgl->gcds[a][p][0], pgl->gcds[a][p][1], g + 1, gca->gcds[g][0], gca->gcds[g][1]);
	  if (pgl->rflag == 0) {
	    if (gca->gcds[g][1] > pgl->gcds[a][p][1] + TOLERANCE)
	      return (0);
	  }
	  else {
	    if (gca->gcds[g][1] < pgl->gcds[a][p][1] - TOLERANCE)
	      return (0);
	  }
	}
	else {
	  if (LTFLAG)
	    fprintf(stdout,"case  6: %d (%d %d) vs %d (%d %d)\n",
	      p + 1, pgl->gcds[a][p][0], pgl->gcds[a][p][1], g + 1, gca->gcds[g][0], gca->gcds[g][1]);
	  if (pgl->rflag == 0) {
	    if (gca->gcds[g][1] != pgl->gcds[a][p][1])
	      return (0);
	  }
	  else {
	    if (gca->gcds[g][1] != pgl->gcds[a][p][1])
	      return (0);
	  }
	}
      }
      else {
	if (g == gca->exn - 1) {
	  if (LTFLAG)
	    fprintf(stdout,"case  7: %d (%d %d) vs %d (%d %d)\n",
	      p + 1, pgl->gcds[a][p][0], pgl->gcds[a][p][1], g + 1, gca->gcds[g][0], gca->gcds[g][1]);
	  if (pgl->rflag == 0) {
	    if (gca->gcds[g][0] > pgl->gcds[a][p][0] + TOLERANCE || gca->gcds[g][1] > pgl->gcds[a][p][1] + TOLERANCE)
	      return (0);
	  }
	  else {
	    if (gca->gcds[g][0] < pgl->gcds[a][p][0] - TOLERANCE || gca->gcds[g][1] < pgl->gcds[a][p][1] - TOLERANCE)
	      return (0);
	  }
	}
	else {
	  if (LTFLAG)
	    fprintf(stdout,"case  8: %d (%d %d) vs %d (%d %d)\n",
	      p + 1, pgl->gcds[a][p][0], pgl->gcds[a][p][1], g + 1, gca->gcds[g][0], gca->gcds[g][1]);
	  if (pgl->rflag == 0) {
	    if (gca->gcds[g][0] > pgl->gcds[a][p][0] + TOLERANCE || gca->gcds[g][1] != pgl->gcds[a][p][1])
	      return (0);
	  }
	  else {
	    if (gca->gcds[g][0] < pgl->gcds[a][p][0] - TOLERANCE || gca->gcds[g][1] != pgl->gcds[a][p][1])
	      return (0);
	  }
	}
      }
    }
  }
  else {
    if (p == pgl->exn[a] - 1) {
      if (g == 0) {
	if (g == gca->exn - 1) {
	  if (LTFLAG)
	    fprintf(stdout,"case  9: %d (%d %d) vs %d (%d %d)\n",
	      p + 1, pgl->gcds[a][p][0], pgl->gcds[a][p][1], g + 1, gca->gcds[g][0], gca->gcds[g][1]);
	  if (pgl->rflag == 0) {
	    if (gca->gcds[g][0] < pgl->gcds[a][p][0] - TOLERANCE)
	      return (0);
	  }
	  else {
	    if (gca->gcds[g][0] > pgl->gcds[a][p][0] + TOLERANCE)
	      return (0);
	  }
	}
	else {
	  if (LTFLAG)
	    fprintf(stdout,"case 10: %d (%d %d) vs %d (%d %d)\n",
	      p + 1, pgl->gcds[a][p][0], pgl->gcds[a][p][1], g + 1, gca->gcds[g][0], gca->gcds[g][1]);
	  if (pgl->rflag == 0) {
	    if (gca->gcds[g][0] < pgl->gcds[a][p][0] - TOLERANCE || gca->gcds[g][1] < pgl->gcds[a][p][1] - TOLERANCE)
	      return (0);
	  }
	  else {
	    if (gca->gcds[g][0] > pgl->gcds[a][p][0] + TOLERANCE || gca->gcds[g][1] > pgl->gcds[a][p][1] + TOLERANCE)
	      return (0);
	  }
	}
      }
      else {
	if (g == gca->exn - 1) {
	  if (LTFLAG)
	    fprintf(stdout,"case 11: %d (%d %d) vs %d (%d %d)\n",
	      p + 1, pgl->gcds[a][p][0], pgl->gcds[a][p][1], g + 1, gca->gcds[g][0], gca->gcds[g][1]);
	  if (pgl->rflag == 0) {
	    if (gca->gcds[g][0] != pgl->gcds[a][p][0])
	      return (0);
	  }
	  else {
	    if (gca->gcds[g][0] != pgl->gcds[a][p][0])
	      return (0);
	  }
	}
	else {
	  if (LTFLAG)
	    fprintf(stdout,"case 12: %d (%d %d) vs %d (%d %d)\n",
	      p + 1, pgl->gcds[a][p][0], pgl->gcds[a][p][1], g + 1, gca->gcds[g][0], gca->gcds[g][1]);
	  if (pgl->rflag == 0) {
	    if (gca->gcds[g][0] != pgl->gcds[a][p][0] || gca->gcds[g][1] < pgl->gcds[a][p][1])
	      return (0);
	  }
	  else {
	    if (gca->gcds[g][0] != pgl->gcds[a][p][0] || gca->gcds[g][1] > pgl->gcds[a][p][1])
	      return (0);
	  }
	}
      }
    }
    else {
      if (g == 0) {
	if (g == gca->exn - 1) {
	  if (LTFLAG)
	    fprintf(stdout,"case 13: %d (%d %d) vs %d (%d %d)\n",
	      p + 1, pgl->gcds[a][p][0], pgl->gcds[a][p][1], g + 1, gca->gcds[g][0], gca->gcds[g][1]);
	  if (pgl->rflag == 0) {
	    if (gca->gcds[g][0] < pgl->gcds[a][p][0] - TOLERANCE || gca->gcds[g][1] > pgl->gcds[a][p][1] + TOLERANCE)
	      return (0);
	  }
	  else {
	    if (gca->gcds[g][0] > pgl->gcds[a][p][0] + TOLERANCE || gca->gcds[g][1] < pgl->gcds[a][p][1] - TOLERANCE)
	      return (0);
	  }
	}
	else {
	  if (LTFLAG)
	    fprintf(stdout,"case 14: %d (%d %d) vs %d (%d %d)\n",
	      p + 1, pgl->gcds[a][p][0], pgl->gcds[a][p][1], g + 1, gca->gcds[g][0], gca->gcds[g][1]);
	  if (pgl->rflag == 0) {
	    if (gca->gcds[g][0] < pgl->gcds[a][p][0] - TOLERANCE || gca->gcds[g][1] != pgl->gcds[a][p][1])
	      return (0);
	  }
	  else {
	    if (gca->gcds[g][0] > pgl->gcds[a][p][0] + TOLERANCE || gca->gcds[g][1] != pgl->gcds[a][p][1])
	      return (0);
	  }
	}
      }
      else {
	if (g == gca->exn - 1) {
	  if (LTFLAG)
	    fprintf(stdout,"case 15: %d (%d %d) vs %d (%d %d)\n",
	      p + 1, pgl->gcds[a][p][0], pgl->gcds[a][p][1], g + 1, gca->gcds[g][0], gca->gcds[g][1]);
	  if (pgl->rflag == 0) {
	    if (gca->gcds[g][0] != pgl->gcds[a][p][0] || gca->gcds[g][1] > pgl->gcds[a][p][1] + TOLERANCE)
	      return (0);
	  }
	  else {
	    if (gca->gcds[g][0] != pgl->gcds[a][p][0] || gca->gcds[g][1] < pgl->gcds[a][p][1] - TOLERANCE)
	      return (0);
	  }
	}
	else {
	  if (LTFLAG)
	    fprintf(stdout,"case 16: %d (%d %d) vs %d (%d %d)\n",
	      p + 1, pgl->gcds[a][p][0], pgl->gcds[a][p][1], g + 1, gca->gcds[g][0], gca->gcds[g][1]);
	  if (gca->gcds[g][0] != pgl->gcds[a][p][0] || gca->gcds[g][1] != pgl->gcds[a][p][1])
	    return (0);
	}
      }
    }
  }

  return (1);

}							/* end is_consistent() */



int adj_ags(pgl, a, p, gca, g)
  struct pgl *pgl;
  int a, p, g;
  struct gcalgnmnt *gca;
{
  int i, j, ja, k, n;
  int cds[MAXNEXNS][2], cdsl[MAXNEXNS];
  float exnscr[MAXNEXNS], itrscr[MAXNEXNS][2];

  if (LTFLAG) {
    fprintf(stdout,"MERGE: AGS-%d %d %d and GCA %d %d\n",
      a + 1, pgl->gcds[a][0][0], pgl->gcds[a][pgl->exn[a] - 1][1], gca->gcds[0][0], gca->gcds[gca->exn - 1][1]);

    fprintf(stdout,"\n  AGS (");
    for (i = 0; i < pgl->exn[a] - 1; ++i)
      fprintf(stdout,"%d  %d,", pgl->gcds[a][i][0], pgl->gcds[a][i][1]);
    fprintf(stdout,"%d  %d)", pgl->gcds[a][i][0], pgl->gcds[a][i][1]);
    fprintf(stdout,"\n  GCA (");
    for (i = 0; i < gca->exn - 1; ++i)
      fprintf(stdout,"%d  %d,", gca->gcds[i][0], gca->gcds[i][1]);
    fprintf(stdout,"%d  %d)", gca->gcds[i][0], gca->gcds[i][1]);
  }

  n = 0;
  for (i = p, j = g; i > 0 && j > 0; --i, --j);
  if (LTFLAG)
    fprintf(stdout,"\nIPJG i=%d p=%d j=%d g=%d", i, p, j, g);
  if (i > 0) {
    for (k = 0; k < p; ++k) {			/* ... gca left end contained in pgl left end */
      cds[k][0] = pgl->gcds[a][k][0];
      cds[k][1] = pgl->gcds[a][k][1];
      cdsl[k] = abs(cds[k][1]-cds[k][0]);
      ++n;
    }
  }
  else if (j > 0) {				/* ... pgl left end contained in gca left end */
    for (k = 0; k < g; ++k) {
      cds[k][0] = gca->gcds[k][0];
      cds[k][1] = gca->gcds[k][1];
      cdsl[k] = abs(cds[k][1]-cds[k][0]);
      ++n;
    }
  }
  else if (p > 0) {					/* ... left ends of equal number of exons */
    for (k = 0; k < p; ++k) {
      cds[k][0] = pgl->gcds[a][k][0];
      cds[k][1] = pgl->gcds[a][k][1];
      cdsl[k] = abs(cds[k][1]-cds[k][0]);
      ++n;
    }
    if (pgl->rflag == 0) {			/* ... adjust first exon left ends */
      if (cds[0][0] > gca->gcds[0][0])
	cds[0][0] = gca->gcds[0][0];
    }
    else {
      if (cds[0][0] < gca->gcds[0][0])
	cds[0][0] = gca->gcds[0][0];
    }
  }
  if (LTFLAG && n > 0) {
    fprintf(stdout,"\n  NOW1 (");
    for (i = 0; i < n - 1; ++i)
      fprintf(stdout,"%d  %d,", cds[i][0], cds[i][1]);
    fprintf(stdout,"%d  %d)", cds[i][0], cds[i][1]);
  }

  if (p == 0) {
    if (g == 0) {
      if (p == pgl->exn[a] - 1) {
	if (g == gca->exn -1) {
	  cds[n][0] = pgl->gcds[a][p][0];
	  cds[n][1] = pgl->gcds[a][p][1];
	  cdsl[n] = abs(cds[n][1]-cds[n][0]);
	  if (pgl->rflag == 0) {
	    if (cds[n][0] > gca->gcds[g][0])
	      cds[n][0] = gca->gcds[g][0];
	  }
	  else {
	    if (cds[n][0] < gca->gcds[g][0])
	      cds[n][0] = gca->gcds[g][0];
	  }
	  if (pgl->rflag == 0) {
	    if (cds[n][1] < gca->gcds[g][1])
	      cds[n][1] = gca->gcds[g][1];
	  }
	  else {
	    if (cds[n][1] > gca->gcds[g][1])
	      cds[n][1] = gca->gcds[g][1];
	  }
        }
        else {
	  cds[n][0] = gca->gcds[g][0];
	  cds[n][1] = gca->gcds[g][1];
	  cdsl[n] = abs(cds[n][1]-cds[n][0]);
	  if (pgl->rflag == 0) {
	    if (cds[n][0] > pgl->gcds[a][p][0])
	      cds[n][0] = pgl->gcds[a][p][0];
	  }
	  else {
	    if (cds[n][0] < pgl->gcds[a][p][0])
	      cds[n][0] = pgl->gcds[a][p][0];
	  }
	}
      }
      else {
	if (g == gca->exn -1) {
	  cds[n][0] = pgl->gcds[a][p][0];
	  cds[n][1] = pgl->gcds[a][p][1];
	  cdsl[n] = abs(cds[n][1]-cds[n][0]);
	  if (pgl->rflag == 0) {
	    if (cds[n][0] > gca->gcds[g][0])
	      cds[n][0] = gca->gcds[g][0];
	  }
	  else {
	    if (cds[n][0] < gca->gcds[g][0])
	      cds[n][0] = gca->gcds[g][0];
	  }
	}
	else {
	  cds[n][0] = pgl->gcds[a][p][0];
	  cds[n][1] = pgl->gcds[a][p][1];
	  cdsl[n] = abs(cds[n][1]-cds[n][0]);
	  if (pgl->rflag == 0) {
	    if (cds[n][0] > gca->gcds[g][0])
	      cds[n][0] = gca->gcds[g][0];
	  }
	  else {
	    if (cds[n][0] < gca->gcds[g][0])
	      cds[n][0] = gca->gcds[g][0];
	  }
	}
      }
    }
    else {
      if (p == pgl->exn[a] - 1) {
	if (g == gca->exn -1) {
	  cds[n][0] = gca->gcds[g][0];
	  cds[n][1] = gca->gcds[g][1];
	  cdsl[n] = abs(cds[n][1]-cds[n][0]);
	  if (pgl->rflag == 0) {
	    if (cds[n][1] < pgl->gcds[a][p][1])
	      cds[n][1] = pgl->gcds[a][p][1];
	  }
	  else {
	    if (cds[n][1] > pgl->gcds[a][p][1])
	      cds[n][1] = pgl->gcds[a][p][1];
	  }
	}
	else {
	  cds[n][0] = gca->gcds[g][0];
	  cds[n][1] = gca->gcds[g][1];
	  cdsl[n] = abs(cds[n][1]-cds[n][0]);
	}
      }
      else {
	if (g == gca->exn -1) {
	  cds[n][0] = gca->gcds[g][0];
	  cds[n][1] = pgl->gcds[a][p][1];
	  cdsl[n] = abs(cds[n][1]-pgl->gcds[a][p][0]);
	}
	else {
	  cds[n][0] = gca->gcds[g][0];
	  cds[n][1] = gca->gcds[g][1];
	  cdsl[n] = abs(cds[n][1]-cds[n][0]);
	}
      }
    }
  }
  else {
    if (g == 0) {
      if (p == pgl->exn[a] - 1) {
	if (g == gca->exn -1) {
	  cds[n][0] = pgl->gcds[a][p][0];
	  cds[n][1] = pgl->gcds[a][p][1];
	  cdsl[n] = abs(cds[n][1]-cds[n][0]);
	  if (pgl->rflag == 0) {
	    if (cds[n][1] < gca->gcds[g][1])
	      cds[n][1] = gca->gcds[g][1];
	  }
	  else {
	    if (cds[n][1] > gca->gcds[g][1])
	      cds[n][1] = gca->gcds[g][1];
	  }
	}
	else {
	  cds[n][0] = pgl->gcds[a][p][0];
	  cds[n][1] = gca->gcds[g][1];
	  cdsl[n] = abs(pgl->gcds[a][p][1]-cds[n][0]);
	}
      }
      else {
	if (g == gca->exn -1) {
	  cds[n][0] = pgl->gcds[a][p][0];
	  cds[n][1] = pgl->gcds[a][p][1];
	  cdsl[n] = abs(cds[n][1]-cds[n][0]);
	}
	else {
	  cds[n][0] = pgl->gcds[a][p][0];
	  cds[n][1] = pgl->gcds[a][p][1];
	  cdsl[n] = abs(cds[n][1]-cds[n][0]);
	}
      }
    }
    else {
      if (LTFLAG) fprintf(stdout,"\nadj_ags: THIS SHOULD NOT HAPPEN (p %d g %d\n",p,g);
    }
  }
  ++n;
  if (LTFLAG) {
    fprintf(stdout,"\n  NOW2 (");
    for (i = 0; i < n - 1; ++i)
      fprintf(stdout,"%d  %d,", cds[i][0], cds[i][1]);
    fprintf(stdout,"%d  %d)", cds[i][0], cds[i][1]);
  }

  for (i = 1, j = 1; i < pgl->exn[a] - p && j < gca->exn - g; ++i, ++j);
  if (LTFLAG)
    fprintf(stdout,"\nIPJGN i=%d (%d) p=%d j=%d (%d) g=%d n=%d", i, pgl->exn[a] - p, p, j, gca->exn - g, g, n);
  if (i < pgl->exn[a] - p) {			/* ... gca right end contained in pgl right end */
    for (k = 1; k < pgl->exn[a] - p; ++k) {
      cds[n][0] = pgl->gcds[a][p + k][0];
      cds[n][1] = pgl->gcds[a][p + k][1];
      cdsl[n] = abs(cds[n][1]-cds[n][0]);
      ++n;
    }
  }
  else if (j < gca->exn - g) {			/* ... pgl right end contained in gca right end */
    for (k = 1; k < gca->exn - g; ++k) {
      cds[n][0] = gca->gcds[g + k][0];
      cds[n][1] = gca->gcds[g + k][1];
      cdsl[n] = abs(cds[n][1]-cds[n][0]);
      ++n;
    }
  }
  else {					/* ... right ends of equal number of exons */
    for (k = 1; k < pgl->exn[a] - p; ++k) {
      cds[n][0] = pgl->gcds[a][p + k][0];
      cds[n][1] = pgl->gcds[a][p + k][1];
      cdsl[n] = abs(cds[n][1]-cds[n][0]);
      ++n;
    }
    if (pgl->rflag == 0) {			/* ... adjust right end of last exon */
      if (cds[n - 1][1] < gca->gcds[g + k - 1][1])
	cds[n - 1][1] = gca->gcds[g + k - 1][1];
    }
    else {
      if (cds[n - 1][1] > gca->gcds[g + k - 1][1])
	cds[n - 1][1] = gca->gcds[g + k - 1][1];
    }
  }

  if (LTFLAG) {
    fprintf(stdout,"\n  NEW (");
    for (i = 0; i < n - 1; ++i)
      fprintf(stdout,"%d  %d,", cds[i][0], cds[i][1]);
    fprintf(stdout,"%d  %d)", cds[i][0], cds[i][1]);
  }

  ja = 0;
  for (i = 0; i < n; ++i) {
    exnscr[i] = 0.0;
    for (j = ja; j < pgl->exn[a]; ++j) {
      if (pgl->gcds[a][j][0] == cds[i][0] || pgl->gcds[a][j][1] == cds[i][1]) {
	exnscr[i] = pgl->exnscr[a][j];
	if (i < n - 1 && j < pgl->exn[a] - 1) {
	  itrscr[i][0] = pgl->itrscr[a][j][0];
	  itrscr[i][1] = pgl->itrscr[a][j][1];
	}
	ja = j + 1;
	break;
      }
    }
  }
  ja = 0;
  for (i = 0; i < n; ++i) {
    for (j = ja; j < gca->exn; ++j) {
      if (gca->gcds[j][0] == cds[i][0] || gca->gcds[j][1] == cds[i][1]) {
	if ( abs(gca->gcds[j][1] - gca->gcds[j][0]) >  cdsl[i]+30     ||
	    (abs(gca->gcds[j][1] - gca->gcds[j][0]) >= cdsl[i]-30  &&
	     gca->exnscr[j] > exnscr[i]                              )  )
	  exnscr[i] = gca->exnscr[j];
					/* ... assign to a merged or duplicated exon the score of the
					   longest contributing exon, or the highest score if the
					   contributing exons are of similar lengths (equal +- 30 bases) */
	if (i < n - 1 && j < gca->exn - 1) {
	  itrscr[i][0] = gca->itrscr[j][0];
	  itrscr[i][1] = gca->itrscr[j][1];
	}
	ja = j + 1;
	break;
      }
    }
  }

  if (LTFLAG) {
    fprintf(stdout,"\n  SCR (");
    for (i = 0; i < n - 1; ++i)
      fprintf(stdout,"e %5.3f  d %5.3f a %5.3f,", exnscr[i], itrscr[i][0], itrscr[i][1]);
    fprintf(stdout,"e %4.2f)", exnscr[i]);
    fprintf(stdout,"\n\n");
  }

  pgl->exn[a] = n;
  for (i = 0; i < n; ++i) {
    pgl->gcds[a][i][0] = cds[i][0];
    pgl->gcds[a][i][1] = cds[i][1];
    pgl->exnscr[a][i] = exnscr[i];
    if (i < n - 1) {
      pgl->itrscr[a][i][0] = itrscr[i][0];
      pgl->itrscr[a][i][1] = itrscr[i][1];
    }
  }
  gca->ags[a] = 1;

  return (1);

}							/* end adj_ags() */



int new_ags(pgl, a, gca)
  struct pgl *pgl;
  int a;
  struct gcalgnmnt *gca;
{
  int i;

  pgl->nags = a + 1;
  pgl->exn[a] = gca->exn;
  for (i = 0; i < gca->exn; ++i) {
    pgl->gcds[a][i][0] = gca->gcds[i][0];
    pgl->gcds[a][i][1] = gca->gcds[i][1];
    pgl->exnscr[a][i] = gca->exnscr[i];
    if (i < gca->exn - 1) {
      pgl->itrscr[a][i][0] = gca->itrscr[i][0];
      pgl->itrscr[a][i][1] = gca->itrscr[i][1];
    }
  }
  gca->ags[a] = 1;
  if (LTFLAG) fprintf(stdout,"new AGS-%d created\n\n", pgl->nags);
  return 0;
}							/* end new_ags() */



int consolidate_pgl(struct pgl **pglhpp, int *npgl)
{
  struct pgl *ip = *pglhpp, *tmp, *newheadp = NULL;
  int rval = 0;

  *npgl = 0;
  while (ip != NULL) {
    tmp = ip->next;
    ip->next = NULL;
    if (insert_pgl(&newheadp, ip) == 0) {
      if (LTFLAG) fprintf(stdout,"\nDEALLOC_PGL %ld # %d NO-INSERT IN consolidate_pgl\n",(long)ip,++NDEALLOC_PGL);
      free(ip);
      rval = 1;
    }
    else
      ++(*npgl);
    ip = tmp;
  }
  *pglhpp = newheadp;
  return (rval);

}							/* end consolidate_pgl() */



int prt_consensus_pgl(FILE *fp, struct pgl *pgl, char *gname, char *gsgmntn, char *gdna, char *gdnaR, int numbp, int npgl, int bflag)
{
  int i, j, pgln = 0;
  struct pgl *headpgl;
  struct gcalgnmnt *pgca;
  char anchorName[257];
  extern int htmlop;
  extern int frompA, top;

  if (htmlop) {
    strcpy(anchorName,gsgmntn);
#ifdef HTMLWS
    fprintf(imageDataFh,"%d\n",npgl);
    fprintf(imageDataFh,"%d %d\n",frompA,top);
#endif
  }

  if (npgl == 0)   return 0;

  if (htmlop) {
    sprintf(anchorName,"PGL%d@%s",1,gsgmntn);
    ADD_NAME(fp,anchorName);
#ifdef HTMLWS
    fprintf(imageDataFh,"PGL %*d %*d %s\n",ifwdth, pgl->gcrd[0], ifwdth, pgl->gcrd[1],anchorName);
#endif
  }
  fprintf(fp, "\n\nPredicted gene locations (%d) in segment %d to %d:\n\n", npgl, frompA+1, top+1);
  if (htmlop) {
    fprintf(fp, "\nScroll up   to     <A HREF=\"#TOP\">top</A>");
    fprintf(fp, "\nScroll down to PGL");
    headpgl = pgl;
    while (pgl != NULL) {
      sprintf(anchorName,"PGL%d@%s",pgln+1,gsgmntn);
      fprintf(fp, " <A HREF=\"#%s\">%d</A>",anchorName,pgln+1);
      ++pgln;
      pgl = pgl->next;
    }
    fprintf(fp, "\n");
    pgl = headpgl;
    pgln = 0;
  }

  while (pgl != NULL) {
    if (htmlop) {
      if (pgln > 0) {
        sprintf(anchorName,"PGL%d@%s",pgln+1,gsgmntn);
        ADD_NAME(fp,anchorName);
#ifdef HTMLWS
        fprintf(imageDataFh,"PGL %*d %*d %s\n",ifwdth, pgl->gcrd[0], ifwdth, pgl->gcrd[1],anchorName);
#endif
      }
      fprintf(fp, "\n<A HREF=\"#HEAD-PGL-%s\">Scroll back to \"Predicted gene locations\"</A>",gsgmntn);
      fprintf(fp, "\n<A HREF=\"#BOTTOM-PGL-%s\">Scroll down to next segment</A>\n",gsgmntn);
    }

    fprintf(fp, "\nPGL %3d (", ++pgln);
    if (pgl->rflag == 0)
      fprintf(fp, "+ strand");
    else
      fprintf(fp, "- strand");
    fprintf(fp, "):\t%*d  %*d", ifwdth, pgl->gcrd[0], ifwdth, pgl->gcrd[1]);


    for (i = 0; i < pgl->nags; ++i) {

      if (htmlop) {
        sprintf(anchorName,"PGL%d_AGS%d@%s",pgln,i,gsgmntn);
        ADD_NAME(fp,anchorName);
      }

      fprintf(fp, "\nAGS-%d (", i + 1);
      for (j = 0; j < pgl->exn[i] - 1; ++j)
	fprintf(fp, "%d  %d,", pgl->gcds[i][j][0], pgl->gcds[i][j][1]);
      fprintf(fp, "%d  %d)", pgl->gcds[i][j][0], pgl->gcds[i][j][1]);

#ifdef HTMLWS
      if (htmlop) {
        fprintf(imageDataFh,"AGS %*d %*d %s\n",ifwdth, pgl->gcds[i][0][0],
		ifwdth, pgl->gcds[i][j][1],anchorName);
        for (j=0;j<pgl->exn[i]-1;j++)
	  fprintf(imageDataFh, "%d %d ", pgl->gcds[i][j][0], pgl->gcds[i][j][1]);
	fprintf(imageDataFh, "%d %d\n", pgl->gcds[i][j][0], pgl->gcds[i][j][1]);
      }
#endif

      fprintf(fp, "\nSCR   (");
      for (j = 0; j < pgl->exn[i] - 1; ++j)
	fprintf(fp, "e %5.3f  d %5.3f a %5.3f,", pgl->exnscr[i][j], pgl->itrscr[i][j][0], pgl->itrscr[i][j][1]);
      fprintf(fp, "e %5.3f)\n", pgl->exnscr[i][j]);

      fprintf(fp, "\n");
      for (j = 0; j < pgl->exn[i] - 1; ++j) {
	fprintf(fp, "  Exon %2d %*d %*d (%4d n); score: %5.3f\n", j+1, ifwdth, pgl->gcds[i][j][0], ifwdth, pgl->gcds[i][j][1], abs(pgl->gcds[i][j][1]- pgl->gcds[i][j][0])+1, pgl->exnscr[i][j]);
        if (pgl->gcds[i][0][0] < pgl->gcds[i][pgl->exn[i]-1][1])
	  fprintf(fp, "   Intron %2d %*d %*d (%4d n);           Pd: %5.3f  Pa: %5.3f\n", j+1, ifwdth, pgl->gcds[i][j][1]+1, ifwdth, pgl->gcds[i][j+1][0]-1, pgl->gcds[i][j+1][0]-pgl->gcds[i][j][1]-1, pgl->itrscr[i][j][0], pgl->itrscr[i][j][1]);
        else
	  fprintf(fp, "   Intron %2d %*d %*d (%4d n);           Pd: %5.3f  Pa: %5.3f\n", j+1, ifwdth, pgl->gcds[i][j][1]-1, ifwdth, pgl->gcds[i][j+1][0]+1, pgl->gcds[i][j][1]-pgl->gcds[i][j+1][0]-1, pgl->itrscr[i][j][0], pgl->itrscr[i][j][1]);
      }
      fprintf(fp, "  Exon %2d %*d %*d (%4d n); score: %5.3f\n", j+1, ifwdth, pgl->gcds[i][j][0], ifwdth, pgl->gcds[i][j][1], abs(pgl->gcds[i][j][1]- pgl->gcds[i][j][0])+1, pgl->exnscr[i][j]);

      pgca = pgl->gca;
      while (pgca != NULL) {
	if (pgca->ags[i] == 1) {
	  if (htmlop) {
            sprintf(anchorName,"PGS%d@%s",pgca->calln,gsgmntn);
            ADD_PGL_LINK(fp,anchorName);
	  }
	  else
	    fprintf(fp, "\n  PGS (");
	  for (j = 0; j < pgca->exn - 1; ++j)
	    fprintf(fp, "%d  %d,", pgca->gcds[j][0], pgca->gcds[j][1]);
	  if (htmlop) {
#ifdef PLANTGDB
	    fprintf(fp, "%d  %d)\t%s", pgca->gcds[j][0], pgca->gcds[j][1], getDnaLink(2,pgca->cname));
#else
	    fprintf(fp, "%d  %d)\t%s", pgca->gcds[j][0], pgca->gcds[j][1], getDnaLink(1,pgca->cname));
#endif
#ifdef HTMLWS
            if (htmlop) {
              fprintf(imageDataFh,"PGS %*d %*d %s\n",ifwdth, pgca->gcds[0][0],
			ifwdth, pgca->gcds[j][1],anchorName);
              for (j = 0; j < pgca->exn - 1; ++j)
	        fprintf(imageDataFh, "%d %d ", pgca->gcds[j][0], pgca->gcds[j][1]);
	      fprintf(imageDataFh, "%d %d\n", pgca->gcds[j][0],pgca->gcds[j][1]);
            }
#endif
	  }
	  else
	    fprintf(fp, "%d  %d)\t%s", pgca->gcds[j][0], pgca->gcds[j][1], pgca->cname);
	}
	pgca = pgca->link;
      }
      fprintf(fp,"\n\n");
      OutputORF(fp, gdna, gdnaR, gname, gsgmntn, numbp, pgl->rflag, pgl, pgln, i);
      if (bflag == 1  &&  pgl->exn[i] == 1) {
	reverse_pgl(pgl,i,numbp);
	OutputORF(fp, gdna, gdnaR, gname, gsgmntn, numbp, pgl->rflag, pgl, pgln, i);
	reverse_pgl(pgl,i,numbp);
      }
    }
    fprintf(fp, "\n");
    if (htmlop) {
      fprintf(fp, "\n<A HREF=\"#HEAD-PGL-%s\">Scroll back to \"Predicted gene locations\"</A>",gsgmntn);
      fprintf(fp, "\n<A HREF=\"#TOP\">Scroll up to top</A>\n");
    }
    pgl = pgl->next;
  }
  return 0;
}							/* end prt_consensus_pgl() */



void free_pgl(struct pgl *pgl)
{
  struct pgl *tpgl = pgl;

  while (tpgl != NULL) {
    tpgl = pgl->next;
    if (LTFLAG) fprintf(stdout,"\nDEALLOC_PGL %ld # %d",(long)pgl,++NDEALLOC_PGL);
    free((struct pgl *) pgl);
    pgl = tpgl;
  }
 if (LTFLAG) fprintf(stdout,"\n");

}							/* end free_pgl() */



void OutputORF(FILE *fp, char *gdna, char *gdnaR, char *gname, char *gsgmntn, int numbp, int rflag, struct pgl *pgl, int pgln, int a)
{
  int n, ni, nic, i, j, l = 0, lj[MAXNEXNS];
  char *seq;
  int pa;

  if (rflag == 1) {
    for (n = 0; n < pgl->exn[a]; ++n)   l+= (pgl->gcds[a][n][0]-pgl->gcds[a][n][1]+1);
    if ((seq = (char *) calloc((l + 1) , sizeof(char))) == NULL)
      fatal_error("Error: memory allocation failed. Exit.\n");
    l = 0;
    for (n = 0; n < pgl->exn[a]; ++n) {
      for (i = pgl->gcds[a][n][0]; i >= pgl->gcds[a][n][1]; --i)
	seq[l++] = gdnaR[numbp - i];
      lj[n] = l - 1;
    }
    lj[n - 1] = -999;
  }
  else {
    for (n = 0; n < pgl->exn[a]; ++n)   l+= (pgl->gcds[a][n][1]-pgl->gcds[a][n][0]+1);
     if ((seq = (char *) calloc((l + 1) , sizeof(char))) == NULL)
       fatal_error("Error: memory allocation failed. Exit.\n");
    l = 0;
    for (n = 0; n < pgl->exn[a]; ++n) {
      for (i = pgl->gcds[a][n][0] - 1; i < pgl->gcds[a][n][1]; ++i)
	seq[l++] = gdna[i];
      lj[n] = l - 1;
    }
    lj[n - 1] = -999;
  }

  pa = 0;
  n = 0;
  ni = 0;
  fprintf(fp, "3-phase translation of AGS-%d", a + 1);
  if (rflag == 1) fprintf(fp, " (-strand):\n");
  else            fprintf(fp, " (+strand):\n");
  for (i = 0; i <= l; ++i) {
    if (i % 60 == 0 || i == l) {
      if (i > 0) {
	fprintf(fp, "\n  ");
        for (j=0;j<ifwdth;++j) fprintf(fp," ");
	nic = ni;
	for (j = pa + 1; j <= pa + 58 && j + 1 <= l - 1; j += 3) {
	  if (j - 1 == lj[n - nic]) {
	    fprintf(fp, "  : ");
	    --nic;
	  }
	  else if (j - 2 == lj[n - nic]) {
	    fprintf(fp, " :  ");
	    --nic;
	  }
	  else if (j - 3 == lj[n - nic]) {
	    fprintf(fp, ":   ");
	    --nic;
	  }
	  else
	    fprintf(fp, " ");
	  if (seq[j - 1] > 3 || seq[j] > 3 || seq[j + 1] > 3)
	    fprintf(fp, "%c", AAUC[22]);
	  else
	    fprintf(fp, "%c", CODON2CHAR(j));
	  if (j == lj[n - nic]) {
	    fprintf(fp, " :  ");
	    --nic;
	  }
	  else if (j + 1 == lj[n - nic]) {
	    fprintf(fp, "  : ");
	    --nic;
	  }
	  else
	    fprintf(fp, " ");
	}
	fprintf(fp, "\n   ");
        for (j=0;j<ifwdth;++j) fprintf(fp," ");
	nic = ni;
	if (lj[n - nic] == pa) {
	  fprintf(fp, " : ");
	  --nic;
	}
	for (j = pa + 2; j <= pa + 59 && j + 1 <= l - 1; j += 3) {
	  if (j - 1 == lj[n - nic]) {
	    fprintf(fp, "  : ");
	    --nic;
	  }
	  else
	    fprintf(fp, " ");
	  if (seq[j - 1] > 3 || seq[j] > 3 || seq[j + 1] > 3)
	    fprintf(fp, "%c", AAUC[22]);
	  else
	    fprintf(fp, "%c", CODON2CHAR(j));
	  if (j == lj[n - nic]) {
	    fprintf(fp, " :  ");
	    --nic;
	  }
	  else if (j + 1 == lj[n - nic] && j + 1 < pa + 60) {
	    fprintf(fp, "  : ");
	    --nic;
	  }
	  else if (j + 2 == lj[n - nic] && j + 2 < pa + 60) {
	    fprintf(fp,  "   :");
	    --nic;
	  }
	  else
	    fprintf(fp, " ");
	}
	fprintf(fp, "\n ");
        for (j=0;j<ifwdth;++j) fprintf(fp," ");
	nic = ni;
	for (j = pa; j <= pa + 57 && j + 1 <= l - 1; j += 3) {
	  if (j == 0) {
	    if (j == lj[n - nic])
	      fprintf(fp, "   :  ");
	    else
	      fprintf(fp, "   ");
	  }
	  else {
	    if (j - 1 == lj[n - nic]) {
	      fprintf(fp, "  : ");
	      --nic;
	    }
	    else if (j - 2 == lj[n - nic]) {
	      fprintf(fp,  " :  ");
	      --nic;
	    }
	    else
	      fprintf(fp, " ");
	    if (seq[j - 1] > 3 || seq[j] > 3 || seq[j + 1] > 3)
	      fprintf(fp, "%c", AAUC[22]);
	    else
	      fprintf(fp, "%c", CODON2CHAR(j));
	    if (j == lj[n - nic]) {
	      fprintf(fp, " :  ");
	      --nic;
	    }
	    else if (j + 1 == lj[n - nic]) {
	      fprintf(fp, "  : ");
	      --nic;
	    }
	    else if (j + 2 == lj[n - nic]) {
	      fprintf(fp,  "   :");
	      --nic;
	    }
	    else
	      fprintf(fp, " ");
	  }
	}
	fprintf(fp, "\n\n");
	ni = 0;
      }
      if (i == l)
	break;
      pa = i;
      fprintf(fp, "\n  ");
      for (j=0;j<ifwdth;++j) fprintf(fp," ");
      nic = n;
      for (j = i; j < i + 60 && j < l; ++j) {
	if ((j - i) % 10 == 0)
	  fprintf(fp, ".");
	else
	  fprintf(fp, " ");
	if (j == lj[nic] && j < l - 2) {
	  fprintf(fp, " : ");
	  ++nic;
	}
      }
      fprintf(fp, "\n%*d  ", ifwdth, crdfct(i, pgl->exn[a], pgl->gcds[a], numbp, rflag));
    }
    fprintf(fp, "%c", NAUC[(int)seq[i]]);
    if (i == lj[n] && i < l - 1) {
      fprintf(fp, " : ");
      ++n;
      ++ni;
    }
  }
  fprintf(fp, "\n");
  find_orfs(fp, seq, gname, gsgmntn, l, 64, numbp, rflag, pgl, pgln, a);
  fprintf(fp, "\n\n");
  free(seq);

}							/* end OutputORF() */



void reverse_pgl(struct pgl *pgl, int a, int numbp)
{
  int i, tmpa, tmpb;

  pgl->rflag = MOD2(pgl->rflag + 1);
  for (i = pgl->exn[a] - 1; i >= DIV2(pgl->exn[a] + 1); --i) {
    tmpa= pgl->gcds[a][i][1];   tmpb= pgl->gcds[a][i][0];
    pgl->gcds[a][i][0]= pgl->gcds[a][pgl->exn[a] - 1 - i][1];
    pgl->gcds[a][i][1]= pgl->gcds[a][pgl->exn[a] - 1 - i][0];
    pgl->gcds[a][pgl->exn[a] - 1 - i][0]= tmpa;
    pgl->gcds[a][pgl->exn[a] - 1 - i][1]= tmpb;

  }
  if (MOD2(pgl->exn[a] ) == 1) {
    tmpa= pgl->gcds[a][DIV2(pgl->exn[a])][0];
    pgl->gcds[a][DIV2(pgl->exn[a])][0]= pgl->gcds[a][DIV2(pgl->exn[a])][1];
    pgl->gcds[a][DIV2(pgl->exn[a])][1]= tmpa;
  }

}				/* end reverse_pgl() */



void find_orfs(FILE *fp, char *seq, char *sname, char *sgmntn, int slgth, int minorfl, int numbp, int rflag, struct pgl *pgl, int pgln, int a)
 /* finds the maximal open reading frames >= minorfl in seq */
{
  int i, j, aa, n = 0;
  int last_stp[3];
  int frame, l_orf;
  struct orf *orf, *orfheadp = NULL;

  fprintf(fp, "Maximal non-overlapping open reading frames (>= %d codons):\n", minorfl);

  for (i = 0; i <= 2; ++i)
    last_stp[i] = -1;
  for (i = 0; i <= slgth - minorfl * 3; ++i) {
    if (seq[i] > 3 || seq[i + 1] > 3 || seq[i + 2] > 3)
      aa = 22;
    else
      aa = CODON2NUMBER(i+1);
    frame = (i + 1) % 3;
    if (i <= last_stp[frame])
      continue;
    j = i;
    while (aa != 23) {
      j += 3;
      if (j > slgth - 3)
	break;
      if (seq[j] > 3 || seq[j + 1] > 3 || seq[j + 2] > 3)
	aa = 22;
      else
	aa = CODON2NUMBER(j+1);
    }
    last_stp[frame] = j;
    l_orf = (j - i) / 3;
    if (l_orf >= minorfl) {
      if ((orf = (struct orf *) calloc(1 , sizeof(struct orf))) == NULL)
        fatal_error("Error: memory allocation failed. Exit.\n");
      orf->i = i;
      orf->j= j;
      orf->l = l_orf;
      orf->f = frame;
      orf->n = ++n;
      orf->next = NULL;
      insert_orf(&orfheadp,orf);
    }

  }

  if (n == 0) {
    fprintf(fp, "none\n");
    pgl->npps[a] = 0;
  }
  else {
    consolidate_orfs(&orfheadp);
    orf = orfheadp;
    while (orf != NULL) {
      trl_orf(fp, seq, sname, sgmntn, slgth, orf->i, orf->j, orf->f, numbp, rflag, orf->n, pgl, pgln, a);
      orf = orf->next;
    }
  }
  free_orf(orfheadp);

}							/* end find_orfs() */



int insert_orf(orfhpp, nodep)
  struct orf **orfhpp, *nodep;
{
  if (*orfhpp == NULL) {
	*orfhpp = nodep;
	nodep->next = NULL;
  }
  else {
	while ((*orfhpp)->next != NULL)
	 {
	  if ((*orfhpp)->l > nodep->l)
            orfhpp = &((*orfhpp)->next);
	  else
	    break;
         }
        if ((*orfhpp)->next == NULL) {
		if ((*orfhpp)->l > nodep->l) {(*orfhpp)->next = nodep; nodep->next = NULL;}
		else
		  {nodep->next = *orfhpp; *orfhpp = nodep;}
          }
        else
          { nodep->next = *orfhpp; *orfhpp = nodep; }
  }
  return (1);

}							/* end insert_orf() */



int consolidate_orfs(struct orf **orfhpp)
{
  struct orf *orf = *orfhpp, *tmp, *newheadp = NULL;
  int i, n = 0;
  int crd[100][2];
  int minoverlap = 30, deltares = 30;

  crd[n][0]= orf->i; crd[n][1]= orf->j; ++n;
  tmp = orf->next;
  orf->n = 1;
  orf->next = NULL;
  insert_orf(&newheadp, orf);
  orf = tmp;
  while (orf != NULL) {
    for (i=0;i<n;++i) {
      if (is_xwab(orf->i,crd[i][0],crd[i][1],-minoverlap)              ||
	  is_xwab(orf->j,crd[i][0],crd[i][1],-minoverlap)              ||
	  xy_contains_ab(orf->i,orf->j,crd[i][0],crd[i][1],-minoverlap)  ) 
 	break;
    }
    if (i == n  ||  crd[0][1]-crd[0][0] - 3*orf->l < deltares) {
		/* ... a sub-maximal length ORF is listed if its length
		is within 'deltares' codons of the maximal ORF length or if it
		does does not share at least 'minoverlap' nucleotides with an already
		accepted ORF */
      crd[n][0]= orf->i;
      crd[n][1]= orf->j;
      ++n;
      orf->n = n;
      tmp = orf->next;
      orf->next = NULL;
      insert_orf(&newheadp, orf);
      orf = tmp;
    }
    else{
      tmp = orf;
      orf = orf->next;
      free((struct orf *)tmp);
    }
  }
  *orfhpp = newheadp;
  return 0;
}							/* end consolidate_orfs() */



void trl_orf(FILE *fp, char *seq, char *sname, char *sgmntn, int slgth, int firstbp, int lastbp, int frame, int numbp, int rflag, int orfn, struct pgl *pgl, int pgln, int a)
{
  int i, l, ip, n, psg;
  int olgth, numaa;
  char cod[MAXGLGTH], protein[MAXGLGTH];
  char anchorName[257];
  static int orfCount=0;
  int pos[200], ind=0, addstop=0;

  orfCount++;

  olgth = lastbp - firstbp;
  numaa = olgth / 3;

  psg = 0;
  i = crdfct(firstbp, pgl->exn[a], pgl->gcds[a], numbp, rflag);
  if (orfn <= MAXNPPS)
    pgl->pcds[a][orfn - 1][psg][0] = i;
  if (rflag) sname[strlen(sname) - 1] = '-';
  else       sname[strlen(sname) - 1] = '+';

  if (htmlop) {
    sprintf(anchorName, "PGL%d_ORF%d@%s",pgln,orfCount,sgmntn);
    ADD_NAME(fp,anchorName);
  }
  fprintf(fp, "\n>%s_PGL-%d_AGS-%d_PPS_%d (%d", sname, pgln, a + 1, orfn, i);
  pos[ind++]=i;

  if (lastbp < slgth - 2) {
    lastbp += 3; /*INCLUDE THE STOP CODON*/
    addstop = 1;
  }

  l = crdfct(lastbp - 1, pgl->exn[a], pgl->gcds[a], numbp, rflag);
  if (rflag) {
    for (n = 0; n < pgl->exn[a]; ++n) {
      if (i >= pgl->gcds[a][n][1] && l < pgl->gcds[a][n][1]) {
	fprintf(fp, "  %d,", pgl->gcds[a][n][1]);
        pos[ind++]=pgl->gcds[a][n][1];
	if (orfn <= MAXNPPS)
	  pgl->pcds[a][orfn - 1][psg++][1] = pgl->gcds[a][n][1];
	if (l <= pgl->gcds[a][n + 1][0]) {
	  fprintf(fp, "%d", pgl->gcds[a][n + 1][0]);
          pos[ind++]=pgl->gcds[a][n + 1][0];
	  if (orfn <= MAXNPPS)
	    pgl->pcds[a][orfn - 1][psg][0] = pgl->gcds[a][n + 1][0];
	}
      }
    }
  }
  else {
    for (n = 0; n < pgl->exn[a]; ++n) {
      if (i <= pgl->gcds[a][n][1] && l > pgl->gcds[a][n][1]) {
	fprintf(fp, "  %d,", pgl->gcds[a][n][1]);
        pos[ind++]=pgl->gcds[a][n][1];
	if (orfn <= MAXNPPS)
	  pgl->pcds[a][orfn - 1][psg++][1] = pgl->gcds[a][n][1];
	if (l >= pgl->gcds[a][n + 1][0]) {
	  fprintf(fp, "%d", pgl->gcds[a][n + 1][0]);
          pos[ind++]=pgl->gcds[a][n + 1][0];
	  if (orfn <= MAXNPPS)
	    pgl->pcds[a][orfn - 1][psg][0] = pgl->gcds[a][n + 1][0];
	}
      }
    }
  }
  fprintf(fp, "  %d)", l);

  if (htmlop) {
    pos[ind++]=l;
#ifdef HTMLWS
    fprintf(imageDataFh,"ORF %10d %10d %s\n",i,l,anchorName);
    for (i=0;i<ind-1;i++)
      fprintf(imageDataFh,"%d ",pos[i]);
    fprintf(imageDataFh,"%d\n",pos[i]);
#endif
  }

  if (orfn <= MAXNPPS) {
    pgl->pcds[a][orfn - 1][psg++][1] = l;
    pgl->npsg[a][orfn - 1] = psg;
    pgl->npps[a] = orfn;
  }
  fprintf(fp, "\t(frame '%1d'; %5d bp, %4d residues)", frame, olgth, numaa);

  for (i = firstbp; i < lastbp; ++i)
    cod[i - firstbp] = seq[i];
  if (addstop) ++numaa;

  i = 0;
  for (ip = 0; ip < numaa; ++ip) {
    if (cod[i] > 3 || cod[i + 1] > 3 || cod[i + 2] > 3)
      protein[ip] = 22;
    else
      protein[ip] = codtoaa[(int)cod[i]][(int)cod[i + 1]][(int)cod[i + 2]];
    i += 3;
  }

  for (i = 0; i < numaa; ++i) {
    if (i % 10 == 0)
      fprintf(fp, " ");
    if (i % 60 == 0)
      fprintf(fp, "\n%*d  ", ifwdth, i + 1);
    fprintf(fp, "%c", AAUC[(int)protein[i]]);
  }
  fprintf(fp, "\n");
  if (htmlop) {
    fprintf(fp, "\n<A HREF=\"%s",BLASTP_HEAD);
    for(i=0;i<numaa;++i) if (protein[i]<24) fprintf(fp, "%c", AAUC[(int)protein[i]]);
    fprintf(fp, "%s\" TARGET=\"NCBI Blastp\">NCBI Blastp</A>",
        BLASTP_TAIL);
  }

}							/* end trl_orf() */



void free_orf(struct orf *orf)
{
  struct orf *torf = orf;

  while (torf != NULL) {
    torf = orf->next;
    free((struct orf *) orf);
    orf = torf;
  }

}							/* end free_orf() */



int crdfct(int n, int nsgmts, int iab[MAXNEXNS][2], int numbp, int rflag)
{
  int i, r = 0;

  if (rflag) {
    for (i = 0; i < nsgmts; ++i) {
      if (n < r + (iab[i][0] - iab[i][1] + 1))
	break;
      r += (iab[i][0] - iab[i][1] + 1);
      if (n == r)
	break;
    }
    if (n>0  &&  n == r  &&  i<nsgmts-1)
      return (iab[i+1][0]);
    else
      return (iab[i][0] - n + r);
  }
  else {
    for (i = 0; i < nsgmts; ++i) {
      if (n < r + (iab[i][1] - iab[i][0] + 1))
	break;
      r += (iab[i][1] - iab[i][0] + 1);
      if (n == r)
	break;
    }
    if (n>0  &&  n == r  &&  i<nsgmts-1)
      return (iab[i+1][0]);
    else
      return (iab[i][0] + n - r);
  }

}							/* end crdfct() */
