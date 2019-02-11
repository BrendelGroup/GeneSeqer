/* sahmtp.c
   Sequence Alignment Hidden Markov Tool - target Protein version

   Last update:  February 11, 2019. (VB)
*/


/* Corresponding author:                                                      */

/*   Volker Brendel, Department of Biology                                    */
/*   Indiana University, Bloomington, IN 47405                                */
/*   (812) 855-7074, vbrendel@indiana.edu                                     */

/* Past contributing authors:                                                 */
/*   Jonathan Usuka, Department of Chemistry, Stanford University             */
/*   Wei Zhu, Department of Zoology & Genetics, Iowa State University         */
/*   Fred Goodman, VisualMetrics Corporation                                  */
/*   George Juras, VisualMetrics Corporation                                  */
/*   Michael Sparks, Department of Genetics, Development and Cell Biology,    */
/*    Iowa State University                                                   */

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



/*      INCLUSIONS      */

#define MATHFUNC					/* Gives access to real math functions */
#define MEMOFUNC					/* Gives access to memory allocation */
#define STRGFUNC					/* Gives access to string manipulation */
#include "platform.h"

#include "sahmt.h"
#include "minmax.h"
#include "html.h"
#include "strapp.h"


void prt_seq_segment_to_str(char *str, char *seq, int numbp, int ia, int ib, int rflag, char ABC[], int nflag);

extern char AAUC[], NAUC[];
extern int htmlop, ifwdth;

int cod2aa[4][4][4] = {
 { {13, 13,  0,  0}, { 3,  3,  3,  3}, {15, 15, 23, 23}, {18, 18, 23, 19} },
 { { 0,  0,  0,  0}, {11, 11, 11, 11}, {16, 16, 14, 14}, {10, 10, 10, 10} },
 { { 9,  9,  9, 17}, { 7,  7,  7,  7}, {12, 12,  5,  5}, { 3,  3, 10, 10} },
 { { 4,  4,  4,  4}, { 1,  1,  1,  1}, { 8,  8,  6,  6}, { 2,  2,  2,  2} } };
float SCALEF = 0.4F;


/*      ^^^^^^^^^ FUNCTION DEFINITIONS ^^^^^^^^^ */


int sahmtP (FILE* fp, char *gdna, char *gname, int glgth, int ia,
   int numbp, int rflag, float *pd, float *pa, char *protein,
   char *pname, int plgth, int smat[23][23],
   struct gpalgnmnt *gpa, struct sprmtr sprm, int calln);

static void CompletePathMatrix(int glgth, int plgth, char *gdna, char *protein, float *pd, float *pa);
static void C_0m(int plgth);
static void I_0m(int plgth);
static void C_nm(int n, int m, char *gdna, char *protein, int glgth, int plgth, float *pd, float *pa);
static void I_nm(int n, int m, char *gdna, float *pd, float *pa);
static void IA_nm(int n, int m, float *pd, float *pa);
static void IB_nm(int n, int m, char *gdna, float *pd, float *pa);
static void IC_nm(int n, int m, char *gdna, float *pd, float *pa);


static void OptimalPath(int glgth, int plgth, char *gdna, char *protein);
static void TracePath(int state, char *gdna, char *protein, int plgth);
static void WritePath(int g1, int g2, int g3, int aa, int state);
static void WriteCodon(int g1, int g2, int g3, int aa);
static void WriteIntron(int g, int state);

static void OutputEnds(int ia);
static void OutputParse(float *pd, float *pa, int ia, int numbp,
  int rflag, char *gname, int glgth, char *pname, int plgth, struct gpalgnmnt *gpa);
static void OutputAlignment(int numbp, int rflag, struct gpalgnmnt *gpa);
static void OutputBase(int base);
static void OutputResidue(int aa);

static void allocate_Pspace(int glgth, int plgth);
static void free_Pspace(int glgth, int plgth);

int insert_gpa(struct gpalgnmnt **gpahpp,struct gpalgnmnt *nodep);
int prt_gpa_list(FILE *fp,struct gpalgnmnt *gpap, char *gsgmntn);


static int IndexMax(float *x, int n);
static float GetScore(int g1, int g2, int g3, int aa);
static float weight(int state, int path, int n, int glgth, int m,
  int plgth, int g1, int g2, int g3, int aa, float *pd, float *pa);
extern void fatal_error(char *buf);
int compare_gpas(struct gpalgnmnt *gpa1p, struct gpalgnmnt *gpa2p);




/*      ^^^^^^^^^ FUNCTION IMPLEMENTATIONS ^^^^^^^^^ */

static float *Log_Pd,*Log_1_Pd,*Log_Pa,*Log_1_Pa;

int sahmtP(FILE *outfp, char *gdna, char *gname, int glgth, int ia,
  int numbp, int rflag, float *pd, float *pa, char *protein,
  char *pname, int plgth, int smat[23][23],
  struct gpalgnmnt *gpa, struct sprmtr sprm, int calln)
{
  int i, j;
  extern int logPia;
  extern float *logV[8];
  int offset=ia-logPia;

  if (glgth > MAXGLGTH) {
    fprintf(outfp, "\nWARNING: %s gDNA length (%d) exceeds allowed maximum (%d).\n", gname, glgth, MAXGLGTH);
    return (0);
  }
  if (plgth > MAXPLGTH) {
    fprintf(outfp, "\nWARNING: %s protein length (%d) exceeds allowed maximum (%d).\n", pname, plgth, MAXPLGTH);
    return (0);
  }

  nextFree=0;
  AppendChar('\n');

  for (i = 0; i < 80; ++i)
    AppendChar('*');
  if (htmlop) {
    CHECKALIGNLEN; nextFree+=sprintf(global_algnmnt+nextFree, "\nQuery protein sequence %4d (File: %s)\n", calln, getProteinLink(pname));
  }
  else {
    CHECKALIGNLEN; nextFree+=sprintf(global_algnmnt+nextFree, "\nQuery protein sequence %4d (File: %s)\n", calln, pname);
  }
  prt_seq_segment_to_str(global_algnmnt, protein, plgth, 0, plgth - 1, 0, AAUC, 1);
  gpa->calln = calln;

  nextFree=strlen(global_algnmnt);

  allocate_Pspace(glgth, plgth);

  JNLGTH = sprm.join_length;
  for (i = 0; i < 23; ++i)
    for (j = 0; j < 23; ++j)
      ssmat[i][j] = smat[i][j];
  MIN_INTRONLENGTH = sprm.min_intron_length;
  SIWEIGHT = -100.0;

  if (rflag == 0) {
    Log_Pd=logV[0]+offset;
    Log_1_Pd=logV[1]+offset;
    Log_Pa=logV[2]+offset;
    Log_1_Pa=logV[3]+offset;
  }
  else {
    Log_Pd=logV[4]+offset;
    Log_1_Pd=logV[5]+offset;
    Log_Pa=logV[6]+offset;
    Log_1_Pa=logV[7]+offset;
  }


  /*      FILL OUT THE STATE AND PATH MATRICES: */

  CompletePathMatrix(glgth, plgth, gdna, protein, pd, pa);

  /*       BACK TRACE AN OPTIMAL PATH: */

  OptimalPath(glgth, plgth, gdna, protein);

  /*      PRINT OUT AN OPTIMAL ALIGNMENT: */

  OutputEnds(ia);

  OutputParse(pd, pa, ia, numbp, rflag, gname, glgth, pname, plgth, gpa);
  OutputAlignment(numbp, rflag, gpa);
  if (gpa->score <= 0.0) {				/*THE FOLLOWING PREVENTS DISPLAY OF UNSUCCESSFUL ALIGNMENTS */
    free_Pspace(glgth, plgth);
    return (0);
  }
  gpa->next = NULL;


  free_Pspace(glgth, plgth);
  return (1);

}							/* end sahmtP */



static void CompletePathMatrix(int glgth, int plgth, char *gdna, char *protein, float *pd, float *pa)
{
  int i, m, n;

/* FILL OUT THE M = 0 ZERO POSITION FOR ALL SCORING ARRAYS      */
  for (i = C_STATE; i < IC_STATE; i++)
    scoreP[i][0][0] = 0.0;
  for (n = 0; n < 3; n++) {
    scoreP[C_STATE][n][0] = 0.0;
    pathP[C_STATE][n][0] = C_N;
    scoreP[IA_STATE][n][0] = 0.0;
    pathP[IA_STATE][n][0] = IA_N;
    scoreP[IB_STATE][n][0] = 0.0;
    pathP[IB_STATE][n][0] = IB_N;
    scoreP[IC_STATE][n][0] = 0.0;
    pathP[IC_STATE][n][0] = IC_N;
  }
  for (n = 0; n < 4; ++n)
    for (m = 0; m <= plgth; ++m)
      intronstart_IA[n][m] = intronstart_IB[n][m] = intronstart_IC[n][m] = 0;


/* FILL OUT THE FIRST CODON ROWS OF THE EXON PATH MATRIX: */
  C_0m(plgth);

/* ENSURE THAT THERE CANNOT BE AN AMINO ACID WITHOUT A CODON */
  I_0m(plgth);

/* COMPLETE THE REST OF THE PATH MATRIX FOR n > 0;  */
/* STEPPING ALONG THE GENOME */
  for (n = 3; n <= glgth; n++) {
    pathP[C_STATE][n][0] = C_N;
    pathP[IA_STATE][n][0] = IA_N;
    pathP[IB_STATE][n][0] = IB_N;
    pathP[IC_STATE][n][0] = IC_N;
/* STEPPING ALONG THE PROTEIN */
    for (m = 1; m <= plgth; m++) {
      C_nm(n, m, gdna, protein, glgth, plgth, pd, pa);
      I_nm(n, m, gdna, pd, pa);
    }
  }

}							/* end CompletePathMatrix() */




static void C_0m(int plgth)
{
  int n, m;

  for (n = 0; n < 3; n++) {
    for (m = 1; m <= plgth; m++) {
      scoreP[C_STATE][n][m] = 0.0;
      pathP[C_STATE][n][m] = C_M;
    }
  }
  return;

}							/* end C_0m() */



static void I_0m(int plgth)
{
  int m, n;

  for (n = 0; n < 3; n++) {
    for (m = 1; m <= plgth; m++) {
      scoreP[IA_STATE][n][m] = -99999.;
      pathP[IA_STATE][n][m] = IA_N;
      scoreP[IB_STATE][n][m] = -99999.;
      pathP[IB_STATE][n][m] = IB_N;
      scoreP[IC_STATE][n][m] = -99999.;
      pathP[IC_STATE][n][m] = IC_N;
    }
  }
  return;

}							/* end I_0m() */





static void C_nm(int n, int m, char *gdna, char *protein, int glgth, int plgth, float *pd, float *pa)
{
  float x[11];
  int xpathP[11], largest;

  x[0] = scoreP[C_STATE][MOD4(n-3)][m - 1] +
    weight(C_STATE, C_N3M, n, glgth, m, plgth, gdna[n - 3], gdna[n - 2], gdna[n - 1], protein[m - 1], pd, pa);
  xpathP[0] = C_N3M;

  x[1] = scoreP[C_STATE][MOD4(n-2)][m - 1] +
    weight(C_STATE, C_N2M, n, glgth, m, plgth, gdna[n - 2], gdna[n - 1], DASH, protein[m - 1], pd, pa);
  xpathP[1] = C_N2M;

  x[2] = scoreP[C_STATE][MOD4(n-1)][m - 1] +
    weight(C_STATE, C_N1M, n, glgth, m, plgth, gdna[n - 1], DASH, DASH, protein[m - 1], pd, pa);
  xpathP[2] = C_N1M;

  x[3] = scoreP[C_STATE][MOD4(n)][m - 1] +
    weight(C_STATE, C_M, n, glgth, m, plgth, DASH, DASH, DASH, protein[m - 1], pd, pa);
  xpathP[3] = C_M;

  x[4] = scoreP[C_STATE][MOD4(n-3)][m] +
    weight(C_STATE, C_N3, n, glgth, m, plgth, gdna[n - 3], gdna[n - 2], gdna[n - 1], PROT_DASH, pd, pa);
  xpathP[4] = C_N3;

  x[5] = scoreP[C_STATE][MOD4(n-2)][m] +
    weight(C_STATE, C_N2, n, glgth, m, plgth, gdna[n - 2], gdna[n - 1], DASH, PROT_DASH, pd, pa);
  xpathP[5] = C_N2;

  x[6] = scoreP[C_STATE][MOD4(n-1)][m] +
    weight(C_STATE, C_N, n, glgth, m, plgth, gdna[n - 1], DASH, DASH, PROT_DASH, pd, pa);
  xpathP[6] = C_N;

  x[7] = scoreP[IA_STATE][MOD4(n-3)][m] +
    weight(C_STATE, I_N3, n, glgth, m, plgth, gdna[n - 3], gdna[n - 2], gdna[n - 1], PROT_DASH, pd, pa);
  if (n - 2 - intronstart_IA[MOD4(n-3)][m] < MIN_INTRONLENGTH)
    x[7] += SIWEIGHT;
  xpathP[7] = I_N3;

  x[8] = scoreP[IA_STATE][MOD4(n-3)][m - 1] +
    weight(C_STATE, I_N3M, n, glgth, m, plgth, gdna[n - 3], gdna[n - 2], gdna[n - 1], protein[m - 1], pd, pa);
  if (n - 2 - intronstart_IA[MOD4(n-3)][m - 1] < MIN_INTRONLENGTH)
    x[8] += SIWEIGHT;
  xpathP[8] = I_N3M;

  x[9] = scoreP[IB_STATE][MOD4(n-2)][m - 1] +
    weight(C_STATE, I_N2M, n, glgth, m, plgth,
    splitcodon_IB[MOD4(n-2)][m - 1], gdna[n - 2], gdna[n - 1], protein[m - 1], pd, pa);
  if (n - 1 - intronstart_IB[MOD4(n-2)][m - 1] < MIN_INTRONLENGTH)
    x[9] += SIWEIGHT;
  xpathP[9] = I_N2M;

  x[10] = scoreP[IC_STATE][MOD4(n-1)][m - 1] +
    weight(C_STATE, I_N1M, n, glgth, m, plgth,
    splitcodon_IC1[MOD4(n-1)][m - 1], splitcodon_IC2[MOD4(n-1)][m - 1], gdna[n - 1], protein[m - 1], pd, pa);
  if (n - intronstart_IC[MOD4(n-1)][m - 1] < MIN_INTRONLENGTH)
    x[10] += SIWEIGHT;
  xpathP[10] = I_N1M;

  largest = IndexMax(x, 11);
  scoreP[C_STATE][MOD4(n)][m] = x[largest];
  pathP[C_STATE][n][m] = (short) (xpathP[largest]);

  return;

}							/* end C_nm() */



static void I_nm(int n, int m, char *gdna, float *pd, float *pa)
{
  IA_nm(n, m, pd, pa);
  IB_nm(n, m, gdna, pd, pa);
  IC_nm(n, m, gdna, pd, pa);
  return;

}							/* end I_nm() */



static void IA_nm(int n, int m, float *pd, float *pa)
{
  float x[2];
  int xpathP[2], largest;

  x[0] = (float) (scoreP[IA_STATE][MOD4(n-1)][m] + Log_1_Pa[n-2]);
  xpathP[0] = IA_N;

  x[1] = (float) (scoreP[C_STATE][MOD4(n-1)][m] + Log_Pd[n-1]);
  xpathP[1] = C_N;

  largest = (x[0]>=x[1])?0:1;
  scoreP[IA_STATE][MOD4(n)][m] = x[largest];
  pathP[IA_STATE][n][m] = (short) (xpathP[largest]);

  switch (largest) {
  case 0:
    intronstart_IA[MOD4(n)][m] = intronstart_IA[MOD4(n-1)][m];
    break;
  case 1:
    intronstart_IA[MOD4(n)][m] = n;
    break;
  default:
    printf("\n\n  bug alert in IB_nm function! \n\n");
  }

  return;

}							/* end IA_nm() */



static void IB_nm(int n, int m, char *gdna, float *pd, float *pa)
{
  float x[2];
  int xpathP[2], largest;

  x[0] = (float) (scoreP[IB_STATE][MOD4(n-1)][m] + Log_1_Pa[n-2]);
  xpathP[0] = IB_N;

  x[1] = (float) (scoreP[C_STATE][MOD4(n-2)][m] + Log_Pd[n-1]);
  xpathP[1] = C_N2;

  largest = (x[0]>=x[1])?0:1;
  scoreP[IB_STATE][MOD4(n)][m] = x[largest];
  pathP[IB_STATE][n][m] = (short) (xpathP[largest]);

  switch (largest) {
  case 0:
    intronstart_IB[MOD4(n)][m] = intronstart_IB[MOD4(n-1)][m];
    splitcodon_IB[MOD4(n)][m] = splitcodon_IB[MOD4(n-1)][m];
    break;
  case 1:
    intronstart_IB[MOD4(n)][m] = n;
    splitcodon_IB[MOD4(n)][m] = gdna[n - 2];
    break;
  default:
    printf("\n\n  bug alert in IB_nm function! \n\n");
  }

  return;

}							/* end IB_nm() */



static void IC_nm(int n, int m, char *gdna, float *pd, float *pa)
{
  float x[2];
  int xpathP[2], largest;

  x[0] = (float) (scoreP[IC_STATE][MOD4(n-1)][m] + Log_1_Pa[n-2]);
  xpathP[0] = IC_N;

  x[1] = (float) (scoreP[C_STATE][MOD4(n-3)][m] + Log_Pd[n-1]);
  xpathP[1] = C_N3;
  largest = (x[0]>=x[1])?0:1;
  scoreP[IC_STATE][MOD4(n)][m] = x[largest];
  pathP[IC_STATE][n][m] = (short) (xpathP[largest]);

  switch (largest) {
  case 0:
    intronstart_IC[MOD4(n)][m] = intronstart_IC[MOD4(n-1)][m];
    splitcodon_IC1[MOD4(n)][m] = splitcodon_IC1[MOD4(n-1)][m];
    splitcodon_IC2[MOD4(n)][m] = splitcodon_IC2[MOD4(n-1)][m];
    break;
  case 1:
    intronstart_IC[MOD4(n)][m] = n;
    splitcodon_IC1[MOD4(n)][m] = gdna[n - 3];
    splitcodon_IC2[MOD4(n)][m] = gdna[n - 2];
    break;
  }

  return;

}							/* end IC_nm() */



static void OptimalPath(int glgth, int plgth, char *gdna, char *protein)
{

  int largest;
  float x[4];


  /*FIRST, FIND WHICH OF THE FOUR MATRICES HAS THE HIGHEST FINAL VALUE */
  x[0] = scoreP[C_STATE][MOD4(glgth)][plgth];
  x[1] = scoreP[IA_STATE][MOD4(glgth)][plgth];
  x[2] = scoreP[IB_STATE][MOD4(glgth)][plgth];
  x[3] = scoreP[IC_STATE][MOD4(glgth)][plgth];

  largest = IndexMax(x, 4);


  optN = glgth;
  optM = plgth;
  optcounter = 0;

  TracePath(largest, gdna, protein, plgth);

}							/* end OptimalPath() */



static void TracePath(int state, char *gdna, char *protein, int plgth)
{
  if (optN == 0 && optM == 0)
    return;

  switch (pathP[state][optN][optM]) {
  case (C_N3M):
    WritePath(gdna[optN - 1], gdna[optN - 2], gdna[optN - 3], protein[optM - 1], state);
    optN -= 3;
    optM--;
    TracePath(C_STATE, gdna, protein, plgth);
    return;

  case (C_N2M):
    WritePath(gdna[optN - 1], DASH, gdna[optN - 2], protein[optM - 1], state);
    optN -= 2;
    optM--;
    TracePath(C_STATE, gdna, protein, plgth);
    return;

  case (C_N1M):
    WritePath(DASH, gdna[optN - 1], DASH, protein[optM - 1], state);
    optN--;
    optM--;
    TracePath(C_STATE, gdna, protein, plgth);
    return;

  case (C_M):
    WritePath(DASH, DASH, DASH, protein[optM - 1], state);
    optM--;
    TracePath(C_STATE, gdna, protein, plgth);
    return;

  case (C_N3):
    if (state == C_STATE) {
      WritePath(gdna[optN - 1], gdna[optN - 2], gdna[optN - 3], PROT_DASH, state);
    }
    else {
      WritePath(gdna[optN - 1], DUMMY, DUMMY, STAR, state);
      if (optM == plgth)
	WritePath(DUMMY, gdna[optN - 2], gdna[optN - 3], PROT_DASH, C_STATE);
      else
	WritePath(DUMMY, gdna[optN - 2], gdna[optN - 3], protein[optM], C_STATE);
    }
    optN -= 3;
    TracePath(C_STATE, gdna, protein, plgth);
    return;

  case (C_N2):
    if (state == C_STATE) {
      WritePath(gdna[optN - 1], DASH, gdna[optN - 2], PROT_DASH, state);
    }
    else {
      WritePath(gdna[optN - 1], DUMMY, DUMMY, STAR, state);
      WritePath(DUMMY, gdna[optN - 2], DUMMY, BLANK, C_STATE);
    }
    optN -= 2;
    TracePath(C_STATE, gdna, protein, plgth);
    return;

  case (C_N):
    if (state == C_STATE) {
      WritePath(DASH, gdna[optN - 1], DASH, PROT_DASH, state);
    }
    else {
      WritePath(gdna[optN - 1], DUMMY, DUMMY, STAR, state);
    }
    optN--;
    TracePath(C_STATE, gdna, protein, plgth);
    return;

  case (I_N3):
    WritePath(gdna[optN - 1], gdna[optN - 2], gdna[optN - 3], PROT_DASH, state);
    optN -= 3;
    TracePath(IA_STATE, gdna, protein, plgth);
    return;

  case (I_N3M):
    WritePath(gdna[optN - 1], gdna[optN - 2], gdna[optN - 3], protein[optM - 1], state);
    optN -= 3;
    optM--;
    TracePath(IA_STATE, gdna, protein, plgth);
    return;

  case (I_N2M):
    WritePath(gdna[optN - 1], gdna[optN - 2], DUMMY, protein[optM - 1], state);
    optN -= 2;
    optM--;
    TracePath(IB_STATE, gdna, protein, plgth);
    return;

  case (I_N1M):
    WritePath(gdna[optN - 1], DUMMY, DUMMY, BLANK, state);
    optN--;
    optM--;
    TracePath(IC_STATE, gdna, protein, plgth);
    return;

  case (IA_N):
    WritePath(gdna[optN - 1], DUMMY, DUMMY, DUMMY, state);
    optN--;
    TracePath(IA_STATE, gdna, protein, plgth);
    return;

  case (IB_N):
    WritePath(gdna[optN - 1], DUMMY, DUMMY, DUMMY, state);
    optN--;
    TracePath(IB_STATE, gdna, protein, plgth);
    return;

  case (IC_N):
    WritePath(gdna[optN - 1], DUMMY, DUMMY, DUMMY, state);
    optN--;
    TracePath(IC_STATE, gdna, protein, plgth);
    return;
  }

}							/* end TracePath() */



static void WritePath(int g1, int g2, int g3, int aa, int state)
{
  if (state == C_STATE) {
    WriteCodon(g1, g2, g3, aa);
  }

  else {
    WriteIntron(g1, state);
  }
  return;

}							/* end WritePath() */



static void WriteCodon(int g1, int g2, int g3, int aa)
{

  /* position 1 */
  if (g1 != DUMMY) {
    optGDNA[optcounter] = g1;
    optProtein[optcounter] = BLANK;
    optState[optcounter] = C_STATE;
    optcounter++;
  }

  /* position 2 */
  if (g2 != DUMMY) {
    optGDNA[optcounter] = g2;
    optProtein[optcounter] = aa;
    optState[optcounter] = C_STATE;
    optcounter++;
  }

  /* position 0 */
  if (g3 != DUMMY) {
    optGDNA[optcounter] = g3;
    optProtein[optcounter] = BLANK;
    optState[optcounter] = C_STATE;
    optcounter++;
  }

  return;

}							/* end WriteCodon() */



static void WriteIntron(int g, int state)
{

  optGDNA[optcounter] = g;
  optProtein[optcounter] = STAR;
  optState[optcounter] = state;
  optcounter++;

  return;

}							/* end WriteIntron() */



int Obeg, Gbeg, Pbeg, Oend;

static void OutputEnds(int ia)
{
  int i;

  Gbeg = ia;
  Pbeg = 0;

  for (i = 0; i < optcounter; i++) {
    if (optGDNA[optcounter - 1 - i] < 11) {
      ++Gbeg;
      if (optProtein[optcounter - 1 - i] < 24) {
	++Pbeg;
	if (optGDNA[optcounter - 1 - i - 1] < 11) {
	  --Gbeg;
	  --i;
	  break;
	}
      }
    }
    else if (optProtein[optcounter - 1 - i] < 24)
      ++Pbeg;
  }
  Obeg = i;

  for (i = 0; i < optcounter; i++) {
    if (optGDNA[i] < 11) {
      if (optProtein[i] < 24)
	if (optGDNA[i + 1] < 11) {
	  if (optProtein[i - 1] == 31)
	    --i;
	  break;
	}
    }
  }
  Oend = optcounter - 1 - i;

}							/* end OutputEnds() */



static void OutputParse(float *pd, float *pa, int ia, int numbp,
  int rflag, char *gname, int glgth, char *pname, int plgth, struct gpalgnmnt *gpa)
{
  int i, j, exn = 0;
  int currentState, GBegin, GEnd, PBegin, PEnd;
  int mlgth;
  float escr, oescr, gescr, ogescr, cvrge, simscore;
#ifdef HTMLWS
  char anchorName[257];
#endif

  mlgth = 0;
  escr = oescr = gescr = ogescr = simscore = cvrge = 0.0;
  
  if (rflag) {
    CHECKALIGNLEN; nextFree+=sprintf(global_algnmnt+nextFree, "\nPredicted gene structure (within gDNA segment %d to %d):\n\n",numbp - gpa->gia,numbp - gpa->gib);
  }
  else {
    CHECKALIGNLEN; nextFree+=sprintf(global_algnmnt+nextFree, "\nPredicted gene structure (within gDNA segment %d to %d):\n\n",gpa->gia + 1,gpa->gib + 1);
  }

  strcpy(gpa->gname, gname);
  strcpy(gpa->pname, pname);
  currentState = optState[optcounter - 1 - Obeg];
  GBegin = GEnd = Gbeg;
  PBegin = Pbeg;
  PEnd = Pbeg - 1;
  for (i = Obeg + 1; i <= Oend; i++) {
    if (optState[optcounter - 1 - i] != currentState) {
      if (currentState == 0) {
	if (rflag) {
	  gpa->gcds[exn][0] = numbp - GBegin + 1;
	  gpa->gcds[exn][1] = numbp - GEnd + 1;
	  CHECKALIGNLEN; nextFree+=sprintf(global_algnmnt+nextFree, " Exon %2d %*d %*d ", ++exn, ifwdth, numbp - GBegin + 1, ifwdth, numbp - GEnd + 1);
	}
	else {
	  gpa->gcds[exn][0] = GBegin;
	  gpa->gcds[exn][1] = GEnd;
	  CHECKALIGNLEN; nextFree+=sprintf(global_algnmnt+nextFree, " Exon %2d %*d %*d ", ++exn, ifwdth, GBegin, ifwdth, GEnd);
	}
	CHECKALIGNLEN; nextFree+=sprintf(global_algnmnt+nextFree, "(%4d n);  Protein %5d %5d (%4d aa);", GEnd - GBegin + 1, PBegin, PEnd, PEnd - PBegin + 1);
	if (oescr > 0) {
	  CHECKALIGNLEN; nextFree+=sprintf(global_algnmnt+nextFree, " score: %5.3f\n", escr / oescr);
	}
	else {
	  CHECKALIGNLEN; nextFree+=sprintf(global_algnmnt+nextFree, " score:     0\n");
	}
	gpa->exnscr[exn - 1] = escr / oescr;
	PBegin = PEnd + 1;
	if (GEnd - GBegin + 1 >= MINELGTHfS) {
	  gescr += escr;
	  ogescr += oescr;
	}
	mlgth += GEnd - GBegin + 1;
	escr = oescr = 0.0;
      }
      else if (exn > 0) {
	if (rflag) {
	  CHECKALIGNLEN; nextFree+=sprintf(global_algnmnt+nextFree, "  Intron %2d %*d %*d ", exn, ifwdth, numbp - GBegin + 1, ifwdth, numbp - GEnd + 1);
	}
	else {
	  CHECKALIGNLEN; nextFree+=sprintf(global_algnmnt+nextFree, "  Intron %2d %*d %*d ", exn, ifwdth, GBegin, ifwdth, GEnd);
	}
	CHECKALIGNLEN; nextFree+=sprintf(global_algnmnt+nextFree, "(%4d n);  ", GEnd - GBegin + 1);
	CHECKALIGNLEN; nextFree+=sprintf(global_algnmnt+nextFree, "Pd: %5.3f   ", pd[GBegin - 1 - ia]);
	gpa->itrscr[exn - 1][0] = pd[GBegin - 1 - ia];
	CHECKALIGNLEN; nextFree+=sprintf(global_algnmnt+nextFree, "Pa: %5.3f", pa[GEnd - 1 - ia]);
	gpa->itrscr[exn - 1][1] = pa[GEnd - 1 - ia];
	if (GEnd - GBegin + 1 <= MIN_INTRONLENGTH) {
	  sprintf(global_algnmnt+nextFree, " ??");
	}
	AppendChar('\n');
      }
      GBegin = GEnd + 1;
    }
/*SCORING: */
    if (optProtein[optcounter - 1 - i] < 23) {
      oescr += SCALEF * ssmat[optProtein[optcounter - 1 - i]][optProtein[optcounter - 1 - i]];
      if (optState[optcounter - i] == 0) {
	if (optState[optcounter - i - 2] == 0) {
	  escr += GetScore(optGDNA[optcounter - i],
	    optGDNA[optcounter - i - 1], optGDNA[optcounter - i - 2], optProtein[optcounter - 1 - i]);
	}
	else {						/* state 3 */
	  for (j = 1;; ++j) {
	    if (optState[optcounter - i - 2 - j] == 0) {
	      escr += GetScore(optGDNA[optcounter - i],
		optGDNA[optcounter - i - 1], optGDNA[optcounter - i - 2 - j], optProtein[optcounter - 1 - i]);
	      break;
	    }
	  }
	}
      }
      else {						/* state 2 */
	for (j = 1;; ++j) {
	  if (optState[optcounter - i + j] == 0) {
	    escr += GetScore(optGDNA[optcounter - i + j],
	      optGDNA[optcounter - i - 1], optGDNA[optcounter - i - 2], optProtein[optcounter - 1 - i]);
	    break;
	  }
	}
      }
    }
    else if (optProtein[optcounter - 1 - i] == 24)
      oescr += (-SCALEF * (float) (INDEL_PEN));
    currentState = optState[optcounter - 1 - i];
    if (optGDNA[optcounter - 1 - i] < 11)
      ++GEnd;
    if (optProtein[optcounter - 1 - i] < 23)
      ++PEnd;
  }
  if (currentState == 0) {
    if (rflag) {
      gpa->gcds[exn][0] = numbp - GBegin + 1;
      gpa->gcds[exn][1] = numbp - GEnd + 1;
      CHECKALIGNLEN; nextFree+=sprintf(global_algnmnt+nextFree, " Exon %2d %*d %*d ", ++exn, ifwdth, numbp - GBegin + 1, ifwdth, numbp - GEnd + 1);
    }
    else {
      gpa->gcds[exn][0] = GBegin;
      gpa->gcds[exn][1] = GEnd;
      CHECKALIGNLEN; nextFree+=sprintf(global_algnmnt+nextFree, " Exon %2d %*d %*d ", ++exn, ifwdth, GBegin, ifwdth, GEnd);
    }
    CHECKALIGNLEN; nextFree+=sprintf(global_algnmnt+nextFree,"(%4d n);  Protein %5d %5d (%4d aa);", GEnd - GBegin + 1, PBegin, PEnd, PEnd - PBegin + 1);
    if (oescr > 0) {
      CHECKALIGNLEN; nextFree+=sprintf(global_algnmnt+nextFree, " score: %5.3f\n", escr / oescr);
    }
    else {
      CHECKALIGNLEN; nextFree+=sprintf(global_algnmnt+nextFree, " score:     0\n");
    }
    gpa->exnscr[exn - 1] = escr / oescr;
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
    CHECKALIGNLEN; nextFree+=sprintf(global_algnmnt+nextFree, "(%4d n);  ", GEnd - GBegin + 1);
    CHECKALIGNLEN; nextFree+=sprintf(global_algnmnt+nextFree, "Pd: %5.3f   ", pd[GBegin - 1 - ia]);
    gpa->itrscr[exn - 1][0] = pd[GBegin - 1 - ia];
    CHECKALIGNLEN; nextFree+=sprintf(global_algnmnt+nextFree, "Pa: %5.3f", pa[GEnd - 1 - ia]);
    gpa->itrscr[exn - 1][1] = pa[GEnd - 1 - ia];
    if (GEnd - GBegin + 1 <= MIN_INTRONLENGTH) {
      CHECKALIGNLEN; nextFree+=sprintf(global_algnmnt+nextFree, " ??");
    }
    CHECKALIGNLEN; nextFree+=sprintf(global_algnmnt+nextFree, "\n");
  }

  cvrge = (float) mlgth / (float) glgth > (float) mlgth / (float) (3 * plgth) ?
    (float) mlgth / (float) glgth : (float) mlgth / (float) (3 * plgth);
  if (ogescr > 0.0  ||  cvrge >= MINCVRGEfS) {
    if (ogescr > 0.0) simscore = gescr / ogescr;
    CHECKALIGNLEN; nextFree+=sprintf(global_algnmnt+nextFree, "\nMATCH\t%s\t%s\t%5.3f\t%d\t%5.3f", gname, pname, simscore, mlgth, cvrge);
    if ((float) mlgth / (float) glgth > (float) mlgth / (float) (3 * plgth)) {
      CHECKALIGNLEN; nextFree+= sprintf(global_algnmnt+nextFree, "\tG\n");
    }
    else {
      CHECKALIGNLEN; nextFree+=sprintf(global_algnmnt+nextFree, "\tP\n");
    }

    CHECKALIGNLEN; nextFree+=sprintf(global_algnmnt+nextFree, "PGS_%s_%s\t(", gname, pname);
    for (i = 0; i < exn - 1; ++i) {
      CHECKALIGNLEN; nextFree+=sprintf(global_algnmnt+nextFree, "%d  %d,", gpa->gcds[i][0], gpa->gcds[i][1]);
      
    }
    CHECKALIGNLEN; nextFree+=sprintf(global_algnmnt+nextFree, "%d  %d)\n", gpa->gcds[exn - 1][0], gpa->gcds[exn - 1][1]);

#ifdef HTMLWS
    if (htmlop && simscore > 0.0) {	/* ... simscore condition - to prevent display of unsuccessful alignments */
      sprintf(anchorName,"PPGS%d@%s",gpa->calln,gpa->gsgmntn);
      fprintf(imageDataFh,"PPGS %*d %*d %s\n",ifwdth, gpa->gcds[0][0], ifwdth, gpa->gcds[exn-1][1], anchorName);
      for (i = 0; i < exn - 1; i++) {
	fprintf(imageDataFh, "%d  %d  ",gpa->gcds[i][0],gpa->gcds[i][1]);
      }
      fprintf(imageDataFh, "%d  %d\n", gpa->gcds[i][0],gpa->gcds[i][1]);
    }
#endif
  }
  else {
    CHECKALIGNLEN; nextFree+=sprintf(global_algnmnt+nextFree, "\nNOMATCH\t%s\t%s\t%5.3f\t%5.3f\t%d\n", gname, pname, gescr, ogescr, mlgth);
  }

  gpa->exn = exn;
  gpa->score = simscore;

}							/* end OutputParse() */



static void OutputAlignment(int numbp, int rflag, struct gpalgnmnt *gpa)
{
  int i, j, k;
  int base = 0, aacid, gcount = Gbeg - 1, acount = Pbeg - 1;

  CHECKALIGNLEN; nextFree+=sprintf(global_algnmnt+nextFree, "\nAlignment:\n\n");

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
/*PRINT OUT CODON TRANSLATION: */
      aacid = 0;
      for (j = i - base + 1; j <= i; ++j) {
	if (optProtein[optcounter - 1 - j] < 25) {
	  if (optGDNA[optcounter - j - 1] <= 3) {
	    if (optState[optcounter - j] == 0) {
	      if (optGDNA[optcounter - j] <= 3) {
		if (optState[optcounter - j - 2] == 0) {
		  if (optGDNA[optcounter - j - 2] <= 3)
		    OutputResidue(cod2aa[optGDNA[optcounter - j]]
		      [optGDNA[optcounter - j - 1]][optGDNA[optcounter - j - 2]]);
		  else
		    AppendChar(' ');
		}
		else {
		  for (k = 1;; ++k) {
		    if (optState[optcounter - j - 2 - k] == 0) {
		      if (optGDNA[optcounter - j - 2 - k] <= 3)
			OutputResidue(cod2aa[optGDNA[optcounter - j]]
			  [optGDNA[optcounter - j - 1]][optGDNA[optcounter - j - 2 - k]]);
		      else
			AppendChar(' ');
		      break;
		    }
		  }
		}
	      }
	      else
		AppendChar(' ');
	    }
	    else {					/*state 2 */
	      if (optGDNA[optcounter - j - 2] > 3)
		        AppendChar(' ');
	      else {
		    for (k = 1;; ++k) {
		        if (optState[optcounter - j + k] == 0) {
		            if (optGDNA[optcounter - j + k] <= 3)
		                 OutputResidue(cod2aa[optGDNA[optcounter - j + k]]
			                [optGDNA[optcounter - j - 1]][optGDNA[optcounter - j - 2]]);
		            else
		                AppendChar(' ');
		            break;
		        }
		    }
	      }
	    }
	  }
	  else
	    AppendChar(' ');
	}
	else
	  AppendChar(' ');
	aacid++;
	if (aacid % 10 == 0)
	  AppendChar(' ');
      }
      AppendChar('\n');
/*PRINT OUT SUBSTITUTION SCORE SYMBOLS: */
      aacid = 0;
      for (j = i - base + 1; j <= i; ++j) {
	if (optProtein[optcounter - 1 - j] < 23) {
	  if (optGDNA[optcounter - j - 1] <= 3) {
	    if (optState[optcounter - j] == 0) {
	      if (optGDNA[optcounter - j] <= 3) {
		if (optState[optcounter - j - 2] == 0) {
		  if (optGDNA[optcounter - j - 2] <= 3) {
		    if (cod2aa[optGDNA[optcounter - j]]
		      [optGDNA[optcounter - j - 1]]
		      [optGDNA[optcounter - j - 2]] == optProtein[optcounter - 1 - j])
		      AppendChar('|');
		    else if (cod2aa[optGDNA[optcounter - j]]
		      [optGDNA[optcounter - j - 1]][optGDNA[optcounter - j - 2]] == 23)
		      AppendChar('*');
		    else if (ssmat[cod2aa[optGDNA[optcounter - j]]
			  [optGDNA[optcounter - j - 1]]
			[optGDNA[optcounter - j - 2]]][optProtein[optcounter - 1 - j]] > 0)
		      AppendChar('+');
		    else if (ssmat[cod2aa[optGDNA[optcounter - j]]
			  [optGDNA[optcounter - j - 1]]
			[optGDNA[optcounter - j - 2]]][optProtein[optcounter - 1 - j]] == 0)
		     AppendChar('.');
		    else
		      AppendChar(' ');
		  }
		  else
		    AppendChar(' ');
		}
		else {
		  for (k = 1;; ++k) {
		    if (optState[optcounter - j - 2 - k] == 0) {
		      if (optGDNA[optcounter - j - 2 - k] <= 3) {
			if (cod2aa[optGDNA[optcounter - j]]
			  [optGDNA[optcounter - j - 1]]
			  [optGDNA[optcounter - j - 2 - k]] == optProtein[optcounter - 1 - j])
			  AppendChar('|');
			else if (cod2aa[optGDNA[optcounter - j]]
			  [optGDNA[optcounter - j - 1]][optGDNA[optcounter - j - 2 - k]] == 23)
			  AppendChar('*');
			else if (ssmat[cod2aa[optGDNA[optcounter - j]]
			      [optGDNA[optcounter - j - 1]]
			    [optGDNA[optcounter - j - 2 - k]]][optProtein[optcounter - 1 - j]] > 0)
			  AppendChar('+');
			else if (ssmat[cod2aa[optGDNA[optcounter - j]]
			      [optGDNA[optcounter - j - 1]]
			    [optGDNA[optcounter - j - 2 - k]]][optProtein[optcounter - 1 - j]] == 0)
			  AppendChar('.');
			else
			  AppendChar(' ');
		      }
		      else
			AppendChar(' ');
		      break;
		    }
		  }
		}
	      }
	      else
		AppendChar(' ');
	    }
	    else {					/*state 2 */
	      if (optGDNA[optcounter - j - 2] > 3)
		AppendChar(' ');
	      else {
		for (k = 1;; ++k) {
		  if (optState[optcounter - j + k] == 0) {
		    if (optGDNA[optcounter - j + k] <= 3) {
		      if (cod2aa[optGDNA[optcounter - j + k]]
			[optGDNA[optcounter - j - 1]]
			[optGDNA[optcounter - j - 2]] == optProtein[optcounter - 1 - j])
			AppendChar('|');
		      else if (cod2aa[optGDNA[optcounter - j + k]]
			[optGDNA[optcounter - j - 1]][optGDNA[optcounter - j - 2]] == 23)
			AppendChar('*');
		      else if (ssmat[cod2aa[optGDNA[optcounter - j + k]]
			    [optGDNA[optcounter - j - 1]]
			  [optGDNA[optcounter - j - 2]]][optProtein[optcounter - 1 - j]] > 0)
			AppendChar('+');
		      else if (ssmat[cod2aa[optGDNA[optcounter - j + k]]
			    [optGDNA[optcounter - j - 1]]
			  [optGDNA[optcounter - j - 2]]][optProtein[optcounter - 1 - j]] == 0)
			AppendChar('.');
		      else
			AppendChar(' ');
		    }
		    else
		      AppendChar(' ');
		    break;
		  }
		}
	      }
	    }
	  }
	  else
	    AppendChar(' ');
	}
	else
	 AppendChar(' ');
	aacid++;
	if (aacid % 10 == 0)
	 AppendChar(' ');
      }
      AppendChar('\n');
/*PRINT OUT ALIGNED TARGET PROTEIN SEQUENCE: */
      aacid = 0;
      for (j = i - base + 1; j <= i; ++j) {
	OutputResidue(optProtein[optcounter - 1 - j]);
	if (optProtein[optcounter - 1 - j] < 24)
	  acount++;
	aacid++;
	if (aacid % 10 == 0)
	  AppendChar(' ');
      }
      if (j == Oend + 1 && aacid % 10 != 0)
	AppendChar(' ');
      CHECKALIGNLEN; nextFree+=sprintf(global_algnmnt+nextFree, "  %*d\n\n\n", ifwdth, acount);
      base = 0;
    }
    else if (base % 10 == 0)
      AppendChar(' ');
  }
  AppendChar('\0');
  if ((gpa->algnmnt = (char *) calloc(nextFree+1, sizeof(char))) == NULL)
    fatal_error("Error: memory allocation failed. Exit.\n");
  strcpy(gpa->algnmnt, global_algnmnt);

}							/* end OutputAlignment() */



static char num2nu[]="TCAGYRSWKMN-.????????????????? ";

static void OutputBase(int base)
{

    if (base>=0 && base <=30)
      AppendChar(num2nu[base]);
    else {
      CHECKALIGNLEN; nextFree+=sprintf(global_algnmnt+nextFree, "\n\n  Unknown symbol in OutputBase function! \n\n");
    }

}							/* end OutputBase() */



static char num2aa[]="LAGSVKETDIRPNFQYHMCWBZX*-?????. ";

static void OutputResidue(int aa)
{

   if (aa>=0 && aa <=31)
     AppendChar(num2aa[aa]);
    else {
      CHECKALIGNLEN; nextFree+=sprintf(global_algnmnt+nextFree, "\n\n  Unknown symbol in OutputResidue function! \n\n");
    }

}							/* end OutputResidue() */



static float GetScore(int g1, int g2, int g3, int aa)
{

/* 1) INDEL_PEN for deletions from and insertions into genomic DNA of lengths
   1, 2, or 3, irrespective of indel size:
 */
  if (g1 == DASH || g2 == DASH || g3 == DASH || aa == PROT_DASH) {
    return (SCALEF * (float) (INDEL_PEN));
  }
  else if (g1 <= 3 && g2 <= 3 && g3 <= 3) {
/* 2) (-)/2*INDEL_PEN for matching/mismatching a stop codon:
 */
    if (cod2aa[g1][g2][g3] == 23) {
      if (aa == 23)
	return (SCALEF * (float) (-2 * INDEL_PEN));
      else
	return (SCALEF * (float) (2 * INDEL_PEN));
    }
    else
/* 3) amino acid substitution score:
 */
    {
      if (aa < 23)
	return (SCALEF * (float) ssmat[cod2aa[g1][g2][g3]][aa]);
      else
	return (SCALEF * (float) (2 * INDEL_PEN));
    }
  }
  else
/* 4) neutral score in case of wild-card characters in the genomic DNA:
 */
    return (0.0);

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



static float weight(int state, int path, int n, int glgth, int m, int plgth,
  int g1, int g2, int g3, int aa, float *pd, float *pa)
{
  float rval = 0.0;

  if (state == C_STATE) {
    /* Transition weights: */

    if (path == C_N3M) {
      rval = (float) (rval + Log_1_Pd[n-3]);
    }
    else if (path == C_N2M) {
      if (n < glgth || m < 20)
	rval = (float) (rval + Log_1_Pd[n-2]);
    }
    else if (path == C_N1M) {
      if (n < glgth || m < 20)
	rval = (float) (rval + Log_1_Pd[n-1]);
    }
    else if (path == C_M) {
      if (n < glgth || m < 20)
	rval = (float) (rval + Log_1_Pd[n]);
    }
    else if (path == C_N3) {
      if (m < plgth || n < 60)
	rval = (float) (rval + Log_1_Pd[n-3]);
    }
    else if (path == C_N2) {
      if (m < plgth || n < 60)
	rval = (float) (rval + Log_1_Pd[n-2]);
    }
    else if (path == C_N) {
      if (m < plgth || n < 60)
	rval = (float) (rval + Log_1_Pd[n-1]);
    }
    else if (path == I_N3  ||  path == I_N3M) {
      if (n > 3)
	rval = (float) (rval + Log_Pa[n-4]);
    }
    else if (path == I_N2M) {
      rval = (float) (rval + Log_Pa[n-3]);
    }
    else if (path == I_N1M) {
      rval = (float) (rval + Log_Pa[n-2]);
    }

    /* Output weights: */

    if (path == C_N2M || path == C_N1M || path == C_M) {
      if (n < glgth || m < 20)
	rval += GetScore(g1, g2, g3, aa);
    }
    else if (path == C_N3 || path == C_N2 || path == C_N) {
      if (m < plgth || n < 60)
	rval += GetScore(g1, g2, g3, aa);
    }
    else
      rval += GetScore(g1, g2, g3, aa);
  }

  else {						/* state == Ix_STATE */
    {
      if (path == IA_N || path == IB_N || path == IC_N)
	rval = (float) (rval + Log_1_Pa[n-2]);
      else
	rval = (float) (rval + Log_Pd[n-1]);
    }
  }

  return (rval);

}							/* end weight() */



static void allocate_Pspace(int glgth, int plgth)
{
  int i, j;

  for (i = 0; i < 4; ++i) {
    if ((intronstart_IA[i] = (int *) calloc((plgth + 1) , sizeof(int))) == NULL)
       fatal_error("Error: memory allocation failed. Exit.\n");
    if ((intronstart_IB[i] = (int *) calloc((plgth + 1) , sizeof(int))) == NULL)
       fatal_error("Error: memory allocation failed. Exit.\n");
    if ((intronstart_IC[i] = (int *) calloc((plgth + 1) , sizeof(int))) == NULL)
       fatal_error("Error: memory allocation failed. Exit.\n");

    if ((splitcodon_IB[i] = (int *) calloc((plgth + 1) , sizeof(int))) == NULL)
       fatal_error("Error: memory allocation failed. Exit.\n");
    if ((splitcodon_IC1[i] = (int *) calloc((plgth + 1) , sizeof(int))) == NULL)
       fatal_error("Error: memory allocation failed. Exit.\n");
    if ((splitcodon_IC2[i] = (int *) calloc((plgth + 1) , sizeof(int))) == NULL)
       fatal_error("Error: memory allocation failed. Exit.\n");

    if ((pathP[i] = (short int **) calloc((glgth + 1) , sizeof(short int *))) == NULL)
       fatal_error("Error: memory allocation failed. Exit.\n");
    for (j = 0; j < glgth + 1; ++j)
      if ((pathP[i][j] = (short int *) calloc((plgth + 1) , sizeof(short int))) == NULL)
	 fatal_error("Error: memory allocation failed. Exit.\n");

    for (j = 0; j < 4; ++j)
      if ((scoreP[i][j] = (float *) calloc((plgth + 1) , sizeof(float))) == NULL)
	 fatal_error("Error: memory allocation failed. Exit.\n");
  }

  if ((optGDNA = (int *) calloc((3 * (plgth + 1 + glgth)) , sizeof(int))) == NULL)
     fatal_error("Error: memory allocation failed. Exit.\n");
  if ((optProtein = (int *) calloc((3 * (plgth + 1 + glgth)) , sizeof(int))) == NULL)
     fatal_error("Error: memory allocation failed. Exit.\n");
  if ((optState = (int *) calloc((3 * (plgth + 1 + glgth)) , sizeof(int))) == NULL)
     fatal_error("Error: memory allocation failed. Exit.\n");

}							/* end allocate_Pspace() */



static void free_Pspace(int glgth, int plgth)
{
  int i, j;

  for (i = 0; i < 4; ++i) {
    free(intronstart_IA[i]);
    free(intronstart_IB[i]);
    free(intronstart_IC[i]);

    free(splitcodon_IB[i]);
    free(splitcodon_IC1[i]);
    free(splitcodon_IC2[i]);

    for (j = 0; j < glgth + 1; ++j)
      free(pathP[i][j]);
    free(pathP[i]);

    for (j = 0; j < 4; ++j)
      free(scoreP[i][j]);
  }
  free(optGDNA);
  free(optProtein);
  free(optState);

}							/* end free_Pspace() */



int insert_gpa(gpahpp, nodep)
  struct gpalgnmnt **gpahpp, *nodep;
{
int c;
struct gpalgnmnt *tmp;

  if (*gpahpp == NULL) {
	*gpahpp = nodep;
	nodep->next = NULL;
  }
  else {
	while ((*gpahpp)->next != NULL)
	 {c = compare_gpas(*gpahpp,nodep);
	  if (c == 0) return(0);
	  if (c == 1) {
	    nodep->next = (*gpahpp)->next;
	    tmp = *gpahpp;
	    *gpahpp= nodep;
	    free_gpa(tmp);
	    return(1);
	  }
	  if (c == 2)
            gpahpp = &((*gpahpp)->next);
	  else
	    break;
         }
        if ((*gpahpp)->next == NULL) {
        	c = compare_gpas(*gpahpp,nodep);
		if (c == 0) return(0);
		if (c == 1) {
		  nodep->next = NULL;
		  tmp = *gpahpp;
		  *gpahpp= nodep;
		  free_gpa(tmp);
		  return(1);
		}
		if (c == 2) {(*gpahpp)->next = nodep; nodep->next = NULL;}
		else
		  {nodep->next = *gpahpp; *gpahpp = nodep;}
          }
        else
          { nodep->next = *gpahpp; *gpahpp = nodep; }
  }
  return (1);

}							/* end insert_gpa() */



int compare_gpas(gpa1p, gpa2p)
  struct gpalgnmnt *gpa1p, *gpa2p;
{
	/* return 0: discard gpa2p           */
	/* return 1: replace gpa1p by gpa2p  */
	/* return 2: gpa1p "<" gpa2p         */
	/* return 3: gpa1p ">" gpa2p         */

if ( strcmp(gpa1p->gname, gpa2p->gname) == 0   &&
     strcmp(gpa1p->pname, gpa2p->pname) == 0     )
 {if ( MIN(gpa1p->gcds[0][0],gpa1p->gcds[gpa1p->exn - 1][1]) <=
       MIN(gpa2p->gcds[0][0],gpa2p->gcds[gpa2p->exn - 1][1])    &&
       MAX(gpa1p->gcds[0][0],gpa1p->gcds[gpa1p->exn - 1][1]) >=
       MAX(gpa2p->gcds[0][0],gpa2p->gcds[gpa2p->exn - 1][1])      ) return(0);
  if ( MIN(gpa2p->gcds[0][0],gpa2p->gcds[gpa2p->exn - 1][1]) <=
       MIN(gpa1p->gcds[0][0],gpa1p->gcds[gpa1p->exn - 1][1])    &&
       MAX(gpa2p->gcds[0][0],gpa2p->gcds[gpa2p->exn - 1][1]) >=
       MAX(gpa1p->gcds[0][0],gpa1p->gcds[gpa1p->exn - 1][1])      ) return(1);
 }

if ( (MIN(gpa1p->gcds[0][0],gpa1p->gcds[gpa1p->exn - 1][1]) <
      MIN(gpa2p->gcds[0][0],gpa2p->gcds[gpa2p->exn - 1][1])       )  ||
     (MIN(gpa1p->gcds[0][0],gpa1p->gcds[gpa1p->exn - 1][1]) ==
      MIN(gpa2p->gcds[0][0],gpa2p->gcds[gpa2p->exn - 1][1])     &&
      MAX(gpa1p->gcds[0][0],gpa1p->gcds[gpa1p->exn - 1][1]) >=
      MAX(gpa2p->gcds[0][0],gpa2p->gcds[gpa2p->exn - 1][1])       )     )  return(2);

return(3);

}							/* end compare_gpas() */


void free_gpa(struct gpalgnmnt *gpa)
{
  if (gpa->algnmnt!=NULL)
    free((char *)gpa->algnmnt);
  free((struct gpalgnmnt *) gpa);
}							/* end free_gpa() */



int prt_gpa_list(FILE *fp, struct gpalgnmnt *gpa, char *gsgmntn)
{
  int ngpa = 0;
  char anchorName[513];

  while (gpa != NULL) {
    if (htmlop) {
      sprintf(anchorName,"PPGS%d@%s",gpa->calln,gpa->gsgmntn);
      ADD_NAME(fp,anchorName);
    }
    fprintf(fp, "%s", gpa->algnmnt);
    ++ngpa;
    gpa = gpa->next;
  }

  
  return (ngpa);

}							/* end prt_gpa_list() */
