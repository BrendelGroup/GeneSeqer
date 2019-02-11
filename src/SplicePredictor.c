/* SplicePredictor.c;                         Last update: February 11, 2019. */
/* Dependencies:   getlns.c getgbs.c sahmtD.c sahmtP.c                        */
/* Bugs:  add: ATG; STP; assembly check                                       */

/* Corresponding author:                                                      */

/*   Volker Brendel, Department of Biology                                    */
/*   Indiana University, Bloomington, IN 47405                                */
/*   (812) 855-7074, vbrendel@indiana.edu                                     */

/* Past contributing authors:                                                 */
/*   Liqun Xing, Department of Zoology & Genetics, Iowa State University      */
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



#ifdef DAPBM
#define USAGE1 "Usage:\
 %s [-s species] [-c cutoff] [-t pval] [-T sval]\n\
    [-oO] [-n topN] [-u] [-U] [-p pstyle] [-w [nsites sscwl]] [-x] [-eE estdbn]\n\
    [-i prmfile] [-qQ qpfname] [-I matname] [-a from] [-b to] [-rR]\n\
    [-lL libfname] [-g gbfname(s)]\n\n\
  -s species :   Set species to select the most appropriate splice site models.\n\
                  This parameter must be specified. Options:\n\
                  'human', 'mouse', 'rat', 'chicken', 'Drosophila', 'nematode'\n\
                  'yeast', 'Aspergillus', 'Arabidopsis', 'maize', 'rice',\n\
                  'Medicago', 'Populus', 'generic'.\n\
  -c cutoff   :  set log (BF) cutoff\n\
  -t pval     :  Set prediction threshold to 'pval' [overrides -c option]\n\
  -T sval     :  Set prediction threshold to 'sval' [overrides -c option]\n\
  -o          :  Order sites by P-value [default: by position]\n\
  -O          :  Order sites by *-value [default: by position]\n\
  -n topN     :  Display top N splice sites\n\
  -u          :  Disable local pruning of non-optimal sites\n\
  -U          :  Score also non-canonical splice site dinucleotides\n\
  -p pstyle   :  1 (terse=WWW); 2 (default); 3 (very terse=EXDOMINO);\n\
                 4 (verbose); 5 (spreadsheet)\n\
  -w [nsites sscwl]:  report splice site clusters\n\
                     (>= nsites in <= sscwl bases; default: 4/1500,\n\
                      appropriate for -T 14 option)\n\
  -x [from to]: LaTex graphical output in *.tex file(s)\n\
  -eE estdbn :   Read EST sequence data from library file 'estdbn'\n\
                  (FASTA-format: >name).\n\
		  -e: align + strand only, -E: align + and - strands.\n\
  -i prmfile :   Read parameters for EST alignments from file 'prmfile'.\n"
#define USAGE2 "\
  -qQ qpfname:   Read target protein sequence data from library file\n\
                  'qpfname' (FASTA-format; '-q': >gi|name; '-Q': >name).\n\
  -I matname :   Read amino acid substitution scoring matrix from file\n\
                  'matname' [default: PBLO62].\n\
  -a from    :   Analyze genomic sequence from position 'from'.\n\
  -b to      :   Analyze genomic sequence up to position 'to'.\n\
  -r         :   Analyze reverse strand.\n\
  -R         :   Analyze both strands.\n\
  -lL libfname : Read (multiple) sequence data from library file 'libfname'\n\
                  (FASTA-format; '-l': >gi|name; '-L': >name).\n\
  -g gbfname(s): Read nucleic acid sequence data from GenBank file(s)\n\
                  'gbfname(s)'.  If specified, the -g option MUST BE LAST.\n\n"
#else
#define USAGE1 "Usage:\
 %s [-m model] [-s species] [-c cutoff] [-t pval] [-T sval]\n\
    [-oO] [-n topN] [-p pstyle] [-w [nsites sscwl]] [-x] [-eE estdbn]\n\
    [-i prmfile] [-qQ qpfname] [-I matname] [-a from] [-b to] [-rR]\n\
    [-lL libfname] [-g gbfname(s)]\n\n\
  -m model:  set model [0= without/ 1= with (default) sub-classification]\n\
  -s species:  set species [options: 'maize' (default); 'Arabidopsis']\n\
  -c level:  set prediction threshold level\n\
              [0 = all GU (AG) sites with 50 base flanks\n\
               1 = threshold at 100%% sensitivity for training set [default]\n\
               2 = threshold at  95%% sensitivity for training set\n\
               3 = threshold at maximal Tau value for training set\n\
              ]\n\
  -t pval:  set prediction threshold to 'pval' [overrides -c option]\n\
  -T sval:  set prediction threshold to 'sval' [overrides -c option]\n\
  -o     :  order sites by P-value [default: by position]\n\
  -O     :  order sites by *-value [default: by position]\n\
  -n topN:  display top N splice sites\n\
  -p pstyle:  1 (terse=WWW); 2 (default); 3 (very terse=EXDOMINO);\n\
              4 (verbose); 5 (spreadsheet)\n\
  -w [nsites sscwl]:  report splice site clusters\n\
                     (>= nsites in <= sscwl bases; default: 4/1500,\n\
                      appropriate for -T 14 option)\n\
  -x [from to]: LaTex graphical output in *.tex file(s)\n\
  -eE estdb  :   Read EST sequence data from library file 'estdbn';\n\
		  -e: align + strand only, -E: align + and - strands.\n\
  -i prmfile :   Read parameters for EST alignments from file 'prmfile'.\n"
#define USAGE2 "\
  -qQ qpfname:   Read target protein sequence data from library file\n\
                  'qpfname' (FASTA-format; '-q': >gi|name; '-Q': >name).\n\
  -I matname :   Read amino acid substitution scoring matrix from file\n\
                  'matname' [default: PBLO62].\n\
  -a from    :   Analyze genomic sequence from position 'from'.\n\
  -b to      :   Analyze genomic sequence up to position 'to'.\n\
  -r         :   Analyze reverse strand.\n\
  -R         :   Analyze both strands.\n\
  -lL libfname : Read (multiple) sequence data from library file 'libfname'\n\
                  (FASTA-format; '-l': >gi|name; '-L': >name).\n\
  -g gbfname(s): Read nucleic acid sequence data from GenBank file(s)\n\
                  'gbfname(s)'.  If specified, the -g option MUST BE LAST.\n\n"
#endif


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include "space.h"
#include "smat.h"
#include "def.h"
#include "abc.h"
#include "sahmt.h"
#include "platform.h"
#include "RawTextFile.h"	/* Header, raw text file services */
#include "html.h"

#define CHECKFLAG 0

#define MAXNGENES 500
#define MAXADL     50
#define MAXRANGE  101
int ADLGTH= 15;

int pstyle= 2, estop= 0, qpop= 0, qnoffset = 4, prmop= 0, fsflag= 0;
int DELTAL= 0, DELTAR= 0, lflag;
int obsset= 0, mpsflag, pvalset= 0, svalset= 0, prtclstrs= 0;
int lxflag= 0, dispset= 0, rflag= 0, bflag= 0;
int TOPN= 0;
int pruning= 1, allownn= 0;
int nsites= 4, sscwl= 1500;
double pval;
int sval;
char libname[257], estdbn[257], qpfname[257], prmfile[257], sfname[257],
	lxfname[257], ftstrg[257], dstrg[257];
FILE *infp, *estfp, *qpfp, *prmfp, *matfp, *latexfp;
FILE *outfp;

#define INPUT_FILE   1
#define EST_FILE     2
#define QP_FILE      3

int mflag= 0, minsc, maxsc, smat[23][23];
int WIDTH= 400, HEIGHT= 100, VLGTH= 4;
float SCALEF;
int gbfcount= 0, fcount= 0;
char *gdna, *gdnaR;
int af[23], pf[9];
float aq[20][4], pq[9][4];

float *logV[8]={NULL};
char  *global_algnmnt=NULL;
int   logPia=0; /* to pass info */

struct sprmtr sprm;

struct bcvct {
  int bcnt[11], ccnt[6];
  float bfrq[11], cfrq[6];
 } tbcv;

struct splsite {
  int stype, loc, ADstrg[MAXADL], glgth, istr;
  double scr, cval, Lscr, delU, delS;
  double rho, gamma;
  int pstar, rstar, gstar, sstar, vote;
  struct splsite *nexts, *prevs;
 };

double pscr[4];
int MAXRWstrg[4][MAXADL];

#include "sequence.h"

#ifdef DAPBM
#ifdef DAPBM7
#include "daPbm7.h"
#else
#include "daPbm2.h"
#endif
#else
#include "daPll.h"
#endif

struct eicrta {
  int mincolgth, minielgth, minmelgth, mintelgth, minilgth, maxilgth;
  float minapop, mindscr, minascr;
 } gcrta;

int DISPFR, DISPTO;

int SGMNTSZE = 300000;
int SGMNTSHFT =   102;
int fromp, frompA, frompp, top, htmlop, ifwdth;




void reverse_seq(char *seq, int numbp)
{
  int i, tmpn;

  for (i = numbp - 1; i >= (numbp + 1) / 2; --i) {
    if (seq[i] <= 3)
      tmpn = (seq[i] + 2) % 4;
    else
      tmpn = 10;
    if (seq[numbp - 1 - i] <= 3)
      seq[i] = (char) ((seq[numbp - 1 - i] + 2) % 4);
    else
      seq[i] = 10;
    seq[numbp - 1 - i] = (char) (tmpn);
  }
  if (numbp % 2 == 1) {
    if (seq[numbp / 2] <= 3)
      seq[numbp / 2] = (char) ((seq[numbp / 2] + 2) % 4);
    else
      seq[numbp / 2] = 10;
  }

}				/* end reverse_seq() */



void complement_seq(char *seq, int numbp, char *seqR)
{
  int i;

  for (i = numbp - 1; i >= 0; --i) {
    if (seq[i] <= 3)
      seqR[numbp - 1 - i] = (char) ((seq[i] + 2) % 4);
    else
      seqR[numbp - 1 - i] = 10;
  }

}				/* end complement_seq() */



void reverse_ftc(int ftc[MAXNBAF][2], int nbaf, int numbp)
{
  int i;

  for (i = 0; i < nbaf; ++i) {
    ftc[i][0] = numbp - ftc[i][0] + 1;
    ftc[i][1] = numbp - ftc[i][1] + 1;
  }

}				/* end reverse_ftc() */



int base_composition(char *seq, int ia, int ib, struct bcvct *bcv)
{
  int i, n;

  for (i = 0; i < 11; ++i)
    bcv->bcnt[i] = 0;
  for (i = 0; i < 6; ++i)
    bcv->ccnt[i] = 0;
  n = ib - ia + 1;

  for (i = ia; i <= ib; ++i)
    ++bcv->bcnt[(int)seq[i]];

  bcv->ccnt[0] = bcv->bcnt[2] + bcv->bcnt[3];	/* frequency of R= A or G */
  bcv->ccnt[1] = bcv->bcnt[0] + bcv->bcnt[1];	/* frequency of Y= C or T */
  bcv->ccnt[2] = bcv->bcnt[1] + bcv->bcnt[3];	/* frequency of S= C or G */
  bcv->ccnt[3] = bcv->bcnt[0] + bcv->bcnt[2];	/* frequency of W= A or T */
  bcv->ccnt[4] = bcv->bcnt[0] + bcv->bcnt[3];	/* frequency of K= G or T */
  bcv->ccnt[5] = bcv->bcnt[1] + bcv->bcnt[2];	/* frequency of M= A or C */

  for (i = 0; i < 11; ++i)
    bcv->bfrq[i] = (float) bcv->bcnt[i] / n;
  for (i = 0; i < 6; ++i)
    bcv->cfrq[i] = (float) bcv->ccnt[i] / n;

  if (!bcv->bcnt[10])
    return (1);
  else
    return (0);

}				/* end base_composition() */



void prt_base_composition(FILE *fp, struct bcvct *bcv)
{

  fprintf(fp, "\n");
  fprintf(fp, "T:%6d(%4.1f%%); C:%6d(%4.1f%%); A:%6d(%4.1f%%); G:%6d(%4.1f%%)",
	  bcv->bcnt[0], 100. * bcv->bfrq[0], bcv->bcnt[1], 100. * bcv->bfrq[1],
	  bcv->bcnt[2], 100. * bcv->bfrq[2], bcv->bcnt[3], 100. * bcv->bfrq[3]);
  if (bcv->bcnt[10])
    fprintf(fp, "; N:%3d(%4.1f%%)", bcv->bcnt[10], 100. * bcv->bfrq[10]);

  fprintf(fp, "\nR:%6d(%4.1f%%); Y:%6d(%4.1f%%);",
	  bcv->ccnt[0], 100. * bcv->cfrq[0], bcv->ccnt[1], 100. * bcv->cfrq[1]);
  fprintf(fp, " S:%6d(%4.1f%%); W:%6d(%4.1f%%)\n",
	  bcv->ccnt[2], 100. * bcv->cfrq[2], bcv->ccnt[3], 100. * bcv->cfrq[3]);
  fprintf(fp, "K:%6d(%4.1f%%); M:%6d(%4.1f%%)\n",
	  bcv->ccnt[4], 100. * bcv->cfrq[4], bcv->ccnt[5], 100. * bcv->cfrq[5]);

}				/* end prt_base_composition() */



void word_composition(FILE *fp, char *seq, int crd[MAXNBAF][2], int ia, int ib,
		      int pflag, int cflag, int tflag, float MINgcPCNT, float MAXgcPCNT)
{
  int i, j, k, l, m, n, o, phs = 0, ctflag;
  int m_pcnt[3][4], m_pcntS[3], m_ccnt[4], m_ccntS = 0;
  int d_pcnt[3][4][4], d_pcntS[3], d_ccnt[4][4], d_ccntS = 0;
  int t_pcnt[3][4][4][4], t_pcntS[3], t_ccnt[4][4][4], t_ccntS = 0;
  int r_pcnt[3][4][4][4][4], r_pcntS[3], r_ccnt[4][4][4][4], r_ccntS = 0;
  int p_pcnt[3][4][4][4][4][4], p_pcntS[3], p_ccnt[4][4][4][4][4], p_ccntS = 0;
  int h_pcnt[3][4][4][4][4][4][4], h_pcntS[3], h_ccnt[4][4][4][4][4][4], h_ccntS = 0;
  struct bcvct bcv;

  for (i = 0; i < 3; ++i) {
    m_pcntS[i] = d_pcntS[i] = t_pcntS[i] = r_pcntS[i] = p_pcntS[i] = h_pcntS[i] = 0;
    for (j = 0; j < 4; ++j) {
      m_pcnt[i][j] = 0;
      for (k = 0; k < 4; ++k) {
	d_pcnt[i][j][k] = 0;
	for (l = 0; l < 4; ++l) {
	  t_pcnt[i][j][k][l] = 0;
	  for (m = 0; m < 4; ++m) {
	    r_pcnt[i][j][k][l][m] = 0;
	    for (n = 0; n < 4; ++n) {
	      p_pcnt[i][j][k][l][m][n] = 0;
	      for (o = 0; o < 4; ++o) {
		h_pcnt[i][j][k][l][m][n][o] = 0;
	      }
	    }
	  }
	}
      }
    }
  }
  for (i = 0; i < 4; ++i) {
    m_ccnt[i] = 0;
    for (j = 0; j < 4; ++j) {
      d_ccnt[i][j] = 0;
      for (k = 0; k < 4; ++k) {
	t_ccnt[i][j][k] = 0;
	for (l = 0; l < 4; ++l) {
	  r_ccnt[i][j][k][l] = 0;
	  for (m = 0; m < 4; ++m) {
	    p_ccnt[i][j][k][l][m] = 0;
	    for (n = 0; n < 4; ++n) {
	      h_ccnt[i][j][k][l][m][n] = 0;
	    }
	  }
	}
      }
    }
  }
  ctflag = 0;
  for (i = ia; i <= ib; ++i) {
    if (tflag == 0)
      ctflag = 0;
    else {
      base_composition(seq, crd[i][0] - 1, crd[i][1] - 1, &bcv);
      if (100. * bcv.cfrq[2] >= MINgcPCNT && 100. * bcv.cfrq[2] <= MAXgcPCNT)
	ctflag = 1;
      else
	ctflag = 0;
    }
    if (!cflag)
      phs = 0;
    for (j = crd[i][0] - 1; j < crd[i][1] - 5; ++j) {
      phs = (phs + 1) % 3;
      if (seq[j] != 10) {
	++m_pcntS[phs];
	++m_pcnt[phs][(int)seq[j]];
	if ((int)seq[j + 1] != 10) {
	  ++d_pcntS[phs];
	  ++d_pcnt[phs][(int)seq[j]][(int)seq[j + 1]];
	  if ((int)seq[j + 2] != 10) {
	    ++t_pcntS[phs];
	    ++t_pcnt[phs][(int)seq[j]][(int)seq[j + 1]][(int)seq[j + 2]];
	    if ((int)seq[j + 3] != 10) {
	      ++r_pcntS[phs];
	      ++r_pcnt[phs][(int)seq[j]][(int)seq[j + 1]][(int)seq[j + 2]][(int)seq[j + 3]];
	      if ((int)seq[j + 4] != 10) {
		++p_pcntS[phs];
		++p_pcnt[phs][(int)seq[j]][(int)seq[j + 1]][(int)seq[j + 2]][(int)seq[j + 3]][(int)seq[j + 4]];
		if ((int)seq[j + 5] != 10) {
		  ++h_pcntS[phs];
		  ++h_pcnt[phs][(int)seq[j]][(int)seq[j + 1]][(int)seq[j + 2]][(int)seq[j + 3]][(int)seq[j + 4]]
		    [(int)seq[j + 5]];
		}
	      }
	    }
	  }
	}
      }
    }
    if (i < ib) {
      if (crd[i][1] - crd[i][0] > 3) {
	phs = (phs + 1) % 3;
	if ((int)seq[crd[i][1] - 5] != 10) {
	  ++m_pcntS[phs];
	  ++m_pcnt[phs][(int)seq[crd[i][1] - 5]];
	  if ((int)seq[crd[i][1] - 4] != 10) {
	    ++d_pcntS[phs];
	    ++d_pcnt[phs][(int)seq[crd[i][1] - 5]][(int)seq[crd[i][1] - 4]];
	    if ((int)seq[crd[i][1] - 3] != 10) {
	      ++t_pcntS[phs];
	      ++t_pcnt[phs][(int)seq[crd[i][1] - 5]][(int)seq[crd[i][1] - 4]][(int)seq[crd[i][1] - 3]];
	      if ((int)seq[crd[i][1] - 2] != 10) {
		++r_pcntS[phs];
		++r_pcnt[phs][(int)seq[crd[i][1] - 5]][(int)seq[crd[i][1] - 4]]
		  [(int)seq[crd[i][1] - 3]][(int)seq[crd[i][1] - 2]];
		if ((int)seq[crd[i][1] - 1] != 10) {
		  ++p_pcntS[phs];
		  ++p_pcnt[phs][(int)seq[crd[i][1] - 5]][(int)seq[crd[i][1] - 4]]
		    [(int)seq[crd[i][1] - 3]][(int)seq[crd[i][1] - 2]][(int)seq[crd[i][1] - 1]];
		  if (cflag && (int)seq[crd[i + 1][0] - 1] != 10) {
		    ++h_pcntS[phs];
		    ++h_pcnt[phs][(int)seq[crd[i][1] - 5]][(int)seq[crd[i][1] - 4]]
		      [(int)seq[crd[i][1] - 3]][(int)seq[crd[i][1] - 2]][(int)seq[crd[i][1] - 1]]
		      [(int)seq[crd[i + 1][0] - 1]];
		  }
		}
	      }
	    }
	  }
	}
      }
      if (crd[i][1] - crd[i][0] > 2) {
	phs = (phs + 1) % 3;
	if ((int)seq[crd[i][1] - 4] != 10) {
	  ++m_pcntS[phs];
	  ++m_pcnt[phs][(int)seq[crd[i][1] - 4]];
	  if ((int)seq[crd[i][1] - 3] != 10) {
	    ++d_pcntS[phs];
	    ++d_pcnt[phs][(int)seq[crd[i][1] - 4]][(int)seq[crd[i][1] - 3]];
	    if ((int)seq[crd[i][1] - 2] != 10) {
	      ++t_pcntS[phs];
	      ++t_pcnt[phs][(int)seq[crd[i][1] - 4]][(int)seq[crd[i][1] - 3]][(int)seq[crd[i][1] - 2]];
	      if ((int)seq[crd[i][1] - 1] != 10) {
		++r_pcntS[phs];
		++r_pcnt[phs][(int)seq[crd[i][1] - 4]][(int)seq[crd[i][1] - 3]]
		  [(int)seq[crd[i][1] - 2]][(int)seq[crd[i][1] - 1]];
		if (cflag && (int)seq[crd[i + 1][0] - 1] != 10) {
		  ++p_pcntS[phs];
		  ++p_pcnt[phs][(int)seq[crd[i][1] - 4]][(int)seq[crd[i][1] - 3]]
		    [(int)seq[crd[i][1] - 2]][(int)seq[crd[i][1] - 1]]
		    [(int)seq[crd[i + 1][0] - 1]];
		  if (cflag && crd[i + 1][1] - crd[i + 1][0] > 0 && (int)seq[crd[i + 1][0]] != 10) {
		    ++h_pcntS[phs];
		    ++h_pcnt[phs][(int)seq[crd[i][1] - 4]][(int)seq[crd[i][1] - 3]]
		      [(int)seq[crd[i][1] - 2]][(int)seq[crd[i][1] - 1]][(int)seq[crd[i + 1][0] - 1]]
		      [(int)seq[crd[i + 1][0]]];
		  }
		}
	      }
	    }
	  }
	}
      }
      if (crd[i][1] - crd[i][0] > 1) {
	phs = (phs + 1) % 3;
	if ((int)seq[crd[i][1] - 3] != 10) {
	  ++m_pcntS[phs];
	  ++m_pcnt[phs][(int)seq[crd[i][1] - 3]];
	  if ((int)seq[crd[i][1] - 2] != 10) {
	    ++d_pcntS[phs];
	    ++d_pcnt[phs][(int)seq[crd[i][1] - 3]][(int)seq[crd[i][1] - 2]];
	    if ((int)seq[crd[i][1] - 1] != 10) {
	      ++t_pcntS[phs];
	      ++t_pcnt[phs][(int)seq[crd[i][1] - 3]][(int)seq[crd[i][1] - 2]][(int)seq[crd[i][1] - 1]];
	      if (cflag && (int)seq[crd[i + 1][0] - 1] != 10) {
		++r_pcntS[phs];
		++r_pcnt[phs][(int)seq[crd[i][1] - 3]][(int)seq[crd[i][1] - 2]]
		  [(int)seq[crd[i][1] - 1]][(int)seq[crd[i + 1][0] - 1]];
		if (cflag && crd[i + 1][1] - crd[i + 1][0] > 0 && (int)seq[crd[i + 1][0]] != 10) {
		  ++p_pcntS[phs];
		  ++p_pcnt[phs][(int)seq[crd[i][1] - 3]][(int)seq[crd[i][1] - 2]]
		    [(int)seq[crd[i][1] - 1]][(int)seq[crd[i + 1][0] - 1]]
		    [(int)seq[crd[i + 1][0]]];
		  if (cflag && crd[i + 1][1] - crd[i + 1][0] > 1 && (int)seq[crd[i + 1][0] + 1] != 10) {
		    ++h_pcntS[phs];
		    ++h_pcnt[phs][(int)seq[crd[i][1] - 3]][(int)seq[crd[i][1] - 2]]
		      [(int)seq[crd[i][1] - 1]][(int)seq[crd[i + 1][0] - 1]][(int)seq[crd[i + 1][0]]]
		      [(int)seq[crd[i + 1][0] + 1]];
		  }
		}
	      }
	    }
	  }
	}
      }
      if (crd[i][1] - crd[i][0] > 0) {
	phs = (phs + 1) % 3;
	if ((int)seq[crd[i][1] - 2] != 10) {
	  ++m_pcntS[phs];
	  ++m_pcnt[phs][(int)seq[crd[i][1] - 2]];
	  if ((int)seq[crd[i][1] - 1] != 10) {
	    ++d_pcntS[phs];
	    ++d_pcnt[phs][(int)seq[crd[i][1] - 2]][(int)seq[crd[i][1] - 1]];
	    if (cflag && (int)seq[crd[i + 1][0] - 1] != 10) {
	      ++t_pcntS[phs];
	      ++t_pcnt[phs][(int)seq[crd[i][1] - 2]][(int)seq[crd[i][1] - 1]]
		[(int)seq[crd[i + 1][0] - 1]];
	      if (cflag && crd[i + 1][1] - crd[i + 1][0] > 0 && (int)seq[crd[i + 1][0]] != 10) {
		++r_pcntS[phs];
		++r_pcnt[phs][(int)seq[crd[i][1] - 2]][(int)seq[crd[i][1] - 1]]
		  [(int)seq[crd[i + 1][0] - 1]][(int)seq[crd[i + 1][0]]];
		if (cflag && crd[i + 1][1] - crd[i + 1][0] > 1 && (int)seq[crd[i + 1][0] + 1] != 10) {
		  ++p_pcntS[phs];
		  ++p_pcnt[phs][(int)seq[crd[i][1] - 2]][(int)seq[crd[i][1] - 1]]
		    [(int)seq[crd[i + 1][0] - 1]][(int)seq[crd[i + 1][0]]]
		    [(int)seq[crd[i + 1][0] + 1]];
		  if (cflag && crd[i + 1][1] - crd[i + 1][0] > 2 && (int)seq[crd[i + 1][0] + 2] != 10) {
		    ++h_pcntS[phs];
		    ++h_pcnt[phs][(int)seq[crd[i][1] - 2]][(int)seq[crd[i][1] - 1]]
		      [(int)seq[crd[i + 1][0] - 1]][(int)seq[crd[i + 1][0]]]
		      [(int)seq[crd[i + 1][0] + 1]][(int)seq[crd[i + 1][0] + 2]];
		  }
		}
	      }
	    }
	  }
	}
      }
      phs = (phs + 1) % 3;
      if ((int)seq[crd[i][1] - 1] != 10) {
	++m_pcntS[phs];
	++m_pcnt[phs][(int)seq[crd[i][1] - 1]];
	if (cflag && (int)seq[crd[i + 1][0] - 1] != 10) {
	  ++d_pcntS[phs];
	  ++d_pcnt[phs][(int)seq[crd[i][1] - 1]][(int)seq[crd[i + 1][0] - 1]];
	  if (cflag && crd[i + 1][1] - crd[i + 1][0] > 0 && (int)seq[crd[i + 1][0]] != 10) {
	    ++t_pcntS[phs];
	    ++t_pcnt[phs][(int)seq[crd[i][1] - 1]][(int)seq[crd[i + 1][0] - 1]][(int)seq[crd[i + 1][0]]];
	    if (cflag && crd[i + 1][1] - crd[i + 1][0] > 1 && (int)seq[crd[i + 1][0] + 1] != 10) {
	      ++r_pcntS[phs];
	      ++r_pcnt[phs][(int)seq[crd[i][1] - 1]][(int)seq[crd[i + 1][0] - 1]]
		[(int)seq[crd[i + 1][0]]][(int)seq[crd[i + 1][0] + 1]];
	      if (cflag && crd[i + 1][1] - crd[i + 1][0] > 2 && (int)seq[crd[i + 1][0] + 2] != 10) {
		++p_pcntS[phs];
		++p_pcnt[phs][(int)seq[crd[i][1] - 1]][(int)seq[crd[i + 1][0] - 1]]
		  [(int)seq[crd[i + 1][0]]][(int)seq[crd[i + 1][0] + 1]]
		  [(int)seq[crd[i + 1][0] + 2]];
		if (cflag && crd[i + 1][1] - crd[i + 1][0] > 3 && (int)seq[crd[i + 1][0] + 3] != 10) {
		  ++h_pcntS[phs];
		  ++h_pcnt[phs][(int)seq[crd[i][1] - 1]][(int)seq[crd[i + 1][0] - 1]]
		    [(int)seq[crd[i + 1][0]]][(int)seq[crd[i + 1][0] + 1]]
		    [(int)seq[crd[i + 1][0] + 2]][(int)seq[crd[i + 1][0] + 3]];
		}
	      }
	    }
	  }
	}
      }
    }
    else {
      phs = (phs + 1) % 3;
      if (crd[i][1] - crd[i][0] > 3 && (int)seq[crd[i][1] - 5] != 10) {
	++m_pcntS[phs];
	++m_pcnt[phs][(int)seq[crd[i][1] - 5]];
	if ((int)seq[crd[i][1] - 4] != 10) {
	  ++d_pcntS[phs];
	  ++d_pcnt[phs][(int)seq[crd[i][1] - 5]][(int)seq[crd[i][1] - 4]];
	  if ((int)seq[crd[i][1] - 3] != 10) {
	    ++t_pcntS[phs];
	    ++t_pcnt[phs][(int)seq[crd[i][1] - 5]][(int)seq[crd[i][1] - 4]][(int)seq[crd[i][1] - 3]];
	    if ((int)seq[crd[i][1] - 2] != 10) {
	      ++r_pcntS[phs];
	      ++r_pcnt[phs][(int)seq[crd[i][1] - 5]][(int)seq[crd[i][1] - 4]][(int)seq[crd[i][1] - 3]]
		[(int)seq[crd[i][1] - 2]];
	      if ((int)seq[crd[i][1] - 1] != 10) {
		++p_pcntS[phs];
		++p_pcnt[phs][(int)seq[crd[i][1] - 5]][(int)seq[crd[i][1] - 4]]
		  [(int)seq[crd[i][1] - 3]][(int)seq[crd[i][1] - 2]][(int)seq[crd[i][1] - 1]];
	      }
	    }
	  }
	}
      }
      phs = (phs + 1) % 3;
      if (crd[i][1] - crd[i][0] > 2 && (int)seq[crd[i][1] - 4] != 10) {
	++m_pcntS[phs];
	++m_pcnt[phs][(int)seq[crd[i][1] - 4]];
	if ((int)seq[crd[i][1] - 3] != 10) {
	  ++d_pcntS[phs];
	  ++d_pcnt[phs][(int)seq[crd[i][1] - 4]][(int)seq[crd[i][1] - 3]];
	  if ((int)seq[crd[i][1] - 2] != 10) {
	    ++t_pcntS[phs];
	    ++t_pcnt[phs][(int)seq[crd[i][1] - 4]][(int)seq[crd[i][1] - 3]][(int)seq[crd[i][1] - 2]];
	    if ((int)seq[crd[i][1] - 1] != 10) {
	      ++r_pcntS[phs];
	      ++r_pcnt[phs][(int)seq[crd[i][1] - 4]][(int)seq[crd[i][1] - 3]][(int)seq[crd[i][1] - 2]]
		[(int)seq[crd[i][1] - 1]];
	    }
	  }
	}
      }
      phs = (phs + 1) % 3;
      if (crd[i][1] - crd[i][0] > 1 && (int)seq[crd[i][1] - 3] != 10) {
	++m_pcntS[phs];
	++m_pcnt[phs][(int)seq[crd[i][1] - 3]];
	if ((int)seq[crd[i][1] - 2] != 10) {
	  ++d_pcntS[phs];
	  ++d_pcnt[phs][(int)seq[crd[i][1] - 3]][(int)seq[crd[i][1] - 2]];
	  if ((int)seq[crd[i][1] - 1] != 10) {
	    ++t_pcntS[phs];
	    ++t_pcnt[phs][(int)seq[crd[i][1] - 3]][(int)seq[crd[i][1] - 2]][(int)seq[crd[i][1] - 1]];
	  }
	}
      }
      phs = (phs + 1) % 3;
      if (crd[i][1] - crd[i][0] > 0 && (int)seq[crd[i][1] - 2] != 10) {
	++m_pcntS[phs];
	++m_pcnt[phs][(int)seq[crd[i][1] - 2]];
	if ((int)seq[crd[i][1] - 1] != 10) {
	  ++d_pcntS[phs];
	  ++d_pcnt[phs][(int)seq[crd[i][1] - 2]][(int)seq[crd[i][1] - 1]];
	}
      }
      phs = (phs + 1) % 3;
      if ((int)seq[crd[i][1] - 1] != 10) {
	++m_pcntS[phs];
	++m_pcnt[phs][(int)seq[crd[i][1] - 1]];
      }
    }
  }				/* end for (i=ia ... */

  if (pflag) {
    for (i = 1; i < 4; ++i) {
      if (pstyle == 4) {
	fprintf(fp, "\nPhase %1d:\n", i % 3);
	prt_word_composition(fp, m_pcnt[i % 3], m_pcntS[i % 3], d_pcnt[i % 3], d_pcntS[i % 3],
			     t_pcnt[i % 3], t_pcntS[i % 3], r_pcnt[i % 3], r_pcntS[i % 3], p_pcnt[i % 3],
			     p_pcntS[i % 3], h_pcnt[i % 3], h_pcntS[i % 3]);
      }
      if (ctflag) {
	for (j = 0; j < 4; ++j) {
	  gm_pcnt[i % 3][j] += m_pcnt[i % 3][j];
	  for (k = 0; k < 4; ++k) {
	    gd_pcnt[i % 3][j][k] += d_pcnt[i % 3][j][k];
	    for (l = 0; l < 4; ++l) {
	      gt_pcnt[i % 3][j][k][l] += t_pcnt[i % 3][j][k][l];
	      for (m = 0; m < 4; ++m) {
		gr_pcnt[i % 3][j][k][l][m] += r_pcnt[i % 3][j][k][l][m];
		for (n = 0; n < 4; ++n) {
		  gp_pcnt[i % 3][j][k][l][m][n] += p_pcnt[i % 3][j][k][l][m][n];
		  for (o = 0; o < 4; ++o) {
		    gh_pcnt[i % 3][j][k][l][m][n][o] += h_pcnt[i % 3][j][k][l][m][n][o];
		  }
		}
	      }
	    }
	  }
	}
	gm_pcntS[i % 3] += m_pcntS[i % 3];
	gd_pcntS[i % 3] += d_pcntS[i % 3];
	gt_pcntS[i % 3] += t_pcntS[i % 3];
	gr_pcntS[i % 3] += r_pcntS[i % 3];
	gp_pcntS[i % 3] += p_pcntS[i % 3];
	gh_pcntS[i % 3] += h_pcntS[i % 3];
      }
    }
  }
  else {
    for (i = 0; i < 4; ++i) {
      for (j = 0; j < 3; ++j)
	m_ccnt[i] += m_pcnt[j][i];
      if (ctflag)
	gm_ccnt[i] += m_ccnt[i];
      for (k = 0; k < 4; ++k) {
	for (j = 0; j < 3; ++j)
	  d_ccnt[i][k] += d_pcnt[j][i][k];
	if (ctflag)
	  gd_ccnt[i][k] += d_ccnt[i][k];
	for (l = 0; l < 4; ++l) {
	  for (j = 0; j < 3; ++j)
	    t_ccnt[i][k][l] += t_pcnt[j][i][k][l];
	  if (ctflag)
	    gt_ccnt[i][k][l] += t_ccnt[i][k][l];
	  for (m = 0; m < 4; ++m) {
	    for (j = 0; j < 3; ++j)
	      r_ccnt[i][k][l][m] += r_pcnt[j][i][k][l][m];
	    if (ctflag)
	      gr_ccnt[i][k][l][m] += r_ccnt[i][k][l][m];
	    for (n = 0; n < 4; ++n) {
	      for (j = 0; j < 3; ++j)
		p_ccnt[i][k][l][m][n] += p_pcnt[j][i][k][l][m][n];
	      if (ctflag)
		gp_ccnt[i][k][l][m][n] += p_ccnt[i][k][l][m][n];
	      for (o = 0; o < 4; ++o) {
		for (j = 0; j < 3; ++j)
		  h_ccnt[i][k][l][m][n][o] += h_pcnt[j][i][k][l][m][n][o];
		if (ctflag)
		  gh_ccnt[i][k][l][m][n][o] += h_ccnt[i][k][l][m][n][o];
	      }
	    }
	  }
	}
      }
    }
    for (i = 0; i < 3; ++i) {
      m_ccntS += m_pcntS[i];
      d_ccntS += d_pcntS[i];
      t_ccntS += t_pcntS[i];
      r_ccntS += r_pcntS[i];
      p_ccntS += p_pcntS[i];
      h_ccntS += h_pcntS[i];
    }
    if (pstyle == 4)
      prt_word_composition(fp, m_ccnt, m_ccntS, d_ccnt, d_ccntS, t_ccnt, t_ccntS,
			   r_ccnt, r_ccntS, p_ccnt, p_ccntS, h_ccnt, h_ccntS);
    if (ctflag) {
      gm_ccntS += m_ccntS;
      gd_ccntS += d_ccntS;
      gt_ccntS += t_ccntS;
      gr_ccntS += r_ccntS;
      gp_ccntS += p_ccntS;
      gh_ccntS += h_ccntS;
    }
  }

}				/* end word_composition() */



void prt_word_composition(FILE *fp, int m_cnt[4], int nm, int d_cnt[4][4],
			  int nd, int t_cnt[4][4][4], int nt, int r_cnt[4][4][4][4], int nr,
			  int p_cnt[4][4][4][4][4], int np, int h_cnt[4][4][4][4][4][4], int nh)
{
  int i, j, k, l, m, n;

  fprintf(fp,
	  "prefix  |           T                C                A                G\n");
  fprintf(fp, "        |   ");
  for (i = 0; i < 4; ++i)
    fprintf(fp, "%6d(%9.5f)", m_cnt[i], (float) m_cnt[i] / nm);
  fprintf(fp, "\n\n");

  fprintf(fp,
	  "prefix  |           T                C                A                G\n");
  for (i = 0; i < 4; ++i) {
    fprintf(fp, "%c       |   ", NAUC[i]);
    for (j = 0; j < 4; ++j)
      fprintf(fp, "%6d(%9.5f)", d_cnt[i][j], (float) d_cnt[i][j] / nd);
    fprintf(fp, "\n");
  }
  fprintf(fp, "\n");

  fprintf(fp,
	  "prefix  |           T                C                A                G\n");
  for (i = 0; i < 4; ++i)
    for (j = 0; j < 4; ++j) {
      fprintf(fp, "%c%c      |   ", NAUC[i], NAUC[j]);
      for (k = 0; k < 4; ++k)
	fprintf(fp, "%6d(%9.5f)", t_cnt[i][j][k], (float) t_cnt[i][j][k] / nt);
      fprintf(fp, "\n");
    }
  fprintf(fp, "\n");

  fprintf(fp,
	  "prefix  |           T                C                A                G\n");
  for (i = 0; i < 4; ++i)
    for (j = 0; j < 4; ++j)
      for (k = 0; k < 4; ++k) {
	fprintf(fp, "%c%c%c   |   ", NAUC[i], NAUC[j], NAUC[k]);
	for (l = 0; l < 4; ++l)
	  fprintf(fp, "%6d(%9.5f)", r_cnt[i][j][k][l], (float) r_cnt[i][j][k][l] / nr);
	fprintf(fp, "\n");
      }
  fprintf(fp, "\n");

  fprintf(fp,
	  "prefix  |           T                C                A                G\n");
  for (i = 0; i < 4; ++i)
    for (j = 0; j < 4; ++j)
      for (k = 0; k < 4; ++k)
	for (l = 0; l < 4; ++l) {
	  fprintf(fp, "%c%c%c%c  |   ", NAUC[i], NAUC[j], NAUC[k], NAUC[l]);
	  for (m = 0; m < 4; ++m)
	    fprintf(fp, "%6d(%9.5f)", p_cnt[i][j][k][l][m],
		    (float) p_cnt[i][j][k][l][m] / np);
	  fprintf(fp, "\n");
	}
  fprintf(fp, "\n");

  fprintf(fp,
	  "prefix  |           T                C                A                G\n");
  for (i = 0; i < 4; ++i)
    for (j = 0; j < 4; ++j)
      for (k = 0; k < 4; ++k)
	for (l = 0; l < 4; ++l)
	  for (m = 0; m < 4; ++m) {
	    fprintf(fp, "%c%c%c%c%c |   ", NAUC[i], NAUC[j], NAUC[k], NAUC[l], NAUC[m]);
	    for (n = 0; n < 4; ++n)
	      fprintf(fp, "%6d(%9.5f)", h_cnt[i][j][k][l][m][n],
		      (float) h_cnt[i][j][k][l][m][n] / nh);
	    fprintf(fp, "\n");
	  }
  fprintf(fp, "\n");

}				/* end prt_word_composition() */



void prt_seq_segment(FILE *fp, char *seq, int numbp, int ia, int ib, int rflag,
		     char ABC[], int nflag)
	/* prints specified segment of seq[] in ABC translation */
{
  int i, ni = 0;

  for (i = ia; i <= ib + ni; ++i) {
    if (seq[i] == 24)
      ++ni;
    if ((i - ni - ia) % 10 == 0 && (i - ni - ia) % 60 != 0 && seq[i] != 24)
      fprintf(fp, " ");
    if ((i - ni - ia) % 60 == 0 && seq[i] != 24) {
      if (nflag) {
	if (rflag == 1)
	  fprintf(fp, "\n%*d  ", ifwdth, numbp - i + ni);
	else
	  fprintf(fp, "\n%*d  ", ifwdth, i - ni + 1);
      }
      else
	fprintf(fp, "\n");
    }
    fprintf(fp, "%c", ABC[(int)seq[i]]);
  }
  fprintf(fp, "\n");

}				/* end prt_seq_segment() */



void prt_seq_segment_to_str(char *str, char *seq, int numbp, int ia, int ib,
			    int rflag, char ABC[], int nflag)
	/* prints specified segment of seq[] in ABC translation */
{
  int i, ni = 0;

  for (i = ia; i <= ib + ni; ++i) {
    if (seq[i] == 24)
      ++ni;
    if ((i - ni - ia) % 10 == 0 && (i - ni - ia) % 60 != 0 && seq[i] != 24)
      sprintf(&(str[strlen(str)]), " ");
    if ((i - ni - ia) % 60 == 0 && seq[i] != 24) {
      if (nflag) {
	if (rflag == 1)
	  sprintf(&(str[strlen(str)]), "\n%*d  ", ifwdth, numbp - i + ni);
	else
	  sprintf(&(str[strlen(str)]), "\n%*d  ", ifwdth, i - ni + 1);
      }
      else
	sprintf(&(str[strlen(str)]), "\n");
    }
    sprintf(&(str[strlen(str)]), "%c", ABC[(int)seq[i]]);
  }
  sprintf(&(str[strlen(str)]), "\n");

}				/* end prt_seq_segment_to_str() */



void prt_segment_trl(FILE *fp, int seq[], int numbp, int ia, int ib, int rflag,
		     int phsA, int cl, char ABC[])
	/* prints specified segment of seq[] in ABC translation */
{
  int i, rc, iaa = ia;;

  rc = 0;
  if (phsA == 1) {
    rc = (cl + 3) / 3;
    fprintf(fp, "\n%*d  ", ifwdth, rc);
  }
  if (phsA == 2) {
    ia += 2;
    rc = (cl + 5) / 3;
    fprintf(fp, "\n%*d    ", ifwdth, rc);
  }
  if (phsA == 0) {
    ia += 1;
    rc = (cl + 4) / 3;
    fprintf(fp, "\n%*d   ", ifwdth, rc);
  }
  for (i = ia; i <= ib; ++i) {
    if ((i - iaa) % 10 == 0 && (i - iaa) % 60 != 0)
      fprintf(fp, " ");
    if (i > iaa && (i - iaa) % 60 == 0)
      fprintf(fp, "\n%*d  ", ifwdth, rc);
    if ((i - 1 - ia) % 3 == 0 && i + 1 <= ib) {
      ++rc;
      if (seq[i - 1] > 3 || seq[i] > 3 || seq[i + 1] > 3)
	fprintf(fp, "%c", ABC[22]);	/* N containing codons are printed as X */
      else
	fprintf(fp, "%c", ABC[codtoaa[seq[i - 1]][seq[i]][seq[i + 1]]]);
    }
    else
      fprintf(fp, " ");
  }
  fprintf(fp, "\n");

}				/* end prt_segment_trl() */



int trl_seq_segment(int dseq[], int ia, int ib, int phsA, int pseq[], int imrkr[],
		    int nintrns)
	/* translate specified segment of dseq[] in specified phase */
{
  int i, j = 0, ni = 0;

  if (ia <= ib) {
    if (phsA == 2)
      ia += 2;
    if (phsA == 0)
      ia += 1;
    for (i = ia, j = 0; i + 2 <= ib; i += 3, ++j) {
      if (dseq[i] > 3 || dseq[i + 1] > 3 || dseq[i + 2] > 3)
	pseq[j] = 22;
      else
	pseq[j] = codtoaa[dseq[i]][dseq[i + 1]][dseq[i + 2]];
      if (ni < nintrns && j == imrkr[ni] + ni) {
	pseq[++j] = 24;
	++ni;
      }
    }
  }
  else {
    if (phsA == 2)
      ib -= 2;
    if (phsA == 0)
      ib -= 1;
    for (i = ib, j = 0; i - 2 >= ia; i -= 3, ++j) {
      if (dseq[i] > 3 || dseq[i - 1] > 3 || dseq[i - 2] > 3)
	pseq[j] = 22;
      else
	pseq[j] = codtoaa[(dseq[i] + 2) % 4][(dseq[i - 1] + 2) % 4][(dseq[i - 2] + 2) % 4];
      if (ni < nintrns && j == imrkr[ni]) {
	pseq[++j] = 24;
	++ni;
      }
    }
  }

  return (j - ni);

}				/* end trl_seq_segment() */



int translate(int dseq[], int ia, int ib, int pseq[], int numaa)
{
  int id = ia, ip;

  if (ia <= ib) {
    for (ip = 0; ip < numaa; ++ip) {
      if (dseq[id] > 3 || dseq[id + 1] > 3 || dseq[id + 2] > 3) {
	fprintf(stdout,
		"\nWARNING: Sequence contains ambiguous characters. Not used.\n");
	return (0);
      }
      pseq[ip] = codtoaa[dseq[id]][dseq[id + 1]][dseq[id + 2]];
      id += 3;
    }
  }
  else {
    for (ip = 0; ip < numaa; ++ip) {
      pseq[ip] = codtoaa[(dseq[id] + 2) % 4][(dseq[id - 1] + 2) % 4][(dseq[id - 2] + 2) % 4];
      id += 3;
    }
  }
  return (1);

}				/* end translate() */



int getlps(int Selection, char detsz, char noffset, char *sfname, char *seq);
int getlns(int Selection, char detsz, char noffset, char *sfname, char *seq);
int getgbs(int Selection, char detsz, char *sfname, char *seq);

void fatal_error(char *buf);

void doit(int numbp, int ia, int ib, int rflag);
#ifdef DAPBM
void find_splsignalsBM(char seq[], int numbp, int rflag, int ia, int ib, struct splsite **ss_headpp);
#else
void find_splsignalsLL(char seq[], int numbp, int rflag, int ia, int ib, struct splsite **ss_headpp);
#endif
void makedoublylinked(struct splsite *nodep);
void local_pruning(struct splsite *nodep, struct splsite **ss_headpp);
void rm_splsite(struct splsite *nodep, struct splsite **ss_headpp);
void local_optimality(char seq[], int numbp, int rflag, int ia, int ib, struct splsite *nodep, struct eicrta *crta);
int re_write(int *ADstrg, int *ADloc, double *ADscr, int *RWstrg, int i, int lgth);
void re_write_next(int *ADstrg, int *ADloc, double *ADscr, int *RWstrg, int i, int lgth);
void scr_RWstrg(int *RWstrg, int lgth, int *ADstrg, int *ADloc, double *ADscr);
void order_by_scr(struct splsite **nodepp);
void insert_by_scr(struct splsite **headpp, struct splsite *nodep);
void order_by_sstar(struct splsite **nodepp);
void insert_by_sstar(struct splsite **headpp, struct splsite *nodep);
void order_by_loc(struct splsite **nodepp);
void insert_by_loc(struct splsite **headpp, struct splsite *nodep);
void prt_ss(FILE *fp, struct splsite *nodep, char seq[], int numbp, int ia, int rflag, int nmx);
void prt_clusters(FILE *fp, struct splsite *nodep, char seq[], int numbp, int nsites, int sscwl, int ia, int rflag);
void free_sslst(struct splsite *ssp);
void det_pdpa(struct splsite *nodep, float *pd, float *pa, char *seq, int numbp, int nmx, int unfrm);


#define LOG(f) ((f<=0)?-10000.0:(log(f)))

static void ComputLogValue(float *pd,float *pa,float *pdR,float *paR,int size)
{
    static unsigned int  allocated_log=0;
    float *f[4]={NULL};
    float *p,*v0,*v1;
    int i,j;
    
    f[0]=pd,f[1]=pa,f[2]=pdR,f[3]=paR;
    if(allocated_log<size){
        allocated_log=size;
        for(i=0;i<8;i++){
            ALLOCSPACE(logV[i],float,allocated_log);
        }
    }
    for(i=0;i<4;i++){
        v0=logV[2*i],v1=logV[2*i+1];
        p=f[i];
        for(j=0;j<size;j++){
            v0[j]=LOG(p[j]);
            v1[j]=LOG(1.-p[j]);
        }
        
    }
}    


int NALLOC_GCA = 0, NDEALLOC_GCA = 0;

int main(argc,argv)
int argc;
char *argv[];
{
int i,j,k, numbp, tmp, anum= 0, frmterr= 0;
int libop= 0, lnoffset = 4, gbflag = 0;
int fposset= 0, FROMPOS=0, tposset= 0, TOPOS=0;
static char gbfname[257], gbftstrg[257], dfname[257];
time_t tlc;
int pmt[23], tsmat[23][23];
char ch, ch2[2], aa[23], matname[257];

  sprm.pdg = 0.03F;
  sprm.ids = 2.0F;
  sprm.mms = -2.0F;
  sprm.nns = 0.0F;
  sprm.dls = -5.0F;
  sprm.unfrm = 0;
  sprm.min_intron_length = 30;
  sprm.min_exon_length = 5;
  sprm.tiny_exon = 20;
  sprm.short_exon = 200;
  sprm.long_intron = 300;
  sprm.poor_exon_score = 0.70;
  sprm.poor_donor_score = 0.5;
  sprm.poor_acptr_score = 0.5;
  sprm.join_length = 300;
  outfp = stdout;
  time(&tlc);

if (argc<2)
  {
	fprintf(stderr,USAGE1,argv[0]);
	fprintf(stderr,USAGE2);
	exit(0);
  }

for( i=1; i<argc; ++i )
  {if ( *argv[i]=='-' )   switch ( *(argv[i]+1) )   {
#ifdef DAPBM
      case 's':
	if (strcmp("human", argv[i + 1]) == 0)
	  Species = 0;
	else if (strcmp("mouse", argv[i + 1]) == 0)
	  Species = 1;
	else if (strcmp("rat", argv[i + 1]) == 0)
	  Species = 2;
	else if (strcmp("chicken", argv[i + 1]) == 0)
	  Species = 3;
	else if (strcmp("Drosophila", argv[i + 1]) == 0)
	  Species = 4;
	else if (strcmp("nematode", argv[i + 1]) == 0)
	  Species = 5;
	else if (strcmp("yeast", argv[i + 1]) == 0)
	  Species = 6;
	else if (strcmp("Aspergillus", argv[i + 1]) == 0)
	  Species = 7;
	else if (strcmp("Arabidopsis", argv[i + 1]) == 0)
	  Species = 8;
	else if (strcmp("maize", argv[i + 1]) == 0)
	  Species = 9;
	else if (strcmp("rice", argv[i + 1]) == 0)
	  Species = 10;
	else if (strcmp("Medicago", argv[i + 1]) == 0)
	  Species = 11;
	else if (strcmp("Populus", argv[i + 1]) == 0)
	  Species = 12;
	else if (strcmp("generic", argv[i + 1]) == 0)
	  Species = 99;
	else {
	  fprintf(stderr, "\nINVALID species %s. Exit.\n", argv[i + 1]);
	  exit(0);
	}
	anum += 2;
	break;
	case 'c':
                cutoff= atof(argv[i+1]);
		anum += 2; break;
#else
	case 'm':
                Model= atoi(argv[i+1]);
		if (Model < 0 || Model >= NMDLS)
		 {fprintf(stderr,"\nINVALID model %d. Exit.\n",Model);
		  exit(0);}
		anum+= 2; break;
	case 's':
		if (strcmp("maize", argv[i + 1]) == 0)
	  	Species = 0;
		else if (strcmp("Arabidopsis", argv[i + 1]) == 0)
	  	Species = 1;
		else if (strcmp("generic", argv[i + 1]) == 0)
	  	Species = 2;
		else {
	  	fprintf(stderr, "\nINVALID species %s. Exit.\n", argv[i + 1]);
	  	exit(0);
		}
		anum += 2;
		break;
	case 'c':
                Level= atoi(argv[i+1]);
		if (Level < 0 || Level >= NLEVELS)
		 {fprintf(stderr,"\nINVALID level %d. Exit.\n",Level);
		  exit(0);}
		anum += 2; break;
#endif
	case 't':
                pval= (double)atof(argv[i+1]); pvalset= 1; anum+= 2; break;
	case 'T':
                sval= (double)atof(argv[i+1]); svalset= 1; anum+= 2; break;
	case 'o':
                obsset= 1; ++anum; break;
	case 'O':
                obsset= 2; ++anum; break;
	case 'n':
                TOPN= atoi(argv[i+1]); anum+= 2; break;
#ifdef DAPBM
	case 'u':
                pruning= 0; ++anum; break;
	case 'U':
                allownn= 1; ++anum; break;
#endif
	case 'p':
                pstyle= atoi(argv[i+1]); anum+= 2; break;
	case 'w':
		prtclstrs= 1;
                if ( *argv[i+1]=='-' ) ++anum;
		else
		 {if (*argv[i+1]<49 || *argv[i+1]>57)
		   {
			   fprintf(stderr,USAGE1,argv[0]);
			   fprintf(stderr,USAGE2);
			   exit(0);
		   }
		  nsites= atoi(argv[i+1]); sscwl= atoi(argv[i+2]); anum+=3;
		 }
		break;
	case 'x':
                lxflag= 1;
                if ( *argv[i+1]=='-' ) ++anum;
		else
		 {if (*argv[i+1]<49 || *argv[i+1]>57)
		   {
			   fprintf(stderr,USAGE1,argv[0]);
			   fprintf(stderr,USAGE2);
			   exit(0);
		   }
		  dispset= 1;
		  DISPFR= atoi(argv[i+1])-1; DISPTO= atoi(argv[i+2])-1; anum+=3;
		 }
		break;
	case 'e':
                estop= 1; strcpy(estdbn,argv[i+1]); anum+= 2; break;
	case 'E':
                estop= 2; strcpy(estdbn,argv[i+1]); anum+= 2; break;
	case 'q':
                qpop= 1; qnoffset = 4; strcpy(qpfname,argv[i+1]); anum+= 2; break;
	case 'Q':
                qpop= 1; qnoffset = 1; strcpy(qpfname,argv[i+1]); anum+= 2; break;
	case 'i':
                prmop= 1; strcpy(prmfile,argv[i+1]); anum+= 2; break;
	case 'I':
                mflag= 1; strcpy(matname,argv[i+1]); anum+= 2; break;
	case 'a':
                FROMPOS= atoi(argv[i+1])-1; fposset= 1; anum+= 2; break;
	case 'b':
                TOPOS= atoi(argv[i+1])-1; tposset= 1; anum+= 2; break;
	case 'r':
                rflag= 1; bflag= 0; ++anum; break;
	case 'R':
                rflag= 0; bflag= 1; ++anum; break;
	case 'l':
                libop= 1;
		lnoffset = 4;
                strcpy(libname,argv[i+1]);
                anum+= 2; break;
	case 'L':
                libop= 1;
		lnoffset = 1;
                strcpy(libname,argv[i+1]);
                anum+= 2; break;
	case 'g':
                gbflag = 1; ++anum; break;
        default:
                break;
        }
  }

fprintf(stderr,"\nNOW EXECUTING:  ");
for( i=0; i<argc; ++i )   fprintf(stderr," %s", argv[i]);
fprintf(stderr,"\n");

  if (Species < 0) {
    fprintf(stderr,
            "\nYou must specify a species with the -s option to select the most appropriate\n\
               splice site models for your input.\n\n");
    fprintf(stderr,USAGE1,argv[0]);
    fprintf(stderr,USAGE2);
    exit(0);
  }

  if (!(libop + gbflag)) {
    fprintf(stderr,
            "\nPlease specify an input file by either the -l or the -g option.\n\n");
    fprintf(stderr,USAGE1,argv[0]);
    fprintf(stderr,USAGE2);
    exit(0);
  }

  if (!fsflag) strcpy(ftstrg,"all");
  if (strcmp(ftstrg,"all")==0)
   {strcpy(gbftstrg,"all"); lflag= 0; strcpy(dstrg,"sequence");
   }

  gcrta.mincolgth= 600;
  gcrta.minielgth= 3; gcrta.minmelgth= 30; gcrta.mintelgth= 3;
  gcrta.minilgth= 60; gcrta.maxilgth= 600;
  gcrta.minapop= -0.1;
  gcrta.mindscr= -0.1; gcrta.minascr= -0.1;

  if (pstyle!=3 && pstyle!=5)
#ifdef DAPBM
#ifdef DAPBM7
   {fprintf(outfp,"SplicePredictor.   Version of February 11, 2019.\n");
#else
   {fprintf(outfp,"SplicePredictorB2.   Version of February 11, 2019.\n");
#endif
#else
   {fprintf(outfp,"SplicePredictorLL.   Version of February 11, 2019.\n");
#endif
    fprintf(outfp,"Date run: %s\n", ctime(&tlc) );
#ifdef DAPBM
    if (Species < NMDLS)
      fprintf(outfp,"Species:\t\t\t%s\n",name_model[Species]);
    else
      fprintf(outfp,"Species:\t\t\tgeneric\n");
#ifdef DAPBM7
    if (Species <= 1)
      fprintf(outfp,"Model:\t\t\t\t2-class Bayesian\n");
    else
      fprintf(outfp,"Model:\t\t\t\t7-class Bayesian\n");
#else
    fprintf(outfp,"Model:\t\t\t\t2-class Bayesian\n");
#endif
#else
    fprintf(outfp,"Species:\t\t\t%s\n",species_name[Species]);
    fprintf(outfp,"Model:\t\t\t\t%s\n",model_name[Model]);
#endif
    if (pvalset) fprintf(outfp,"P-value prediction threshold:\t%6.2f\n",pval);
#ifdef DAPBM
    else         fprintf(outfp,"Prediction cutoff (2 ln[BF]):\t%6.2f\n",cutoff);
#else
    else         fprintf(outfp,"Prediction cutoff level:\t%1d\n",Level);
#endif
    if (svalset) fprintf(outfp,"S-value prediction threshold:\t%2d\n",sval);
#ifdef DAPBM
    if (pruning) fprintf(outfp,"Local pruning:\t\t\ton\n");
    else         fprintf(outfp,"Local pruning:\t\t\toff\n");
    if (allownn) fprintf(outfp,"Non-canonical sites:\t\tscored\n");
    else         fprintf(outfp,"Non-canonical sites:\t\tnot scored\n");
#endif
   }

  if (bflag) rflag= 0;

  if (prmop) {
    if ((prmfp = fopen(prmfile, "r")) == NULL) {
      fprintf(stderr, "File %s cannot be opened.\n", prmfile);
      perror(prmfile);
      exit(-1);
    }
    fprintf(stderr, "Parameter file %s has been opened for reading.\n", prmfile);
    read_sahmt_prm(prmfp, &sprm);
    fprintf(outfp, "Alignment parameters:");
    fprintf(outfp, "\n\tPDG	%5.2f", sprm.pdg);
    fprintf(outfp, "\n\tIDS	%5.2f", sprm.ids);
    fprintf(outfp, "\n\tMMS	%5.2f", sprm.mms);
    fprintf(outfp, "\n\tNNS	%5.2f", sprm.nns);
    fprintf(outfp, "\n\tDLS	%5.2f", sprm.dls);
    fprintf(outfp, "\n\tUNFRM %d", sprm.unfrm);
    fprintf(outfp, "\n\tMIN_INTRON_LENGTH	%3d", sprm.min_intron_length);
    fprintf(outfp, "\n\tMIN_EXON_LENGTH	%3d", sprm.min_exon_length);
    fprintf(outfp, "\n\tTINY_EXON	%3d", sprm.tiny_exon);
    fprintf(outfp, "\n\tSHORT_EXON	%3d", sprm.short_exon);
    fprintf(outfp, "\n\tLONG_INTRON	%3d", sprm.long_intron);
    fprintf(outfp, "\n\tPOOR_EXON_SCORE	%4.2f", sprm.poor_exon_score);
    fprintf(outfp, "\n\tPOOR_DONOR_SCORE	%4.2f", sprm.poor_donor_score);
    fprintf(outfp, "\n\tPOOR_ACPTR_SCORE	%4.2f", sprm.poor_acptr_score);
    fprintf(outfp, "\n\tJOIN_LENGTH	%d", sprm.join_length);
    fprintf(outfp, "\n\n");
    fclose(prmfp);
  }

  if (qpop) {
    if (openRawTextFile(QP_FILE, qpfname, 1) != 0) {
      fprintf(stderr, "File %s cannot be opened.\n", qpfname);
      exit(-1);
    }
    fprintf(stderr,
	    "Query protein library file %s has been opened for reading.\n", qpfname);

/******************************************************************************/
/* READ IN SCORING MATRIX:                                                    */
/*                                                                            */
    if (!mflag) {
      strcpy(matname,SMAT_NAME);
      for (i = 0; i < 23; ++i)
        for (j = 0; j < 23; ++j)
	  smat[i][j] = SMAT_PBLO62[i*23+j];
    }else{
      if ((matfp = fopen(matname, "r")) == NULL) {
        fprintf(stderr, "File %s cannot be opened.\n", matname);
        perror(matname);
        exit(-1);
      }
      fprintf(stderr, "Read scoring matrix file %s.\n", matname);

      i = 0;
      while ((ch = (char) (fgetc(matfp))) != '\n')
        if (ch > 62 && ch < 91)
	  aa[i++] = ch;
      if (i != 23) {
        fprintf(outfp, "\nERROR: Matrix %s has less than 23 column labels.\n", matname);
        exit(0);
      }
      for (i = maxsc = minsc = 0; i < 23; ++i) {
        fscanf(matfp, "%s", ch2);
        for (j = 0; j <= i; ++j) {
	  fscanf(matfp, "%d", &k);
	  if (k > maxsc)
	    maxsc = k;
	  if (k < minsc)
	    minsc = k;
	  tsmat[i][j] = tsmat[j][i] = k;
        }
      }
      for (i = 0; i < 23; ++i)
        for (j = 0; j < 23; ++j) {
	  if (aa[j] != AAUC[i])
	    continue;
	  pmt[i] = j;
        }
      for (i = 0; i < 23; ++i)
        for (j = 0; j < 23; ++j)
	  smat[i][j] = tsmat[pmt[i]][pmt[j]];

      if (maxsc - minsc + 1 > MAXRANGE) {
        fprintf(stderr,
	        "\nERROR: Matrix %s has a score range (%3d) exceeding the maximum (%3d).\n",
	        matname, maxsc - minsc + 1, MAXRANGE);
        exit(0);
      }
    }/* end of else */

    if (pstyle == 2 || mflag) {
      fprintf(outfp, "Amino acid substitution scoring matrix: %s \n", matname);
    }
/*                                                                            */
/******************************************************************************/
  }

  if (libop) {
    if (openRawTextFile(INPUT_FILE, libname, 1) != 0) {
      fprintf(stderr, "File %s cannot be opened.\n", libname);
      exit(-1);
    }
    fprintf(stderr, "Library file %s has been opened for reading.\n", libname);
    if (pstyle % 2 == 0)
      fprintf(outfp, "Input file:\t\t\t%s\n", libname);
    frmterr = 1;
    while ((numbp = getlns(INPUT_FILE, 1, lnoffset, sfname, gdna)) != 0) {
      if (numbp == -1)
	continue;
      if ((gdna = (char *) calloc(numbp , sizeof(char))) == NULL)
	 fatal_error("Error: memory allocation failed. Exit.\n");
      numbp = getlns(INPUT_FILE, 0, lnoffset, sfname, gdna);
      if (numbp < 1000000) ifwdth = 6; else ifwdth = 10;
      frmterr = 0;
      if (++fcount % 100 == 0)
	fprintf(stderr,
		" ... now processing nucleotide sequence %4d (%s)\n", fcount, sfname);
      if ((gdnaR = (char *) calloc(numbp , sizeof(char))) == NULL)
	 fatal_error("Error: memory allocation failed. Exit.\n");
      complement_seq(gdna, numbp, gdnaR);
      if (CHECKFLAG) {
	fprintf(outfp, "\n");
	for (i = 0; i < 80; ++i)
	  fprintf(outfp, "*");
	fprintf(outfp, "\nSequence %4d (File: %s)\n", fcount, sfname);
	if (rflag == 1) {
	  fprintf(outfp, "Reverse strand.\n");
	  prt_seq_segment(outfp, gdnaR, numbp, 0, numbp - 1, rflag, NAUC, 1);
	}
	else
	  prt_seq_segment(outfp, gdna, numbp, 0, numbp - 1, rflag, NAUC, 1);
      }
      if (!fposset || FROMPOS < 0)
	FROMPOS = 0;
      if (!tposset || TOPOS >= numbp)
	TOPOS = numbp - 1;
      if (FROMPOS > TOPOS) {
	tmp = TOPOS;
	TOPOS = FROMPOS;
	FROMPOS = tmp;
      }
      if (!fposset || FROMPOS < 0)
	FROMPOS = 0;
      if (!tposset || TOPOS >= numbp)
	TOPOS = numbp - 1;
      if (!dispset) {DISPFR= FROMPOS; DISPTO= TOPOS;}
      strcpy(dfname, sfname);
      for (fromp = FROMPOS; fromp < TOPOS; fromp += (SGMNTSZE - SGMNTSHFT)) {
        top = (fromp + SGMNTSZE - 1 > TOPOS) ? TOPOS : fromp + SGMNTSZE - 1;
        strcpy(sfname, dfname);
        doit(numbp, fromp, top, rflag);
        if (bflag) {
	  rflag = 1;
	  strcpy(sfname, dfname);
	  if (CHECKFLAG && fromp == FROMPOS) {
	    fprintf(outfp, "\n");
	    for (i = 0; i < 80; ++i)
	      fprintf(outfp, "*");
	    fprintf(outfp, "\nSequence %4d (File: %s)\n", fcount, sfname);
	    fprintf(outfp, "Reverse strand.\n");
	    prt_seq_segment(outfp, gdnaR, numbp, 0, numbp - 1, rflag, NAUC, 1);
	  }
          doit(numbp, fromp, top, rflag);
	  rflag = 0;
        }
        if (top == TOPOS)
          break;
      }
      free(gdna);
      if (rflag || bflag)
	free(gdnaR);
      fflush(outfp);
    }
    closeRawTextFile(INPUT_FILE);
    if (frmterr)
      fprintf(outfp,
	      "\nNo FASTA-formatted input sequence data found in file %s.\n",
	      libname);
  }

  if (gbflag)
    for (;;) {
      if (anum + 1 == argc)
	break;
      ++anum;
      if (openRawTextFile(INPUT_FILE, argv[anum], 1) != 0) {
	fprintf(stderr, "GenBank file %s cannot be opened.\n", argv[anum]);
	continue;
      }
      fprintf(stderr, "GenBank file %s has been opened for reading.\n", argv[anum]);
      ++gbfcount;
      frmterr = 1;
      strcpy(gbfname, argv[anum]);
      if (pstyle!=3 && pstyle!=5) {
        if (anum + 1 == argc && gbfcount == 1) {
	  fprintf(outfp, "GenBank file:\t\t\t%s\n", gbfname);
        }
        else {
	  fprintf(outfp, "\n\n");
	  for (i = 0; i < 80; ++i)
	    fprintf(outfp, "*");
	  fprintf(outfp, "\nGenBank file:\t\t\t%s (# %4d)\n", gbfname, gbfcount);
        }
      }
      while ((numbp = getgbs(INPUT_FILE, 1, sfname, gdna)) != 0) {
	if (numbp == -1)
	  continue;
	if ((gdna = (char *) calloc(numbp , sizeof(char))) == NULL) {
	  fatal_error("Error: memory allocation failed. Exit.\n");
	}
	numbp = getgbs(INPUT_FILE, 0, sfname, gdna);
        if (numbp < 1000000) ifwdth = 6; else ifwdth = 10;
	frmterr = 0;
	if (++fcount % 100 == 0) {
	  fprintf(stderr, " ... now processing nucleotide sequence %4d (%s)\n",
		  fcount, sfname);
	}
	time(&tlc);
	fprintf(stderr, "\nSeq #%3d = %s     %s", fcount, sfname, ctime(&tlc));
	if ((gdnaR = (char *) calloc(numbp , sizeof(char))) == NULL)
	   fatal_error("Error: memory allocation failed. Exit.\n");
	complement_seq(gdna, numbp, gdnaR);
	if (CHECKFLAG) {
	  fprintf(outfp, "\n");
	  for (i = 0; i < 80; ++i)
	    fprintf(outfp, "*");
	  fprintf(outfp, "\nSequence %4d (File: %s)\n", fcount, sfname);
	  if (rflag == 1)
	    fprintf(outfp, "Reverse strand.\n");
	  prt_seq_segment(outfp, gdna, numbp, 0, numbp - 1, rflag, NAUC, 1);
	}
	if (!fposset || FROMPOS < 0)
	  FROMPOS = 0;
	if (!tposset || TOPOS >= numbp)
	  TOPOS = numbp - 1;
	if (FROMPOS > TOPOS) {
	  tmp = TOPOS;
	  TOPOS = FROMPOS;
	  FROMPOS = tmp;
	}
	if (!fposset || FROMPOS < 0)
	  FROMPOS = 0;
	if (!tposset || TOPOS >= numbp)
	  TOPOS = numbp - 1;
        if (!dispset) {DISPFR= FROMPOS; DISPTO= TOPOS;}
	strcpy(dfname, sfname);
	for (fromp = FROMPOS; fromp < TOPOS; fromp += (SGMNTSZE - SGMNTSHFT)) {
          top = (fromp + SGMNTSZE - 1 > TOPOS) ? TOPOS : fromp + SGMNTSZE - 1;
          strcpy(sfname, dfname);
          doit(numbp, fromp, top, rflag);
	  if (bflag) {
	    rflag = 1;
	    strcpy(sfname, dfname);
	    if (CHECKFLAG) {
	      fprintf(outfp, "\n");
	      for (i = 0; i < 80; ++i)
	        fprintf(outfp, "*");
	      fprintf(outfp, "\nSequence %4d (File: %s)\n", fcount, sfname);
	      fprintf(outfp, "Reverse strand.\n");
	      prt_seq_segment(outfp, gdnaR, numbp, 0, numbp - 1, rflag, NAUC, 1);
	    }
	    doit(numbp, fromp, top, rflag);
            rflag = 0;
          }
          if (top == TOPOS)
            break;
        }
      free(gdna);
      if (rflag || bflag)
        free(gdnaR);
      fflush(outfp);
    }
    if (frmterr) {
      fprintf(outfp, "\nNo GenBank-formatted input sequence data found in file %s.\n",
	      gbfname);

    }
  closeRawTextFile(INPUT_FILE);
  }

  if (qpop)
    closeRawTextFile(QP_FILE);

for(i=0;i<8;i++)
    free(logV[i]);
exit(0);

} /* end main() */



void doit(numbp,ia,ib,rflag)
int numbp, ia, ib, rflag;
{
  int i, estlg, ecount, numaa, qpcnt;
  float boxX;
  struct splsite *ss_headp= NULL;
  char *cdna=NULL, *protein=NULL;
  char dfname[257], efname[257];
  auto int iaR;
  auto int ibR;
  auto int glgth;
  static float *pd;
  static float *pa;
  static float *pdR;
  static float *paR;
  static struct gcalgnmnt *gca;
  static struct gpalgnmnt *gpa;

  iaR = numbp - 1 - ib;
  ibR = numbp - 1 - ia;


if (pstyle!=3 && pstyle!=5)
 {fprintf(outfp,"\n");
  for (i=0;i<80;++i)   fprintf(outfp,"_");
  fprintf(outfp,"\nSequence %4d:   %s, ",fcount,sfname);
  if (rflag == 0)
    fprintf(outfp,"from %d to %d",ia+1,ib+1);
  else
    fprintf(outfp,"from %d to %d, reverse strand",ia+1,ib+1);
  fprintf(outfp,".\n\n");
 }

if (pstyle%2==0)
 {if (rflag == 0) {
    prt_seq_segment(outfp,gdna,numbp,ia,ib,rflag,NAUC,1);
    base_composition(gdna,ia,ib,&tbcv);
  }
  else {
    prt_seq_segment(outfp,gdnaR,numbp,iaR,ibR,rflag,NAUC,1);
    base_composition(gdnaR,iaR,ibR,&tbcv);
  }
  prt_base_composition(outfp,&tbcv);
 }

if (rflag == 0)
#ifdef DAPBM
  find_splsignalsBM(gdna,numbp,rflag,ia,ib,&ss_headp);
#else
  find_splsignalsLL(gdna,numbp,rflag,ia,ib,&ss_headp);
#endif
else
#ifdef DAPBM
  find_splsignalsBM(gdnaR,numbp,rflag,iaR,ibR,&ss_headp);
#else
  find_splsignalsLL(gdnaR,numbp,rflag,iaR,ibR,&ss_headp);
#endif

makedoublylinked(ss_headp);

#ifdef DAPBM
if (pruning) local_pruning(ss_headp,&ss_headp);
#endif

if (rflag == 0)
  local_optimality(gdna,numbp,rflag,ia,ib,ss_headp,&gcrta);
else
  local_optimality(gdnaR,numbp,rflag,iaR,ibR,ss_headp,&gcrta);
if (obsset==1) order_by_scr(&ss_headp);
else if (obsset==2) order_by_sstar(&ss_headp);

if (lxflag)
 {strcpy(lxfname,sfname); strcat(lxfname,".tex");
  if ( (latexfp = fopen( lxfname, "w" )) == NULL )
   {fprintf(stderr,"File %s cannot be opened.\n", lxfname);
    perror(lxfname); exit(-1);}
  fprintf(latexfp,"\\setlength{\\unitlength}{1.0pt}");
  fprintf(latexfp,"\n\\begin{picture}(%d,%d)(0,%d)\n",
		WIDTH,HEIGHT,-HEIGHT/2);
  fprintf(latexfp,"\\thicklines\n");
  fprintf(latexfp,"\\put(0,0){\\line(1,0){%d}}\n",WIDTH);
  fprintf(latexfp,"\\thinlines\n");
  if (rflag)
    fprintf(latexfp,"\\multiput(0,%.1f)(%d,0){%d}{\\line(1,0){%d}}\n",
		-4*(float)14/15 * (float)HEIGHT/8, WIDTH/80, 80, WIDTH/160);
  else
    fprintf(latexfp,"\\multiput(0,%.1f)(%d,0){%d}{\\line(1,0){%d}}\n",
		4*(float)14/15 * (float)HEIGHT/8, WIDTH/80, 80, WIDTH/160);
  fprintf(latexfp,"\\thicklines\n");
  if (rflag)
    fprintf(latexfp,"\\multiput(0,%.1f)(%d,0){%d}{\\line(1,0){%d}}\n",
		-4*(float)11/15 * (float)HEIGHT/8, WIDTH/40, 40, WIDTH/80);
  else
    fprintf(latexfp,"\\multiput(0,%.1f)(%d,0){%d}{\\line(1,0){%d}}\n",
		4*(float)11/15 * (float)HEIGHT/8, WIDTH/40, 40, WIDTH/80);
  fprintf(latexfp,"\\thinlines\n");
  if (rflag)
    fprintf(latexfp,"\\multiput(0,%.1f)(%d,0){%d}{\\line(1,0){%d}}\n",
		-4*(float)8/15 * (float)HEIGHT/8, WIDTH/80, 80, WIDTH/160);
  else
    fprintf(latexfp,"\\multiput(0,%.1f)(%d,0){%d}{\\line(1,0){%d}}\n",
		4*(float)8/15 * (float)HEIGHT/8, WIDTH/80, 80, WIDTH/160);
  fprintf(latexfp,"\\thicklines\n");
  if (rflag)
    fprintf(latexfp,"\\multiput(0,%.1f)(%d,0){%d}{\\line(1,0){%d}}\n",
		-4*(float)5/15 * (float)HEIGHT/8, WIDTH/40, 40, WIDTH/80);
  else
    fprintf(latexfp,"\\multiput(0,%.1f)(%d,0){%d}{\\line(1,0){%d}}\n",
		4*(float)5/15 * (float)HEIGHT/8, WIDTH/40, 40, WIDTH/80);
  fprintf(latexfp,"\\thinlines\n\n");
  SCALEF= (float)numbp/(float)(DISPTO-DISPFR+1);
  if (rflag==1)
   {boxX= 0.0;
    fprintf(latexfp,"\\put(%.1f,1){\\line(0,1){%.1f}}\n",
			boxX,(float)HEIGHT/8);
    fprintf(latexfp,"\\put(%.1f,%.1f){\\makebox(0,0)[b]{%d}}\n",
			boxX,(float)HEIGHT/8+4,DISPFR+1);
    for (i=DISPFR+100*((int)(2*(DISPTO-DISPFR+1)/1000)-
				(int)(DISPTO-DISPFR+1)/1000);i<=DISPTO+1;
	  i+=100*((int)(2*(DISPTO-DISPFR+1)/1000)-(int)(DISPTO-DISPFR+1)/1000))
     {boxX= (float)(i-DISPFR)/(float)numbp * WIDTH * SCALEF;
      fprintf(latexfp,"\\put(%.1f,1){\\line(0,1){%.1f}}\n",
			boxX,(float)HEIGHT/8);
      fprintf(latexfp,"\\put(%.1f,%.1f){\\makebox(0,0)[b]{%d}}\n",
			boxX,(float)HEIGHT/8+4,i);
     }
   }
  else
   {boxX= 0.0;
    fprintf(latexfp,"\\put(%.1f,%.1f){\\line(0,1){%.1f}}\n",
			boxX,-(float)HEIGHT/8-1,(float)HEIGHT/8);
    fprintf(latexfp,"\\put(%.1f,%.1f){\\makebox(0,0)[b]{%d}}\n",
			boxX,-(float)HEIGHT/8-10,DISPFR+1);
    for (i=DISPFR+100*((int)(2*(DISPTO-DISPFR+1)/1000)-
				(int)(DISPTO-DISPFR+1)/1000);i<=DISPTO+1;
	  i+=100*((int)(2*(DISPTO-DISPFR+1)/1000)-(int)(DISPTO-DISPFR+1)/1000))
     {boxX= (float)(i-DISPFR)/(float)numbp * WIDTH * SCALEF;
      fprintf(latexfp,"\\put(%.1f,%.1f){\\line(0,1){%.1f}}\n",
			boxX,-(float)HEIGHT/8-1,(float)HEIGHT/8);
      fprintf(latexfp,"\\put(%.1f,%.1f){\\makebox(0,0)[b]{%d}}\n",
			boxX,-(float)HEIGHT/8-10,i);
     }
   }
 }

if (obsset==0) TOPN= 0;
if (rflag == 0)
  prt_ss(outfp,ss_headp,gdna,numbp,ia,rflag,TOPN);
else
  prt_ss(outfp,ss_headp,gdnaR,numbp,numbp-ib,rflag,TOPN);
fflush(outfp);

if (prtclstrs)
 {if (rflag == 0)
    prt_clusters(outfp,ss_headp,gdna,numbp,nsites,sscwl,numbp-ib,rflag);
  else
    prt_clusters(outfp,ss_headp,gdnaR,numbp,nsites,sscwl,ia,rflag);
 }
fflush(outfp);

  if (estop || qpop) {
    strcpy(dfname, sfname);
    if (rflag == 0)
      strcat(dfname, "+");
    else
      strcat(dfname, "-");

    if (bflag == 0 || (bflag == 1 && rflag == 0)) {
      if ((pd = (float *) calloc((ib - ia + 2) , sizeof(float))) == NULL) {
        fatal_error("Error: memory allocation failed. Exit.\n");
      }
      if ((pa = (float *) calloc((ib - ia + 2) , sizeof(float))) == NULL) {
        fatal_error("Error: memory allocation failed. Exit.\n");
      }
#ifdef DAPBM
#ifdef DAPBM7
      det_daPbm7(&(gdna[ia]), ib - ia + 1, 0, ib - ia + 1, pd, pa, sprm.unfrm);
#else
      det_daPbm2(&(gdna[ia]), ib - ia + 1, 0, ib - ia + 1, pd, pa, sprm.unfrm);
#endif
#else
      det_daPll(&(gdna[ia]), ib - ia + 1, 0, ib - ia + 1, pd, pa, sprm.unfrm);
#endif
      if ((pdR = (float *) calloc((ib - ia + 2) , sizeof(float))) == NULL) {
        fatal_error("Error: memory allocation failed. Exit.\n");
      }
      if ((paR = (float *) calloc((ib - ia + 2) , sizeof(float))) == NULL) {
        fatal_error("Error: memory allocation failed. Exit.\n");
      }
#ifdef DAPBM
#ifdef DAPBM7
      det_daPbm7(&(gdnaR[iaR]), ib - ia + 1, 0, ib - ia + 1, pdR, paR, sprm.unfrm);
#else
      det_daPbm2(&(gdnaR[iaR]), ib - ia + 1, 0, ib - ia + 1, pdR, paR, sprm.unfrm);
#endif
#else
      det_daPll(&(gdnaR[iaR]), ib - ia + 1, 0, ib - ia + 1, pdR, paR, sprm.unfrm);
#endif
      ComputLogValue(pd,pa,pdR,paR,ib-ia+2);
    }
  }

  if (estop) {
    if (openRawTextFile(EST_FILE, estdbn, 1) != 0) {
  	fprintf(stderr, "File %s cannot be opened.\n", estdbn);
  	exit(-1);
    }
    fprintf(stderr, "EST library file %s has been opened for reading.\n",
  	      estdbn);
    if (pstyle <= 2  ||  pstyle == 4)
  	fprintf(outfp, "\n\n");
    fprintf(outfp, "EST library file:\t%s\n\n", estdbn);
    fflush(outfp);

    if ( (gca= (struct gcalgnmnt *) calloc(1 , sizeof(struct gcalgnmnt))) == NULL)
      fatal_error("Error: memory allocation failed. Exit.\n");
    if (CHECKFLAG) fprintf(stdout,"\nALLOC_GCA %ld # %d",(long)gca,++NALLOC_GCA);

    ecount = 0;
    while ((estlg = getlns(EST_FILE, 1, 1, sfname, cdna)) != 0) {
      if (estlg == -1)
        continue;
      if ((cdna = (char *) calloc(estlg , sizeof(char))) == NULL) {
        fatal_error("Error: memory allocation failed. Exit.\n");
      }
      estlg = getlns(EST_FILE, 0, 1, sfname, cdna);
      strcpy(efname, sfname);
      strcat(efname, "+");
      ++ecount;

      glgth = ib - ia + 1;
      if (glgth > MAXGLGTH) {
        fprintf(outfp,"\nWARNING:\n");
      }

      strcpy(gca->gname, dfname);
      strcpy(gca->cname, efname);
      strcpy(gca->estdbn, estdbn);
      gca->offset = 0;
      if (rflag == 0) {
        gca->gia = ia;
        gca->gib = ib;
      }
      else {
        gca->gia = iaR;
        gca->gib = ibR;
      }
      gca->clgth = estlg;
      if (rflag == 0) {
        if (sahmtD(outfp, gdna, gdnaR, ia, ib, numbp, rflag, pd, pa, pdR, paR, cdna, gca, sprm, ecount))
          fprintf(outfp, "%s", gca->algnmnt);
      }
      else {
        if (sahmtD(outfp, gdna, gdnaR, iaR, ibR, numbp, rflag, pd, pa, pdR, paR, cdna, gca, sprm, ecount))
          fprintf(outfp, "%s", gca->algnmnt);
      }
      if (estop == 2) {	/* ALIGN COMPLEMENTARY STRAND */
        efname[strlen(efname) - 1] = '-';
        reverse_seq(cdna, estlg);
        strcpy(gca->gname, dfname);
        strcpy(gca->cname, efname);
        strcpy(gca->estdbn, estdbn);
        gca->offset = 0;
        if (rflag == 0) {
          gca->gia = ia;
          gca->gib = ib;
        }
        else {
          gca->gia = iaR;
          gca->gib = ibR;
        }
        gca->clgth = estlg;
        if (rflag == 0) {
          if (sahmtD(outfp, gdna, gdnaR, ia, ib, numbp, rflag, pd, pa, pdR, paR, cdna, gca, sprm, ecount))
            fprintf(outfp, "%s", gca->algnmnt);
        }
        else {
          if (sahmtD(outfp, gdna, gdnaR, iaR, ibR, numbp, rflag, pd, pa, pdR, paR, cdna, gca, sprm, ecount))
            fprintf(outfp, "%s", gca->algnmnt);
        }
      }
      free(cdna);
    }
    closeRawTextFile(EST_FILE);
    free_gca(gca);
  }

  if (qpop) {
    setPositionRawTextFile(QP_FILE, 0);
   if ((gpa = (struct gpalgnmnt *) calloc(1 , sizeof(struct gpalgnmnt))) == NULL)
     fatal_error("Error: memory allocation failed. Exit.\n");
    qpcnt = 0;
    while ((numaa = getlps(QP_FILE, 1, qnoffset, sfname, protein)) != 0) {
      if (numaa == -1)
	continue;
      if ((protein = (char *) calloc((numaa + 1) , sizeof(char))) == NULL)
	fatal_error("Error: memory allocation failed. Exit.\n");
      numaa = getlps(QP_FILE, 0, qnoffset, sfname, protein);
      ++qpcnt;
      if (CHECKFLAG) {
	fprintf(outfp, "\n");
	for (i = 0; i < 80; ++i)
	  fprintf(outfp, "*");
	if (htmlop)
	  fprintf(outfp, "\nQuery protein sequence %4d (File: %s)\n", qpcnt, getProteinLink(sfname));
	else
	  fprintf(outfp, "\nQuery protein sequence %4d (File: %s)\n", qpcnt, sfname);
	prt_seq_segment(outfp, protein, numaa, 0, numaa - 1, 0, AAUC, 1);
	fflush(outfp);
      }
      protein[numaa] = 23;
      ++numaa;

      if (rflag == 0) {
        strcpy(gpa->gname, dfname);
        strcpy(gpa->pname, sfname);
        strcpy(gpa->qpdbn, qpfname);
        gpa->gia = ia;
        gpa->gib = ib;
        gpa->plgth = numaa;
        logPia=ia;
        if (sahmtP(outfp, &(gdna[ia]), dfname, ib - ia + 1, ia, numbp, rflag,
		   pd, pa, protein, sfname, numaa, smat, gpa, sprm, qpcnt)) {
          fprintf(outfp, "%s", gpa->algnmnt);
	}
      }
      else {
        if ((gpa = (struct gpalgnmnt *) calloc(1 , sizeof(struct gpalgnmnt))) == NULL)
          fatal_error("Error: memory allocation failed. Exit.\n");
        strcpy(gpa->gname, dfname);
        strcpy(gpa->pname, sfname);
        strcpy(gpa->qpdbn, qpfname);
        gpa->gia = iaR;
        gpa->gib = ibR;
        gpa->plgth = numaa;
        logPia=iaR;
        if (sahmtP(outfp, &(gdnaR[iaR]), dfname, ibR - iaR + 1, iaR, numbp, rflag,
		   pdR, paR, protein, sfname, numaa, smat, gpa, sprm, qpcnt)) {
          fprintf(outfp, "%s", gpa->algnmnt);
	}
      }		
      free(protein);
    }
    free_gpa(gpa);
  }				/* end if (qpop .. */

  if ((bflag == 0 && rflag == 0) || (bflag == 1 && rflag == 1)) {
    free(pd);
    free(pa);
  }
  if (rflag == 1) {
    free(pdR);
    free(paR);
  }

free_sslst(ss_headp);

} /* end doit() */


#ifdef DAPBM
void find_splsignalsBM(seq,numbp,rflag,ia,ib,ss_headpp)
char seq[];
int numbp, rflag, ia, ib;
struct splsite **ss_headpp;
{
int i;
double is_donor(), dscr, is_acptr(), ascr;
double dmin, amin;
struct splsite *nodep, *tmp=NULL;

for (i=ia;i<=ib;++i)
 {if (i-Wlgth >= 0  &&  i+Wlgth+1 < numbp)
   {if (pvalset) dmin= pval; else dmin= 0.0;
    if ((dscr=is_donor(seq,numbp,i,Wlgth,Species,&dsite,allownn)) >= dmin  &&
	 dsite.cval >= cutoff                                                )
     {nodep= (struct splsite *) calloc(1 , sizeof(struct splsite));
      nodep->stype= 1; nodep->loc= i; nodep->scr= dsite.pval;
      nodep->cval= dsite.cval; nodep->delU= dsite.delU; nodep->delS= dsite.delS;
      nodep->istr= 0; nodep->vote= 0;
      nodep->nexts = (struct splsite*) NULL;
      if (*ss_headpp == NULL)
       {*ss_headpp = nodep; tmp = nodep;}
      else
       {tmp->nexts = nodep; tmp = nodep;}
     }
   }
  if (i-Wlgth-1 >= 0  &&  i+Wlgth < numbp)
   {if (pvalset) amin= pval; else amin= 0.0;
    if ((ascr=is_acptr(seq,numbp,i,Wlgth,Species,&asite,allownn)) >= amin  &&
	 asite.cval >= cutoff                                                )
     {nodep= (struct splsite *) calloc(1 , sizeof(struct splsite));
      nodep->stype= 2; nodep->loc= i; nodep->scr= asite.pval;
      nodep->cval= asite.cval; nodep->delU= asite.delU; nodep->delS= asite.delS;
      nodep->istr= 0; nodep->vote= 0;
      nodep->nexts = (struct splsite*) NULL;
      if (*ss_headpp == NULL)
       {*ss_headpp = nodep; tmp = nodep;}
      else
       {tmp->nexts = nodep; tmp = nodep;}
     }
   }
 }

} /* end find_splsignalsBM() */


#else
void find_splsignalsLL(seq,numbp,rflag,ia,ib,ss_headpp)
char seq[];
int numbp, rflag, ia, ib;
struct splsite **ss_headpp;
{
int i, smodel=0;
double is_donor(), dscr, is_acptr(), ascr;
double dmin, amin;
struct splsite *nodep, *tmp=NULL;

for (i=ia;i<=ib;++i)
 {if (i-Wlgth >= 0  &&  i+Wlgth+1 < numbp)
   {if (Model==0) smodel= 0;
    else if (Model==1)
     {if (seq[i-1] == 3) smodel= 1;   else smodel= 2;}
    if (pvalset) dmin= pval; else dmin= dtv[Species][smodel][Level];
    if ((dscr=is_donor(seq,numbp,i,Wlgth,Species,smodel,&dsite)) >= dmin)
     {nodep= (struct splsite *) calloc(1 , sizeof(struct splsite));
      nodep->stype= 1 ; nodep->loc= i; nodep->scr= dscr;
      nodep->Lscr= dsite.Lscr,nodep->delU= dsite.delU,nodep->delS= dsite.delS;
      nodep->istr= 0; nodep->vote= 0;
      nodep->nexts = (struct splsite*) NULL;
      if (*ss_headpp == NULL)
       {*ss_headpp = nodep; tmp = nodep;}
      else
       {tmp->nexts = nodep; tmp = nodep;}
     }
   }
  if (i-Wlgth-1 >= 0  &&  i+Wlgth < numbp)
   {if (Model==0) smodel= 0;
    else if (Model==1)
     {if (seq[i-2] == 1) smodel= 1;   else smodel= 2;}
    if (pvalset) amin= pval; else amin= atv[Species][smodel][Level];
    if ((ascr=is_acptr(seq,numbp,i,Wlgth,Species,smodel,&asite)) >= amin)
     {nodep= (struct splsite *) calloc(1 , sizeof(struct splsite));
      nodep->stype= 2 ; nodep->loc= i; nodep->scr= ascr;
      nodep->Lscr= asite.Lscr,nodep->delU= asite.delU,nodep->delS= asite.delS;
      nodep->istr= 0; nodep->vote= 0;
      nodep->nexts = (struct splsite*) NULL;
      if (*ss_headpp == NULL)
       {*ss_headpp = nodep; tmp = nodep;}
      else
       {tmp->nexts = nodep; tmp = nodep;}
     }
   }
 }

} /* end find_splsignalsLL() */
#endif



void makedoublylinked(nodep)
struct splsite *nodep;
{
struct splsite *lastn= NULL;

while (nodep != NULL)
 {nodep->prevs= lastn;
  lastn= nodep;
  nodep= nodep->nexts;
 }

} /* end makedoublylinked() */



void local_pruning(nodep,headpp)
struct splsite *nodep, **headpp;
{
struct splsite *tmpp, *ntmpp;
float high_cutoff =  6.0, low_cutoff = 4.0;
int ldelta = 100, sdelta = 50;
int uflag;

while (nodep != NULL) {
  if (nodep->cval < high_cutoff) {	/* ... look at high-probability sites */
    nodep= nodep->nexts;
    continue;
  }
  uflag= 1;
  if (CHECKFLAG)
    fprintf(stderr,"\ntype %d at %d %f",nodep->stype,nodep->loc+1,nodep->cval);
  if (nodep->stype == 1) {	/* ... donor site  with cval >= high_cutoff */
    tmpp= nodep->nexts;
    while (tmpp != NULL) {
      if (tmpp->stype == 1) {
	tmpp= tmpp->nexts;
	continue;
      }
      if (tmpp->stype == 2) {
	if (tmpp->loc > nodep->loc + gcrta.minilgth) break;
	if (tmpp->cval >= low_cutoff) {	/* ... good acceptor site nearby
					   downstream */
	  uflag= 0;
	  break;
	}
	tmpp= tmpp->nexts;
      }
    }
    if (uflag == 0) {			/* ... do nothing with high-probability
					   donor sites with good acceptor sites
					   nearby downstream */
      nodep= nodep->nexts;
      continue;
    }
    tmpp= nodep->prevs;
    while (tmpp != NULL) {
      if (tmpp->stype == 2) {
	tmpp= tmpp->prevs;
	continue;
      }
      if (CHECKFLAG)
	fprintf(stderr,"\nprev %d at %d %f",tmpp->stype,tmpp->loc+1,tmpp->cval);
      if (tmpp->loc < nodep->loc - ldelta) break;
      if (tmpp->cval < 0.5 * nodep->cval) {	/* ... remove nearby upstream
						   weak donor sites */
	ntmpp= tmpp->prevs;
        if (CHECKFLAG)
	  fprintf(stderr,"\nrm %d at %d %f",tmpp->stype,tmpp->loc+1,tmpp->cval);
	rm_splsite(tmpp,headpp);
	tmpp= ntmpp;
	continue;
      }
      tmpp= tmpp->prevs;
    }
    tmpp= nodep->nexts;
    if (CHECKFLAG  &&  tmpp != NULL)
      fprintf(stderr,"\nnext %d at %d %f",tmpp->stype,tmpp->loc+1,tmpp->cval);
    while (tmpp != NULL) {
      if (tmpp->stype == 2) {		/* ... remove nearby downstream weak
					   acceptor sites */
	if (tmpp->loc < nodep->loc + sdelta  &&  tmpp->cval < low_cutoff) {
	  ntmpp= tmpp->nexts;
          if (CHECKFLAG)
	    fprintf(stderr,"\nrm %d at %d %f",tmpp->stype,tmpp->loc+1,tmpp->cval);
	  rm_splsite(tmpp,headpp);
	  tmpp= ntmpp;
	  continue;
        }
	tmpp= tmpp->nexts;
	continue;
      }
      if (CHECKFLAG)
	fprintf(stderr,"\nnext %d at %d %f",tmpp->stype,tmpp->loc+1,tmpp->cval);
      if (tmpp->loc > nodep->loc + ldelta) break;
      if (tmpp->cval < 0.5 * nodep->cval) {	/* ... remove nearby downstream
						   weak donor sites */
	ntmpp= tmpp->nexts;
        if (CHECKFLAG)
	  fprintf(stderr,"\nrm %d at %d %f",tmpp->stype,tmpp->loc+1,tmpp->cval);
	rm_splsite(tmpp,headpp);
	tmpp= ntmpp;
	continue;
      }
      tmpp= tmpp->nexts;
    }
  }
  else {			/* acptr with cval >= high_cutoff */
    uflag= 1;
    tmpp= nodep->prevs;
    while (tmpp != NULL) {
      if (tmpp->stype == 2) {
	tmpp= tmpp->prevs;
	continue;
      }
      if (tmpp->stype == 1) {
	if (tmpp->loc < nodep->loc - gcrta.minilgth) break;
	if (tmpp->cval >= low_cutoff) {	/* ... good donor site nearby
					   upstream */
	  uflag= 0;
	  break;
	}
	tmpp= tmpp->prevs;
      }
    }
    if (uflag == 0) {			/* ... do nothing with high-probability
					   acceptor sites with good donor sites
					   nearby upstream */
      nodep= nodep->nexts;
      continue;
    }
    tmpp= nodep->prevs;
    while (tmpp != NULL) {
      if (tmpp->stype == 1) {
	if (tmpp->loc > nodep->loc - sdelta  &&  tmpp->cval < low_cutoff) {
	  ntmpp= tmpp->prevs;		/* ... remove nearby upstream weak
					   donor sites */
	  rm_splsite(tmpp,headpp);
	  tmpp= ntmpp;
	  continue;
        }
	tmpp= tmpp->prevs;
	continue;
      }
      if (tmpp->loc < nodep->loc - ldelta) break;
      if (tmpp->cval < 0.5 * nodep->cval) {	/* ... remove nearby upstream
						   weak acceptor sites */
	ntmpp= tmpp->prevs;
	rm_splsite(tmpp,headpp);
	tmpp= ntmpp;
	continue;
      }
      tmpp= tmpp->prevs;
    }
    tmpp= nodep->nexts;
    while (tmpp != NULL) {
      if (tmpp->stype == 1) {
	tmpp= tmpp->nexts;
	continue;
      }
      if (tmpp->loc > nodep->loc + ldelta) break;
      if (tmpp->cval < 0.5 * nodep->cval) {	/* ... remove nearby downstream
						   weak acceptor sites */
	ntmpp= tmpp->nexts;
	rm_splsite(tmpp,headpp);
	tmpp= ntmpp;
	continue;
      }
      tmpp= tmpp->nexts;
    }
  }
  nodep= nodep->nexts;
}

} /* end local_pruning() */



void rm_splsite(nodep,headpp)
struct splsite *nodep, **headpp;
{

if (nodep->prevs != NULL)
  (nodep->prevs)->nexts = nodep->nexts;
else
  *headpp = nodep->nexts;
if (nodep->nexts != NULL)
  (nodep->nexts)->prevs = nodep->prevs;
free((struct splsite *) nodep);

} /* end rm_splsite() */



void local_optimality(seq,numbp,rflag,ia,ib,nodep,crta)
char seq[]; int numbp, rflag, ia, ib;
struct splsite *nodep;
struct eicrta *crta;
{
int i, aloc=0, dloc=0, mloc=0, sloc=0, delta=0, sstar=0;
struct splsite *ip, *jp;
double ascr, dscr, mscr, oscr=0.0, sscr;
int nls, nrs, lgth, ADstrg[MAXADL], ADloc[MAXADL], RWstrg[MAXADL];
double ADscr[MAXADL];
static char ADsymb[5]= "ADEI";

while (nodep != NULL)
 {pscr[0]= pscr[1]= pscr[2]= pscr[3]= 0.0;
  delta= crta->maxilgth;
  if (nodep->stype == 1)	/* donor site */
   {sscr= 0.0;   ascr= -9.9;   jp= nodep->nexts;
    while (jp != NULL  &&  jp->loc < nodep->loc + delta)
     {if (jp->stype == 2)	/* acceptor */
       {if (jp->loc > nodep->loc + crta->minilgth  &&  jp->scr > ascr)
         {aloc= jp->loc; ascr= jp->scr;}
       }
      else			/* donor */
       {if (jp->loc > nodep->loc + crta->minilgth  &&  jp->scr >= nodep->scr)
	  delta= jp->loc - nodep->loc + crta->minilgth;
       }
      jp= jp->nexts;
     }
    if (ascr > -9.9)
     {oscr= nodep->scr * ascr;   sscr+= oscr;
      if (pstyle==4) fprintf(outfp,
	"\n\nDONOR at %d (%5.3f) MAXACCEPTOR at %d (%5.3f) score %f",
	nodep->loc+1,nodep->scr,aloc+1,ascr,oscr);
      ip= nodep->prevs;		/* to the left */
      delta= crta->maxilgth;
      while (ip != NULL)
       {if (ip->stype == 1)	/* donor */
         {if (ip->loc > aloc - delta)
           {mscr= -9.9;   jp= ip->nexts;
	    while (jp != NULL)
	     {if (jp->stype == 2  &&  jp->loc > ip->loc + crta->minilgth)
               {if (jp->loc <= aloc)
                 {if (jp->scr > mscr) {mscr= jp->scr; mloc= jp->loc;}
                 }
                else
	          break;
	       }
              jp= jp->nexts;
	     }
	    if (mscr > -9.9  &&  mloc > nodep->loc)
	     {if (pstyle==4) fprintf(outfp,
		"\nLONOR at %d (%5.3f) MAXACCEPTOR at %d (%5.3f) score %f",
		ip->loc+1,ip->scr,mloc+1,mscr,ip->scr*mscr);
	      sscr+= (ip->scr * mscr);
	     }
           }
          else
	    break;
         }
	else
	 {if (ip->scr > ascr) delta= aloc - ip->loc + crta->minilgth;
	 }
        ip= ip->prevs;
       }
      jp= nodep->nexts;		/* to the right */
      while (jp != NULL)
       {if (jp->stype == 1)
         {if (jp->loc < aloc - crta->minilgth)
	   {if (pstyle==4) fprintf(outfp,
		"\nRONOR at %d (%5.3f) MAXACCEPTOR at %d (%5.3f) score %f",
		jp->loc+1,jp->scr,aloc+1,ascr,jp->scr*ascr);
	    sscr+= (jp->scr * ascr);
           }
          else if (jp->loc < aloc)
	   {mscr= -9.9;   ip= nodep->nexts;
	    delta= crta->maxilgth;
	    while (ip != NULL &&  ip->loc < jp->loc + delta)
	     {if (ip->stype == 2)	/* acceptor */
	       {if (ip->loc > jp->loc + crta->minilgth  &&  ip->scr > mscr)
	         {sloc= ip->loc; mscr= ip->scr;}
	       }
	      else
	       {if (ip->scr > jp->scr)
		  delta= ip->loc - jp->loc + crta->minilgth;
	       }
	      ip= ip->nexts;
	     }
            if (mscr > -9.9)
	     {if (pstyle==4) fprintf(outfp,
		"\nRRNOR at %d (%5.3f) MAXACCEPTOR at %d (%5.3f) score %f",
		jp->loc+1,jp->scr,sloc+1,mscr,jp->scr*mscr);
	      sscr+= (jp->scr * mscr);
             }
	   }
          else
	    break;
         }
        jp= jp->nexts;
       }
     }
    else
      if (pstyle==4) fprintf(outfp,
		"\n\nDONOR at %d (%5.3f); no suitable acceptor",
		nodep->loc+1,nodep->scr);
    if (sscr > 0)
      nodep->rho= (oscr/sscr)*oscr;
    else
      nodep->rho= 0.0;
   }
  else		/* acceptor site */
   {sscr= 0.0;   dscr= -9.9;   ip= nodep->prevs;	/* <- best donor */
    while (ip != NULL  &&  ip->loc > nodep->loc - delta)
     {if (ip->stype == 1)	/* donor */
       {if (ip->loc < nodep->loc - crta->minilgth  &&  ip->scr > dscr)
	 {dloc= ip->loc; dscr= ip->scr;}
       }
      else			/* acceptor */
       {if (ip->loc < nodep->loc - crta->minilgth  &&  ip->scr >= nodep->scr)
	  delta= nodep->loc - ip->loc + crta->minilgth;
       }
      ip= ip->prevs;
     }
    if (dscr > -9.9)
     {oscr= nodep->scr * dscr;   sscr+= oscr;
      if (pstyle==4) fprintf(outfp,
		"\n\nACPTR at %d (%5.3f) MAXDONOR at %d (%5.3f) score %f",
		nodep->loc+1,nodep->scr,dloc+1,dscr,oscr);
      jp= nodep->nexts;		/* to the right */
      delta= crta->maxilgth;
      while (jp != NULL)
       {if (jp->stype == 2)	/* acceptor */
         {if (jp->loc < dloc + delta)
           {mscr= -9.9;   ip= jp->prevs;
	    while (ip != NULL)
	     {if (ip->stype == 1  &&  ip->loc < jp->loc - crta->minilgth)
	       {if (ip->loc >= dloc)
                 {if (ip->scr > mscr) {mscr= ip->scr; mloc= ip->loc;}
	         }
                else
	          break;
	       }
              ip= ip->prevs;
	     }
	    if (mscr > -9.9  &&  mloc < nodep->loc)
	     {if (pstyle==4) fprintf(outfp,
		"\nRCPTR at %d (%5.3f) MAXDONOR at %d (%5.3f) score %f",
		jp->loc+1,jp->scr,mloc+1,mscr,jp->scr*mscr);
	      sscr+= (jp->scr * mscr);
 	     }
           }
          else
	    break;
         }
        else
         {if (jp->scr > dscr) delta= jp->loc - dloc + crta->minilgth;
         }
        jp= jp->nexts;
       }
      ip= nodep->prevs;		/* to the left */
      while (ip != NULL)
       {if (ip->stype == 2)
         {if (ip->loc > dloc + crta->minilgth)
	   {if (pstyle==4) fprintf(outfp,
		"\nLCPTR at %d (%5.3f) MAXDONOR at %d (%5.3f) score %f",
		ip->loc+1,ip->scr,dloc+1,dscr,ip->scr*dscr);
	    sscr+= (ip->scr * dscr);
 	   }
          else if (ip->loc > dloc)
	   {mscr= -9.9;   jp= nodep->prevs;
	    delta= crta->maxilgth;
	    while (jp != NULL &&  jp->loc > ip->loc - delta)
	     {if (jp->stype == 1)	/* donor */
	       {if (jp->loc < ip->loc - crta->minilgth  &&  jp->scr > mscr)
	         {sloc= jp->loc; mscr= jp->scr;}
	       }
	      else
	       {if (jp->scr > ip->scr)
		  delta= ip->loc - jp->loc - crta->minilgth;
	       }
	      jp= jp->prevs;
	     }
            if (mscr > -9.9)
	     {if (pstyle==4) fprintf(outfp,
		"\nLLPTR at %d (%5.3f) MAXDONOR at %d (%5.3f) score %f",
		ip->loc+1,ip->scr,sloc+1,mscr,ip->scr*mscr);
	      sscr+= (ip->scr * mscr);
             }
	   }
          else
	    break;
         }
        ip= ip->prevs;
       }
     }
    else
      if (pstyle==4) fprintf(outfp,
		"\n\nACPTR at %d (%5.3f); no suitable donor",
		nodep->loc+1,nodep->scr);
    if (sscr > 0)
      nodep->rho= (oscr/sscr)*oscr;
    else
      nodep->rho= 0.0;
   }
  if (pstyle==4) fprintf(outfp,"\n");

  ip= nodep->prevs;   nls= 0;
  while (ip != NULL)
   {if (++nls == ADLGTH/2) break;
    if (ip->prevs == NULL) break;
    ip= ip->prevs;
   }
  jp= nodep->nexts;   nrs= 0;
  while (jp != NULL)
   {if (++nrs == ADLGTH/2) break;
    if (jp->nexts == NULL) break;
    jp= jp->nexts;
   }
  if (nls>nrs) nls= nrs;
  ip= nodep->prevs;   lgth= 0;
  while (ip != NULL)
   {if (++lgth == nls) break;
    ip= ip->prevs;
   }
  jp= ip;   lgth= nrs= 0;
  if (jp == NULL)
   {if (nodep->stype == 1) ADstrg[0]= 1; else ADstrg[0]= 0;}
  else while (jp != NULL)
   {if (jp->stype == 1) ADstrg[lgth]= 1;   else ADstrg[lgth]= 0;
    ADloc[lgth]= jp->loc;   ADscr[lgth]= jp->scr;
    if (++lgth > nls) ++nrs;
    if (nrs ==  nls+1) break;
    jp= jp->nexts;
   }
  nodep->glgth= lgth;
  if (re_write(ADstrg,ADloc,ADscr,RWstrg,0,lgth) == 0)
    nodep->ADstrg[0]= ADstrg[0];
  if (pstyle==4)
    fprintf(outfp,"\tA %5.3f D %5.3f E %5.3f I %5.3f\n",
		pscr[0],pscr[1],pscr[2],pscr[3]);
  if (nodep->stype == 1)
   {if (pscr[1]>pscr[2])
     {if (pscr[2]>pscr[3])
       {nodep->gamma= pscr[1]-pscr[2];
        for (i=0;i<lgth;++i) nodep->ADstrg[i]= MAXRWstrg[1][i];
       }
      else
       {if (pscr[1]>pscr[3])
	 {nodep->gamma= pscr[1]-pscr[3];
          for (i=0;i<lgth;++i) nodep->ADstrg[i]= MAXRWstrg[1][i];
	 }
	else
	 {nodep->gamma= 0.0;
          for (i=0;i<lgth;++i) nodep->ADstrg[i]= MAXRWstrg[3][i];
	 }
       }
     }
    else
     {nodep->gamma= 0.0;
       {if (pscr[2]>pscr[3])
          for (i=0;i<lgth;++i) nodep->ADstrg[i]= MAXRWstrg[2][i];
	else
          for (i=0;i<lgth;++i) nodep->ADstrg[i]= MAXRWstrg[3][i];
       }
     }
   }
  else
   {if (pscr[0]>pscr[2])
     {if (pscr[2]>pscr[3])
       {nodep->gamma= pscr[0]-pscr[2];
        for (i=0;i<lgth;++i) nodep->ADstrg[i]= MAXRWstrg[0][i];
       }
      else
       {if (pscr[0]>pscr[3])
	 {nodep->gamma= pscr[0]-pscr[3];
          for (i=0;i<lgth;++i) nodep->ADstrg[i]= MAXRWstrg[0][i];
	 }
	else
	 {nodep->gamma= 0.0;
          for (i=0;i<lgth;++i) nodep->ADstrg[i]= MAXRWstrg[3][i];
	 }
       }
     }
    else
     {nodep->gamma= 0.0;
       {if (pscr[2]>pscr[3])
          for (i=0;i<lgth;++i) nodep->ADstrg[i]= MAXRWstrg[2][i];
	else
          for (i=0;i<lgth;++i) nodep->ADstrg[i]= MAXRWstrg[3][i];
       }
     }
   }
  if (nodep->stype == 1)
   {for (i=0;i<4;++i)
     {if (nodep->scr >= dqntl[Species][0][i]) break;
     }
    sstar= 5-i; nodep->pstar= 5-i;
    for (i=0;i<4;++i)
     {if (nodep->rho >= dqntl[Species][1][i]) break;
     }
    sstar+= (5-i); nodep->rstar= 5-i;
    for (i=0;i<4;++i)
     {if (nodep->gamma >= dqntl[Species][2][i]) break;
     }
    sstar+= (5-i); nodep->gstar= 5-i;
   }
  else
   {for (i=0;i<4;++i)
     {if (nodep->scr >= aqntl[Species][0][i]) break;
     }
    sstar= 5-i; nodep->pstar= 5-i;
    for (i=0;i<4;++i)
     {if (nodep->rho >= aqntl[Species][1][i]) break;
     }
    sstar+= (5-i); nodep->rstar= 5-i;
    for (i=0;i<4;++i)
     {if (nodep->gamma >= aqntl[Species][2][i]) break;
     }
    sstar+= (5-i); nodep->gstar= 5-i;
   }
  nodep->sstar= sstar;

  if (lgth==ADLGTH)
   {ip= nodep;
    for (i=0;i<lgth/2;++i) ip= ip->prevs;
    for (i=0;i<lgth;++i)
     {if ((ip->stype==1 && nodep->ADstrg[i]==1) ||
          (ip->stype==2 && nodep->ADstrg[i]==0)   ) ++(ip->vote);
      ip= ip->nexts;
     }
   }

  if (pstyle==4)
   {fprintf(outfp,"%2d",sstar);
    if (lgth<ADLGTH) fprintf(outfp," -\t"); else fprintf(outfp," +\t");
#ifdef DAPBM
    fprintf(outfp,"%6.2f\t%5.2f\t%5.2f\t%5.3f\t%6.3f\t%6.3f",
  	nodep->cval,nodep->delU,nodep->delS,nodep->scr,
  	nodep->rho,nodep->gamma);
#else
    fprintf(outfp,"%6.2f\t%5.2f\t%5.2f\t%5.3f\t%6.3f\t%6.3f",
  	nodep->Lscr,nodep->delU,nodep->delS,nodep->scr,
  	nodep->rho,nodep->gamma);
#endif
    if (lgth<ADLGTH) fprintf(outfp," -\t"); else fprintf(outfp," +\t");
    for (i=0;i<lgth/2;++i) fprintf(outfp,"%c",ADsymb[nodep->ADstrg[i]]);
    fprintf(outfp,"-%c-",ADsymb[nodep->ADstrg[lgth/2]]);
    for (i=lgth/2+1;i<lgth;++i) fprintf(outfp,"%c",ADsymb[nodep->ADstrg[i]]);
    fprintf(outfp,"\t%2d",nodep->vote);
    if (lgth<ADLGTH) fprintf(outfp," -\t"); else fprintf(outfp," +\t");
    fprintf(outfp,"%2d",nodep->sstar+nodep->vote);
    if (lgth<ADLGTH) fprintf(outfp," -"); else fprintf(outfp," +");
    fprintf(outfp,"\n");
   }

  nodep= nodep->nexts;
 }

} /* end local_optimality() */



int re_write(ADstrg,ADloc,ADscr,RWstrg,i,lgth)
int *ADstrg, *ADloc, *RWstrg, i, lgth;
double *ADscr;

{
int j;

if (lgth==0) return(0);

if (ADstrg[i]==0)
 {if (i==0)
   {RWstrg[i]= 0;
    re_write_next(ADstrg,ADloc,ADscr,RWstrg,i,lgth);
    RWstrg[i]= 2;
    re_write_next(ADstrg,ADloc,ADscr,RWstrg,i,lgth);
    RWstrg[i]= 3;
    re_write_next(ADstrg,ADloc,ADscr,RWstrg,i,lgth);
   }
  else
   {if (RWstrg[i-1]==0  ||  RWstrg[i-1]==2)
     {RWstrg[i]= 2;
      re_write_next(ADstrg,ADloc,ADscr,RWstrg,i,lgth);
     }
    else
     {for (j=i-1;j>=0;--j) {if (RWstrg[j]==1) break;}
      if (j<0  ||  ADloc[i]-ADloc[j]>60)
       {RWstrg[i]= 0;
        re_write_next(ADstrg,ADloc,ADscr,RWstrg,i,lgth);
       }
      RWstrg[i]= 3;
      re_write_next(ADstrg,ADloc,ADscr,RWstrg,i,lgth);
     }
   }
 }
else
 {if (i==0)
   {RWstrg[i]= 1;
    re_write_next(ADstrg,ADloc,ADscr,RWstrg,i,lgth);
    RWstrg[i]= 2;
    re_write_next(ADstrg,ADloc,ADscr,RWstrg,i,lgth);
    RWstrg[i]= 3;
    re_write_next(ADstrg,ADloc,ADscr,RWstrg,i,lgth);
   }
  else
   {if (RWstrg[i-1]==1  ||  RWstrg[i-1]==3)
     {RWstrg[i]= 3;
      re_write_next(ADstrg,ADloc,ADscr,RWstrg,i,lgth);
     }
    else
     {RWstrg[i]= 1;
      re_write_next(ADstrg,ADloc,ADscr,RWstrg,i,lgth);
      RWstrg[i]= 2;
      re_write_next(ADstrg,ADloc,ADscr,RWstrg,i,lgth);
     }
   }
 }
return(1);

} /* end re_write() */



void re_write_next(ADstrg,ADloc,ADscr,RWstrg,i,lgth)
int ADstrg[], ADloc[], RWstrg[], i, lgth;
double ADscr[];
{
void scr_RWstrg();

if (i<lgth-1)   re_write(ADstrg,ADloc,ADscr,RWstrg,i+1,lgth);
else            scr_RWstrg(RWstrg,lgth,ADstrg,ADloc,ADscr);

} /* end re_write_next() */



void scr_RWstrg(RWstrg,lgth,ADstrg,ADloc,ADscr)
int RWstrg[], lgth, ADstrg[], ADloc[];
double ADscr[];
{
int i;
static char RWsymb[5]= "ADEI";
double scr= 0.0;

for (i=0;i<lgth;++i)
 {if (pstyle==4) fprintf(outfp,"%c",RWsymb[RWstrg[i]]);
  if (RWstrg[i]==ADstrg[i]) scr+= ADscr[i];
 }
if (scr > pscr[RWstrg[lgth/2]])
 {pscr[RWstrg[lgth/2]]= scr;
  for (i=0;i<lgth;++i)
   MAXRWstrg[RWstrg[lgth/2]][i]= RWstrg[i];
 }

if (pstyle==4) fprintf(outfp,"\t%f\n",scr);

} /* end scr_RWstrg() */



void order_by_scr(nodepp)
struct splsite **nodepp;
{
struct splsite *ip= *nodepp, *tmp, *newheadp= NULL;

while (ip != NULL)
 {tmp= ip->nexts;
  ip->nexts = NULL;
  insert_by_scr(&(newheadp),ip);
  ip= tmp;
 }

*nodepp= newheadp;

} /* end order_by_scr() */



void insert_by_scr(headpp,nodep)
struct splsite **headpp, *nodep;
{

if (*headpp == NULL)
 {*headpp= nodep; nodep->nexts= NULL;}
else
 {while ( (*headpp)->nexts != NULL  && (*headpp)->scr >= nodep->scr)
    headpp= &((*headpp)->nexts);
  if ((*headpp)->nexts == NULL)
   {if ((*headpp)->scr >= nodep->scr)
     {(*headpp)->nexts= nodep; nodep->nexts= NULL;}
    else
     {nodep->nexts= *headpp; *headpp= nodep;}
   }
  else
   {nodep->nexts= *headpp; *headpp= nodep;}
 }

} /* end insert_by_scr */



void order_by_sstar(nodepp)
struct splsite **nodepp;
{
struct splsite *ip= *nodepp, *tmp, *newheadp= NULL;

while (ip != NULL)
 {tmp= ip->nexts;
  ip->nexts = NULL;
  insert_by_sstar(&(newheadp),ip);
  ip= tmp;
 }

*nodepp= newheadp;

} /* end order_by_sstar() */



void insert_by_sstar(headpp,nodep)
struct splsite **headpp, *nodep;
{

if (*headpp == NULL)
 {*headpp= nodep; nodep->nexts= NULL;}
else
 {while ( (*headpp)->nexts != NULL  && (*headpp)->sstar >= nodep->sstar)
    headpp= &((*headpp)->nexts);
  if ((*headpp)->nexts == NULL)
   {if ((*headpp)->sstar >= nodep->sstar)
     {(*headpp)->nexts= nodep; nodep->nexts= NULL;}
    else
     {nodep->nexts= *headpp; *headpp= nodep;}
   }
  else
   {nodep->nexts= *headpp; *headpp= nodep;}
 }

} /* end insert_by_sstar */



void order_by_loc(nodepp)
struct splsite **nodepp;
{
struct splsite *ip= *nodepp, *tmp, *newheadp= NULL;

while (ip != NULL)
 {tmp= ip->nexts;
  ip->nexts = NULL;
  insert_by_loc(&(newheadp),ip);
  ip= tmp;
 }

*nodepp= newheadp;

} /* end order_by_loc() */



void insert_by_loc(headpp,nodep)
struct splsite **headpp, *nodep;
{

if (*headpp == NULL)
 {*headpp= nodep; nodep->nexts= NULL;}
else
 {while ( (*headpp)->nexts != NULL  && (*headpp)->loc > nodep->loc )
    headpp= &((*headpp)->nexts);
  if ((*headpp)->nexts == NULL)
   {if ((*headpp)->loc > nodep->loc)
     {(*headpp)->nexts= nodep; nodep->nexts= NULL;}
    else
     {nodep->nexts= *headpp; *headpp= nodep;}
   }
  else
   {nodep->nexts= *headpp; *headpp= nodep;}
 }

} /* end insert_by_loc */



void prt_ss(fp,nodep,seq,numbp,ia,rflag,nmx)
FILE *fp;
struct splsite *nodep;
char seq[]; int numbp, ia, rflag, nmx;
{
int i, nd= 0, na= 0;
static char ADsymb[5]= "ADEI";

if (pstyle==1 || pstyle%2==0)
 {fprintf(fp,"\nPotential splice sites");
  if (obsset==1) fprintf(fp," (ordered by P-value");
  if (obsset==2) fprintf(fp," (ordered by *-value");
  if (obsset)
   {if (nmx) fprintf(fp," ; top %d scoring sites)",nmx); else fprintf(fp,")");}
  fprintf(fp,"\n\n");
  if (pstyle==1)
   {fprintf(fp,"t    q   ");
    for (i=0;i<ifwdth-3;++i) fprintf(fp," ");
    fprintf(fp,"loc     sequence           P      ");
#ifdef DAPBM
    fprintf(fp," c      ");
#endif
    fprintf(fp,"rho   gamma   *  P*R*G*        parse\n\n");
   }
  else
   {fprintf(fp,"t    q   ");
    for (i=0;i<ifwdth-3;++i) fprintf(fp," ");
    fprintf(fp,"loc m      sequence           P         ");
#ifdef DAPBM
    fprintf(fp,"c     U     S      rho   gamma   *  P*R*G*        parse\n\n");
#else
    fprintf(fp,"L     U     S      rho   gamma   *  P*R*G*        parse\n\n");
#endif
   }
 }

while (nodep != NULL)
 {if (svalset  &&  nodep->sstar < sval)
   {nodep= nodep->nexts; continue;}
  if (prtclstrs==0 && lxflag==1)
   {if (rflag==1)
     {if (numbp-nodep->loc >= DISPFR && numbp-nodep->loc <= DISPTO)
       {if (nodep->stype==1)
         {fprintf(latexfp,"\\put(%.1f,-1){\\line(0,-1){%.1f}}\n",
		(float)(numbp-nodep->loc-DISPFR)/(float)numbp * WIDTH * SCALEF,
		4 * (float)nodep->sstar/15 * (float)HEIGHT/8 );
          fprintf(latexfp,"\\put(%.1f,%.1f){\\vector(-1,0){0.0}}\n",
		(float)(numbp-nodep->loc-DISPFR)/(float)numbp * WIDTH * SCALEF
			- VLGTH,
		-4 * (float)nodep->sstar/15 * (float)HEIGHT/8 );
         }
        else
         {fprintf(latexfp,"\\put(%.1f,-1){\\line(0,-1){%.1f}}\n",
		(float)(numbp-nodep->loc-DISPFR)/(float)numbp * WIDTH * SCALEF,
		4 * (float)nodep->sstar/15 * (float)HEIGHT/8 );
          fprintf(latexfp,"\\put(%.1f,%.1f){\\vector(1,0){0.0}}\n",
		(float)(numbp-nodep->loc-DISPFR)/(float)numbp * WIDTH * SCALEF
			+ VLGTH,
		-4 * (float)nodep->sstar/15 * (float)HEIGHT/8 );
         }
       }
     }
    else
     {if (nodep->loc >= DISPFR && nodep->loc <= DISPTO)
       {if (nodep->stype==1)
         {fprintf(latexfp,"\\put(%.1f,1){\\line(0,1){%.1f}}\n",
		(float)(nodep->loc+1-DISPFR)/(float)numbp * WIDTH * SCALEF,
		4 * (float)nodep->sstar/15 * (float)HEIGHT/8 );
          fprintf(latexfp,"\\put(%.1f,%.1f){\\vector(1,0){0.0}}\n",
		(float)(nodep->loc+1-DISPFR)/(float)numbp * WIDTH * SCALEF
			+ VLGTH,
		4 * (float)nodep->sstar/15 * (float)HEIGHT/8 );
         }
        else
         {fprintf(latexfp,"\\put(%.1f,1){\\line(0,1){%.1f}}\n",
		(float)(nodep->loc+1-DISPFR)/(float)numbp * WIDTH * SCALEF,
		4 * (float)nodep->sstar/15 * (float)HEIGHT/8 );
          fprintf(latexfp,"\\put(%.1f,%.1f){\\vector(-1,0){0.0}}\n",
		(float)(nodep->loc+1-DISPFR)/(float)numbp * WIDTH * SCALEF
			- VLGTH,
		4 * (float)nodep->sstar/15 * (float)HEIGHT/8 );
         }
       }
     }
   }
  if (nodep->stype == 1)
   {if (nmx && ++nd>nmx)
     {nodep= nodep->nexts; continue;}
    else
     {if (pstyle==1)
       {if (rflag==1)
         {if (nodep->sstar>=14)      fprintf(fp,"D <----- ");
          else if (nodep->sstar>=11) fprintf(fp,"D  <---- ");
          else if (nodep->sstar>=8)  fprintf(fp,"D   <--- ");
          else if (nodep->sstar>=5)  fprintf(fp,"D    <-- ");
          else                       fprintf(fp,"D     <- ");
          fprintf(fp,"%*d",ifwdth,numbp-nodep->loc);
	 }
        else
         {if (nodep->sstar>=14)      fprintf(fp,"D -----> ");
          else if (nodep->sstar>=11) fprintf(fp,"D ---->  ");
          else if (nodep->sstar>=8)  fprintf(fp,"D --->   ");
          else if (nodep->sstar>=5)  fprintf(fp,"D -->    ");
          else                       fprintf(fp,"D ->     ");
          fprintf(fp,"%*d",ifwdth,nodep->loc +1);
	 }
        fprintf(fp," ");
       }
      if (pstyle%2==0)
       {if (nmx) fprintf(fp,"D#%3d ",nd);
	else     fprintf(fp,"D ");
        if (rflag==1)
         {if (nodep->sstar>=14)      fprintf(fp,"<----- ");
          else if (nodep->sstar>=11) fprintf(fp," <---- ");
          else if (nodep->sstar>=8)  fprintf(fp,"  <--- ");
          else if (nodep->sstar>=5)  fprintf(fp,"   <-- ");
          else                       fprintf(fp,"    <- ");
          fprintf(fp,"%*d %1d  ",ifwdth,numbp-nodep->loc,(numbp-nodep->loc+1)%3);
	 }
        else
         {if (nodep->sstar>=14)      fprintf(fp,"-----> ");
          else if (nodep->sstar>=11) fprintf(fp,"---->  ");
          else if (nodep->sstar>=8)  fprintf(fp,"--->   ");
          else if (nodep->sstar>=5)  fprintf(fp,"-->    ");
          else                       fprintf(fp,"->     ");
          fprintf(fp,"%*d %1d  ",ifwdth,nodep->loc+1,(nodep->loc+1)%3);
	 }
       }
      if (pstyle==3)
       {if (rflag==1)
          fprintf(fp,"D %*d\n",ifwdth,numbp-nodep->loc);
        else
          fprintf(fp,"D %*d\n",ifwdth,nodep->loc +1);
       }
     }
    if (pstyle!=3 && pstyle!=5)
     {fprintf(fp,"          ");
      for (i=nodep->loc-3;i<nodep->loc;++i)
	fprintf(fp,"%c",NALC[(int)seq[i]]);
      for (              ;i<nodep->loc+2;++i)
	fprintf(fp,"%c",NAUC[(int)seq[i]]);
      for (              ;i<nodep->loc+6;++i)
	fprintf(fp,"%c",NALC[(int)seq[i]]);
      fprintf(fp,"  %6.3f",nodep->scr);
#ifdef DAPBM
      if (pstyle==1) fprintf(fp," %6.2f ",nodep->cval);
#endif
      if (pstyle%2==0)
#ifdef DAPBM
        fprintf(fp,"  (%6.2f %5.2f %5.2f)",
		nodep->cval,nodep->delU,nodep->delS);
#else
        fprintf(fp,"  (%6.2f %5.2f %5.2f)",
		nodep->Lscr,nodep->delU,nodep->delS);
#endif
      fprintf(fp," %6.3f %6.3f",nodep->rho,nodep->gamma);
      fprintf(fp,"  %2d (%1d %1d %1d)",nodep->sstar,
			nodep->pstar,nodep->rstar,nodep->gstar);
      fprintf(fp,"  ");
      for (i=0;i<ADLGTH/2-nodep->glgth/2;++i) fprintf(fp," ");
      for (i=0;i<nodep->glgth/2;++i)
	fprintf(fp,"%c",ADsymb[nodep->ADstrg[i]]);
      fprintf(fp,"-%c-",ADsymb[nodep->ADstrg[nodep->glgth/2]]);
      for (i=nodep->glgth/2+1;i<nodep->glgth;++i)
	fprintf(fp,"%c",ADsymb[nodep->ADstrg[i]]);
      fprintf(fp,"\n");
     }
    if (pstyle==5)
     {fprintf(fp,"%1d\tDONOR\t%s\t",nodep->istr,sfname);
      if (rflag==1)
        fprintf(fp,"%*d\t",ifwdth,numbp-nodep->loc);
      else
        fprintf(fp,"%*d\t",ifwdth,nodep->loc+1);
     }
   }
  else
   {if (nmx && ++na>nmx)
     {nodep= nodep->nexts; continue;}
    else
     {if (pstyle==1)
       {if (rflag==1)
         {if (nodep->sstar>=14)      fprintf(fp,"A -----> ");
          else if (nodep->sstar>=11) fprintf(fp,"A ---->  ");
          else if (nodep->sstar>=8)  fprintf(fp,"A --->   ");
          else if (nodep->sstar>=5)  fprintf(fp,"A -->    ");
          else                       fprintf(fp,"A ->     ");
          fprintf(fp,"%*d",ifwdth,numbp-nodep->loc);
	 }
        else
         {if (nodep->sstar>=14)      fprintf(fp,"A <----- ");
          else if (nodep->sstar>=11) fprintf(fp,"A  <---- ");
          else if (nodep->sstar>=8)  fprintf(fp,"A   <--- ");
          else if (nodep->sstar>=5)  fprintf(fp,"A    <-- ");
          else                       fprintf(fp,"A     <- ");
          fprintf(fp,"%*d",ifwdth,nodep->loc +1);
         }
        fprintf(fp," ");
       }
      if (pstyle%2==0)
       {if (nmx) fprintf(fp,"A#%3d ",na);
	else     fprintf(fp,"A ");
        if (rflag==1)
         {if (nodep->sstar>=14)      fprintf(fp,"-----> ");
          else if (nodep->sstar>=11) fprintf(fp,"---->  ");
          else if (nodep->sstar>=8)  fprintf(fp,"--->   ");
          else if (nodep->sstar>=5)  fprintf(fp,"-->    ");
          else                       fprintf(fp,"->     ");
          fprintf(fp,"%*d %1d  ",ifwdth,numbp-nodep->loc,(numbp-nodep->loc+1)%3);
         }
	else
         {if (nodep->sstar>=14)      fprintf(fp,"<----- ");
          else if (nodep->sstar>=11) fprintf(fp," <---- ");
          else if (nodep->sstar>=8)  fprintf(fp,"  <--- ");
          else if (nodep->sstar>=5)  fprintf(fp,"   <-- ");
          else                       fprintf(fp,"    <- ");
          fprintf(fp,"%*d %1d  ",ifwdth,nodep->loc+1,(nodep->loc+1)%3);
         }
       }
      if (pstyle==3)
       {if (rflag==1)
          fprintf(fp,"A %*d\n",ifwdth,numbp-nodep->loc);
        else
          fprintf(fp,"A %*d\n",ifwdth,nodep->loc +1);
       }
     }
    if (pstyle!=3 && pstyle!=5)
     {for (i=nodep->loc-14;i<nodep->loc-1;++i)
	fprintf(fp,"%c",NALC[(int)seq[i]]);
      for (              ;i<nodep->loc+1;++i)
	fprintf(fp,"%c",NAUC[(int)seq[i]]);
      for (              ;i<nodep->loc+3;++i)
	fprintf(fp,"%c",NALC[(int)seq[i]]);
      fprintf(fp,"    %6.3f",nodep->scr);
#ifdef DAPBM
      if (pstyle==1) fprintf(fp," %6.2f ",nodep->cval);
#endif
      if (pstyle%2==0)
#ifdef DAPBM
        fprintf(fp,"  (%6.2f %5.2f %5.2f)",
		nodep->cval,nodep->delU,nodep->delS);
#else
        fprintf(fp,"  (%6.2f %5.2f %5.2f)",
		nodep->Lscr,nodep->delU,nodep->delS);
#endif
      fprintf(fp," %6.3f %6.3f",nodep->rho,nodep->gamma);
      fprintf(fp,"  %2d (%1d %1d %1d)",nodep->sstar,
			nodep->pstar,nodep->rstar,nodep->gstar);
      fprintf(fp,"  ");
      for (i=0;i<ADLGTH/2-nodep->glgth/2;++i) fprintf(fp," ");
      for (i=0;i<nodep->glgth/2;++i)
	fprintf(fp,"%c",ADsymb[nodep->ADstrg[i]]);
      fprintf(fp,"-%c-",ADsymb[nodep->ADstrg[nodep->glgth/2]]);
      for (i=nodep->glgth/2+1;i<nodep->glgth;++i)
	fprintf(fp,"%c",ADsymb[nodep->ADstrg[i]]);
      fprintf(fp,"\n");
     }
    if (pstyle==5)
     {fprintf(fp,"%1d\tACPTR\t%s\t",nodep->istr,sfname);
      if (rflag==1)
        fprintf(fp,"%*d\t",ifwdth,numbp-nodep->loc);
      else
        fprintf(fp,"%*d\t",ifwdth,nodep->loc+1);
     }
   }
  if (pstyle==5)
#ifdef DAPBM
   {fprintf(fp,"%6.2f\t%5.2f\t%5.2f\t%5.3f\t%6.3f\t%6.3f",
	nodep->cval,nodep->delU,nodep->delS,nodep->scr,
	nodep->rho,nodep->gamma);
#else
   {fprintf(fp,"%6.2f\t%5.2f\t%5.2f\t%5.3f\t%6.3f\t%6.3f",
	nodep->Lscr,nodep->delU,nodep->delS,nodep->scr,
	nodep->rho,nodep->gamma);
#endif
    if (nodep->glgth<ADLGTH) fprintf(fp," -\t"); else fprintf(fp," +\t");
    for (i=0;i<nodep->glgth/2;++i) fprintf(fp,"%c",ADsymb[nodep->ADstrg[i]]);
    fprintf(fp,"-%c-",ADsymb[nodep->ADstrg[nodep->glgth/2]]);
    for (i=nodep->glgth/2+1;i<nodep->glgth;++i)
      fprintf(fp,"%c",ADsymb[nodep->ADstrg[i]]);
    fprintf(fp,"\t%2d",nodep->pstar);
    fprintf(fp,"\t%2d",nodep->rstar);
    fprintf(fp,"\t%2d",nodep->gstar);
    if (nodep->glgth<ADLGTH) fprintf(fp," -\t"); else fprintf(fp," +\t");
    fprintf(fp,"%2d",nodep->sstar);
    if (nodep->glgth<ADLGTH) fprintf(fp," -\t"); else fprintf(fp," +\t");
    fprintf(fp,"%2d",nodep->vote);
    if (nodep->glgth<ADLGTH) fprintf(outfp," -\t"); else fprintf(outfp," +\t");
    fprintf(fp,"%2d",nodep->sstar+nodep->vote);
    if (nodep->glgth<ADLGTH) fprintf(outfp," -"); else fprintf(outfp," +");
    fprintf(outfp,"\n");
   }
  nodep= nodep->nexts;
 }
if (pstyle==3) fprintf(fp,"q\n");
if (prtclstrs==0 && lxflag==1)
 {fprintf(latexfp,"\\end{picture}\n");
  fclose(latexfp);
 }

} /* end prt_ss() */



void free_sslst(ssp)
struct splsite *ssp;
{
struct splsite *tmp;

while (ssp != NULL)
 {tmp= ssp->nexts;
  free((struct splsite *) ssp);
  ssp= tmp;
 }

} /* end free_sslst() */



void prt_clusters(fp,nodep,seq,numbp,nsites,sscwl,ia,rflag)
FILE *fp;
struct splsite *nodep;
char seq[]; int numbp, nsites, sscwl, ia, rflag;
{
int ntd, nta, nfd, nfa;
int wbeg= 0, pwend=0, wend= numbp;
struct splsite *ip, *wp=NULL;
float boxX, boxL, boxH;

fprintf(fp,"\n\n\nPotential splice site clusters:\n\n");

sscwl-= 1;
while (nodep != NULL)
 {if (nodep->loc > numbp-(sscwl/2)  ||  (svalset  &&  nodep->sstar < sval))
   {nodep= nodep->nexts; continue;}
  if (lxflag)
   {if (rflag==1)
     {if (nodep->stype==1)
        fprintf(latexfp,"\\put(%.1f,-1){\\vector(-1,-4){%.1f}}\n",
		(float)(numbp-nodep->loc-ia)/(float)numbp * WIDTH * SCALEF,
		(float)nodep->sstar/15 * (float)HEIGHT/8 );
      else
        fprintf(latexfp,"\\put(%.1f,-1){\\vector(1,-4){%.1f}}\n",
		(float)(numbp-nodep->loc-ia)/(float)numbp * WIDTH * SCALEF,
		(float)nodep->sstar/15 * (float)HEIGHT/8 );
     }
    else
     {if (nodep->stype==1)
        fprintf(latexfp,"\\put(%.1f,1){\\vector(1,4){%.1f}}\n",
		(float)(nodep->loc+1-ia)/(float)numbp * WIDTH * SCALEF,
		(float)nodep->sstar/15 * (float)HEIGHT/8 );
      else
        fprintf(latexfp,"\\put(%.1f,1){\\vector(-1,4){%.1f}}\n",
		(float)(nodep->loc+1-ia)/(float)numbp * WIDTH * SCALEF,
		(float)nodep->sstar/15 * (float)HEIGHT/8 );
     }
   }
  ntd= nta= nfd = nfa= 0;
  ip= nodep;
  while (ip != NULL  &&  ip->loc <= nodep->loc + sscwl)
   {if (svalset  &&  ip->sstar < sval)
     {ip= ip->nexts; continue;}
    if (ip->stype == 1)
     {if (ip->istr == 1) ++ntd; else ++nfd;}
    else
     {if (ip->istr == 1) ++nta; else ++nfa;}
    pwend= ip->loc;
    ip= ip->nexts;
   }
  if (ntd+nta+nfa+nfd >= nsites)
   {if (nodep->loc > wend-(sscwl/2)  &&  pwend>wend)
     {ip= wp;   ntd= nta= nfd = nfa= 0;
      while (ip != NULL  &&  ip->loc <= wend)
       {if (svalset  &&  ip->sstar < sval)
         {ip= ip->nexts; continue;}
        if (ip->stype == 1)
         {if (ip->istr == 1) ++ntd; else ++nfd;}
        else
         {if (ip->istr == 1) ++nta; else ++nfa;}
        ip= ip->nexts;
       }
      if (rflag==1)
        fprintf(fp,"W %*d %*d\t",ifwdth,numbp-wbeg,ifwdth,numbp-wend);
      else
        fprintf(fp,"W %*d %*d\t",ifwdth,wbeg+1,ifwdth,wend+1);
      fprintf(fp,"D= %2d A= %2d S= %2d  wl= %4d\n",nfd,nfa,nfd+nfa,wend-wbeg+1);
      if (lxflag)
       {if (rflag==1)
	 {boxX= (float)(numbp-wend-ia)/(float)numbp * WIDTH * SCALEF;
	  boxL= (float)(wend-wbeg)/(float)numbp * WIDTH * SCALEF;
	  boxH= (float)(ntd+nta+nfa+nfd)/(float)nsites *
		(float)(sscwl)/(float)(wend-wbeg+1) * (float)HEIGHT/2;
	  fprintf(latexfp,"\\put(%.1f,%.1f){\\framebox(%.1f,%.1f)}\n",
			boxX,-boxH-1,boxL,boxH);
	 }
        else
	 {boxX= (float)(wbeg+1-ia)/(float)numbp * WIDTH * SCALEF;
	  boxL= (float)(wend-wbeg)/(float)numbp * WIDTH * SCALEF;
	  boxH= (float)(ntd+nta+nfa+nfd)/(float)nsites *
		(float)(sscwl)/(float)(wend-wbeg+1) * (float)HEIGHT/2;
	  fprintf(latexfp,"\\put(%.1f,1){\\framebox(%.1f,%.1f)}\n",
			boxX,boxL,boxH);
	 }
       }
      wp= nodep;   wbeg= nodep->loc;   wend= pwend;
     }
    else
     {if (wbeg==0)
       {wp= nodep;   wbeg= nodep->loc;
       }
      wend= pwend;
     }
   }
  nodep= nodep->nexts;
 }
if (wbeg>0)
 {ip= wp;   ntd= nta= nfd = nfa= 0;
  while (ip != NULL  &&  ip->loc <= wend)
   {if (svalset  &&  ip->sstar < sval)
     {ip= ip->nexts; continue;}
    if (ip->stype == 1)
     {if (ip->istr == 1) ++ntd; else ++nfd;}
    else
     {if (ip->istr == 1) ++nta; else ++nfa;}
    ip= ip->nexts;
   }
  if (rflag==1)
    fprintf(fp,"W %*d %*d\t",ifwdth,numbp-wbeg,ifwdth,numbp-wend);
  else
    fprintf(fp,"W %*d %*d\t",ifwdth,wbeg+1,ifwdth,wend+1);
  fprintf(fp,"D= %2d A= %2d S= %2d  wl= %4d\n",nfd,nfa,nfd+nfa,wend-wbeg+1);
  if (lxflag)
   {if (rflag==1)
     {boxX= (float)(numbp-wend-ia)/(float)numbp * WIDTH * SCALEF;
      boxL= (float)(wend-wbeg)/(float)numbp * WIDTH * SCALEF;
      boxH= (float)(ntd+nta+nfa+nfd)/(float)nsites *
		(float)(sscwl)/(float)(wend-wbeg+1) * (float)HEIGHT/2;
      fprintf(latexfp,"\\put(%.1f,%.1f){\\framebox(%.1f,%.1f)}\n",
		boxX,-boxH-1,boxL,boxH);
     }
    else
     {boxX= (float)(wbeg+1-ia)/(float)numbp * WIDTH * SCALEF;
      boxL= (float)(wend-wbeg)/(float)numbp * WIDTH * SCALEF;
      boxH= (float)(ntd+nta+nfa+nfd)/(float)nsites *
		(float)(sscwl)/(float)(wend-wbeg+1) * (float)HEIGHT/2;
      fprintf(latexfp,"\\put(%4.1f,1){\\framebox(%4.1f,%4.1f)}\n",
		boxX,boxL,boxH);
     }
   }
 }
if (lxflag)
 {fprintf(latexfp,"\\end{picture}\n");
  fclose(latexfp);
 }

} /* end prt_clusters() */



void det_pdpa(nodep,pd,pa,seq,numbp,nmx,unfrm)
struct splsite *nodep;
float pd[], pa[];
char seq[]; int numbp, nmx, unfrm;
{
int i, nd= 0, na= 0;

pd[0]= pa[0]= 0.000001;
for (i=1;i<numbp-1;++i)
 {if (seq[i]==3)
   {if (seq[i+1]==0)      pd[i]= 0.001;
    else if (seq[i+1]==1) pd[i]= 0.0001;
    else                  pd[i]= 0.000001;
    if (seq[i-1]==2)      pa[i]= 0.001;
    else                  pa[i]= 0.000001;
   }
  else
   {pd[i]= pa[i]= 0.000001;}
 }
pd[numbp-1]= pa[numbp-1]= pd[numbp]= pa[numbp]= 0.000001;
if (unfrm) return;

while (nodep != NULL)
 {if (svalset  &&  nodep->sstar < sval)
   {nodep= nodep->nexts; continue;}
  if (nodep->stype == 1)
   {if (nmx && ++nd>nmx)
     {nodep= nodep->nexts; continue;}
    pd[nodep->loc]= (float)nodep->scr;
   }
  else
   {if (nmx && ++na>nmx)
     {nodep= nodep->nexts; continue;}
    pa[nodep->loc]= (float)nodep->scr;
   }
  nodep= nodep->nexts;
 }

} /* end det_pdpa() */



void fatal_error(char *buf)
{

	fprintf(stderr,"%s\n",buf);
	exit(1);

} /* end fatal_error() */
