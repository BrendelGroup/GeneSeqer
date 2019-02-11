/* GeneSeqer.c;                                Last update: February 11, 2019. */
/* Dependencies:   getgbs.c getlns.c getlps.c sahmtd.c sahmtp.c               */
/* Bugs:                                                                      */

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



#define DAPBM

#ifdef DAPBM
#define USAGE1 "Usage:\n\n\
    %s [-hv] -s species [-dD dbest(s)] [-eE estdbn] [-k] [-m maxnest]\n\
    [-M sgmntsze] [-x wsize] [-y minqHSP] [-z minqHSPc] [-w minESTc] [-p prmfile] [-qQ qpdbn(s)] [-I matname]\n\
    [-a from] [-b to] [-frR] [-oO outfname] [-lL libfname] [-g gbfname(s)]\n\n\
  -h           :   HTML output [default: simple text]\n\
  -v           :   Verbose output\n\
  -s species   :   Set species to select the most appropriate splice site models.\n\
                    This parameter must be specified. Options:\n\
                    'human', 'mouse', 'rat', 'chicken', 'Drosophila', 'Daphnia', 'nematode'\n\
                    'yeast', 'Aspergillus', 'Arabidopsis', 'maize', 'rice',\n\
                    'Medicago', 'Populus', 'generic'.\n\
  -dD dbest(s) :   Read EST sequence data from pre-processed library file(s)\n\
                    'dbest(s)' (FASTA-format; '-d': >gi|name; '-D': >name).\n\
  -eE estdbn   :   Read EST sequence data from library file 'estdbn'.\n\
                    (FASTA-format; '-e': >gi|name; '-E': >name).\n\
  -k           :   Assume fixed 5' to 3' transcript orientation of supplied EST sequences.\n\
  -m maxnest   :   Display at most 'maxnest' spliced EST alignments per genomic\n\
                    DNA input [default: 500].\n\
  -M sgmntsze  :   Process input sequence in segments of size 'sgmntsze' [default: 300000;\n\
                    increased by np * 20000 for the MPI version with np processors; maximal\n\
		    value: 750000]. The entire input sequence is tiled with overlapping segments,\n\
                    thus all significant alignments will be reported independent of the setting.\n\
                    Decrease or increase the segment size depending on how much memory is available.\n\
  -x wsize     :   Specify the word size for the initial exact match search step.\n\
                    Default: 'wsize' = 12 to 16, depending on the species.  Must\n\
                    be set to >= 12.  Increasing 'wsize' will improve search\n\
                    speed and selectivity but reduce sensitivity.\n\
  -y minqHSP   :   Specify the minimum quality value of alignment HSPs to be pursued.\n\
                    Default: 'minqHSP' = 12 to 16, depending on the species.  Must\n\
                    be set to >= 'wsize'.  Increasing 'minqHSP' will improve search\n\
                    speed and selectivity but reduce sensitivity.\n\
  -z minqHSPc  :   Specify the minimum quality value of alignment HSP-chains to be pursued.\n\
                    Default: 'minqHSPc' = 30 to 40, depending on the species.  Must\n\
                    be set to >= 'minqHSP'.  Increasing 'minHSPc' will improve search\n\
                    speed and selectivity but reduce sensitivity.\n\
  -w minESTc   :   Specify the minimum EST coverage for an alignment HSP-chain to be pursued.\n\
                    Default: 'minESTc' = 0.0 [minimum value; maximum value: 1.0]\n\
  -p prmfile   :   Read parameters for EST alignments from file 'prmfile'.\n\
  -qQ qpdbn(s) :   Read target protein sequence data from library file(s)\n\
                    'qpdbn(s)' (FASTA-format; '-q': >gi|name; '-Q': >name).\n\
  -I matname   :   Read amino acid substitution scoring matrix from file\n\
                    'matname' [default: PBLO62].\n"
#define USAGE2 "\
  -a from      :   Analyze genomic sequence from position 'from'.\n\
  -b to        :   Analyze genomic sequence up to position 'to'.\n\
  -f           :   Analyze forward strand.\n\
  -r           :   Analyze reverse strand.\n\
  -R           :   Analyze both strands [default].\n\
  -oO outfname :   Redirect sorted output to file 'outfname' [default: stdout;\n\
                    if 'outfname' is set with the -o option, then on-the-fly output\n\
		    will be directed to stdout].\n\
  -c fcountA   :   If multiple genomic input sequences are supplied with the\n\
                    -lLg options, skip the first sequences and process starting from\n\
                    sequence number 'fcountA'.\n\
  -C fcountB   :   If multiple genomic input sequences are supplied with the\n\
                    -lLg options, skip any sequences after the sequence numbered\n\
                    'fcountB'.\n\
  -lL libfname :   Read (multiple) sequence data from library file 'libfname'\n\
                    (FASTA-format; '-l': >gi|name; '-L': >name).\n\
  -g gbfname(s):   Read nucleic acid sequence data from GenBank file(s)\n\
                    'gbfname(s)'.  If specified, the -g option MUST BE LAST.\n\n"
#else
#define USAGE1 "Usage:\n\n\
    %s [-hv] -s species [-dD dbest(s)] [-eE estdbn] [-k] [-m maxnest]\n\
    [-M sgmntsze] [-x wsize] [-y minqHSP] [-z minqHSPc] [-w minESTc] [-p prmfile] [-qQ qpdbn(s)] [-I matname]\n\
    [-a from] [-b to] [-frR] [-oO outfname] [-lL libfname] [-g gbfname(s)]\n\n\
  -h           :   HTML output [default: simple text]\n\
  -v           :   Verbose output\n\
  -s species   :   Set species to select the most appropriate splice site models.\n\
                   This parameter must be specified. Options:\n\
                    'Arabidopsis', 'maize', 'generic'.\n\
  -dD dbest(s) :   Read EST sequence data from pre-processed library file(s)\n\
                    'dbest(s)' (FASTA-format; '-d': >gi|name; '-D': >name).\n\
  -eE estdbn   :   Read EST sequence data from library file 'estdbn'.\n\
                    (FASTA-format; '-e': >gi|name; '-E': >name).\n\
  -k           :   Assume fixed 5' to 3' transcript orientation of supplied EST sequences.\n\
  -m maxnest   :   Display at most 'maxnest' spliced EST alignments per genomic\n\
                    DNA input [default: 500].\n\
  -M sgmntsze  :   Process input sequence in segments of size 'sgmntsze' [default: 300000;\n\
                    increased by np * 20000 for the MPI version with np processors; maximal\n\
		    value: 750000]. The entire input sequence is tiled with overlapping segments,\n\
                    thus all significant alignments will be reported independent of the setting.\n\
                    Decrease or increase the segment size depending on how much memory is available.\n\
  -x wsize     :   Specify the word size for the initial exact match search step.\n\
                    Default: 'wsize' = 12 to 16, depending on the species.  Must\n\
                    be set to >= 12.  Increasing 'wsize' will improve search\n\
                    speed and selectivity but reduce sensitivity.\n\
  -y minqHSP   :   Specify the minimum quality value of alignment HSPs to be pursued.\n\
                    Default: 'minqHSP' = 12 to 16, depending on the species.  Must\n\
                    be set to >= 'wsize'.  Increasing 'minqHSP' will improve search\n\
                    speed and selectivity but reduce sensitivity.\n\
  -z minqHSPc  :   Specify the minimum quality value of alignment HSP-chains to be pursued.\n\
                    Default: 'minqHSPc' = 30 to 40, depending on the species.  Must\n\
                    be set to >= 'minqHSP'.  Increasing 'minHSPc' will improve search\n\
                    speed and selectivity but reduce sensitivity.\n\
  -w minESTc   :   Specify the minimum EST coverage for an alignment HSP-chain to be pursued.\n\
                    Default: 'minESTc' = 0.0 [minimum value; maximum value: 1.0]\n\
  -p prmfile   :   Read parameters for EST alignments from file 'prmfile'.\n\
  -qQ qpdbn(s) :   Read target protein sequence data from library file(s)\n\
                    'qpdbn(s)' (FASTA-format; '-q': >gi|name; '-Q': >name).\n\
  -I matname   :   Read amino acid substitution scoring matrix from file\n\
                    'matname' [default: PBLO62].\n"
#define USAGE2 "\
  -a from      :   Analyze genomic sequence from position 'from'.\n\
  -b to        :   Analyze genomic sequence up to position 'to'.\n\
  -f           :   Analyze forward strand.\n\
  -r           :   Analyze reverse strand.\n\
  -R           :   Analyze both strands [default].\n\
  -oO outfname :   Redirect sorted output to file 'outfname' [default: stdout;\n\
                    if 'outfname' is set with the -o option, then on-the-fly output\n\
		    will be directed to stdout].\n\
  -c fcountA   :   If multiple genomic input sequences are supplied with the\n\
                    -lLg options, skip the first sequences and process starting from\n\
                    sequence number 'fcountA'.\n\
  -C fcountB   :   If multiple genomic input sequences are supplied with the\n\
                    -lLg options, skip any sequences after the sequence numbered\n\
                    'fcountB'.\n\
  -lL libfname :   Read (multiple) sequence data from library file 'libfname'\n\
                    (FASTA-format; '-l': >gi|name; '-L': >name).\n\
  -g gbfname(s):   Read nucleic acid sequence data from GenBank file(s)\n\
                    'gbfname(s)'.  If specified, the -g option MUST BE LAST.\n\n"
#endif



#define DATETIME		/* Gives access to date and time functions */
#define MATHFUNC		/* Gives access to real math functions */
#define OSYSFUNC		/* Gives access to the operating system */
#define PRIOFUNC		/* Gives access to PROMULA I/O functions */
#define STRGFUNC		/* Gives access to string manipulation */
#include "platform.h"
#include "minmax.h"
#include "EstData.h"
#include "SufArray.h"
#include "RawTextFile.h"	/* Header, raw text file services */

#include "html.h"


#ifdef MPI

#include <mpi.h>
int MPIV_rank, MPIV_nprc;
int MPIV_source;
int MPIV_tag;
#ifdef MPITIMING
double MPIV_begin, MPIV_finish, MPIV_firststart, MPIV_start, MPIV_begec,
	MPIV_bsend, MPIV_asend, MPI_Wtime(void);
#endif
MPI_Status MPIV_status;
MPI_Datatype MPI_gca;

#endif


struct bcvct {
  int bcnt[11], ccnt[6];
  float bfrq[11], cfrq[6];
} tbcv;

#include "def.h"
#include "abc.h"
#include "smat.h"
#include "mytype.h"
#include "sahmt.h"
#include "sequence.h"
#include "space.h"
#include "result.h"

static int CHECKFLAG = 0;

#define MAXRANGE  101
#define MAXDBNUM  400

#define INPUT_FILE   1
#define EST_FILE     2
#define QP_FILE      3
int MinMatchLen      = 12;  /* ... word size for exact match search step */
int MinQualityHSP    = 12;  /* ... minimum quality value for alignment HSPs to be pursued */
int MinQualityCHAIN  = 30;  /* ... minimum quality value for alignment HSP-chains to be pursued */
float MinESTcoverage = 0.0; /* ... minimum EST coverage for an alignment HSP-chain to be pursued */

int pstyle = 1, estop = 0, enoffset = 4, qpop = 0, qnoffset = 4;
int prmop = 0, dbop = 0, oflag = 0;
int rflag = 0, bflag = 1, revestallowed = 1;
char estdbn[257], qpfname[257], prmfile[257], sfname[257], outfname[257];
FILE *prmfp, *matfp;
FILE *outfp;

int mflag = 0, minsc, maxsc, smat[23][23];


struct sprmtr sprm;


int fcount = 0, fcountA = 0, fcountB = 999999, gbfcount = 0;
char *gdna, *gdnaR;
int af[23], pf[9];
float aq[20][4], pq[9][4];
int ngenes;

float *logV[8]={NULL};
char  *global_algnmnt=NULL;
int   logPia=0; /* to pass info */

UBYTE *Tstring;
int DataSize;
char estName[257];
char suffixName[257];
char lcptreeName[257];
char indexName[257];
char dataName[257];


#ifdef DAPBM
#include "daPbm7.h"
#else
#include "daPll.h"
#endif

static int estdbnum = 0;
static char estdbfnames[MAXDBNUM][257];
#define MAKENAMES(s) \
	strcpy(dataName,(s)); \
        strcpy(indexName,(s)); \
        strcpy(suffixName,(s)); \
        strcpy(lcptreeName,(s)); \
        strcat(dataName,".dat"); \
        strcat(indexName,".ind"); \
        strcat(suffixName,".suf"); \
        strcat(lcptreeName,".tre")

int frompA, top, htmlop, ifwdth;
char gdnafname[257];
int SGMNTSHFT = MAXGLGTH - 3;
int SGMNTSZE = 6 * (MAXGLGTH - 3);
int INCRPP = 20000;

int getlps(int Selection, char detsz, char noffset, char *sfname, char *seq);
int getlns(int Selection, char detsz, char noffset, char *sfname, char *seq);
int getgbs(int Selection, char detsz, char *sfname, char *seq);

void doit(int numbp, int ia, int ib, int rflag, int maxnest);
void fatal_error(char *buf);

#ifdef MPI
void buf2gca(struct gcalgnmnt_buf *buf, struct gcalgnmnt *gca);
void gca2buf(struct gcalgnmnt *gca, struct gcalgnmnt_buf *buf);
#endif

int FROMPOS, TOPOS;

#define LOG(f) ((f<=0)?-10000.0:(log(f)))

static void ComputeLogValue(float *pd,float *pa,float *pdR,float *paR,int size)
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

static void free_gca_list(struct gcalgnmnt *gca)
{
  struct gcalgnmnt *tgca = gca;

  while (tgca != NULL) {
    tgca = gca->next;
    free_gca(gca);
    gca = tgca;
  }

}						/* end free_gca_list() */



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
	if (seq[j + 1] != 10) {
	  ++d_pcntS[phs];
	  ++d_pcnt[phs][(int)seq[j]][(int)seq[j + 1]];
	  if (seq[j + 2] != 10) {
	    ++t_pcntS[phs];
	    ++t_pcnt[phs][(int)seq[j]][(int)seq[j + 1]][(int)seq[j + 2]];
	    if (seq[j + 3] != 10) {
	      ++r_pcntS[phs];
	      ++r_pcnt[phs][(int)seq[j]][(int)seq[j + 1]][(int)seq[j + 2]][(int)seq[j + 3]];
	      if (seq[j + 4] != 10) {
		++p_pcntS[phs];
		++p_pcnt[phs][(int)seq[j]][(int)seq[j + 1]][(int)seq[j + 2]][(int)seq[j + 3]][(int)seq[j + 4]];
		if (seq[j + 5] != 10) {
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
	if (seq[crd[i][1] - 5] != 10) {
	  ++m_pcntS[phs];
	  ++m_pcnt[phs][(int)seq[crd[i][1] - 5]];
	  if (seq[crd[i][1] - 4] != 10) {
	    ++d_pcntS[phs];
	    ++d_pcnt[phs][(int)seq[crd[i][1] - 5]][(int)seq[crd[i][1] - 4]];
	    if (seq[crd[i][1] - 3] != 10) {
	      ++t_pcntS[phs];
	      ++t_pcnt[phs][(int)seq[crd[i][1] - 5]][(int)seq[crd[i][1] - 4]][(int)seq[crd[i][1] - 3]];
	      if (seq[crd[i][1] - 2] != 10) {
		++r_pcntS[phs];
		++r_pcnt[phs][(int)seq[crd[i][1] - 5]][(int)seq[crd[i][1] - 4]]
		  [(int)seq[crd[i][1] - 3]][(int)seq[crd[i][1] - 2]];
		if (seq[crd[i][1] - 1] != 10) {
		  ++p_pcntS[phs];
		  ++p_pcnt[phs][(int)seq[crd[i][1] - 5]][(int)seq[crd[i][1] - 4]]
		    [(int)seq[crd[i][1] - 3]][(int)seq[crd[i][1] - 2]][(int)seq[crd[i][1] - 1]];
		  if (cflag && seq[crd[i + 1][0] - 1] != 10) {
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
	if (seq[crd[i][1] - 4] != 10) {
	  ++m_pcntS[phs];
	  ++m_pcnt[phs][(int)seq[crd[i][1] - 4]];
	  if (seq[crd[i][1] - 3] != 10) {
	    ++d_pcntS[phs];
	    ++d_pcnt[phs][(int)seq[crd[i][1] - 4]][(int)seq[crd[i][1] - 3]];
	    if (seq[crd[i][1] - 2] != 10) {
	      ++t_pcntS[phs];
	      ++t_pcnt[phs][(int)seq[crd[i][1] - 4]][(int)seq[crd[i][1] - 3]][(int)seq[crd[i][1] - 2]];
	      if (seq[crd[i][1] - 1] != 10) {
		++r_pcntS[phs];
		++r_pcnt[phs][(int)seq[crd[i][1] - 4]][(int)seq[crd[i][1] - 3]]
		  [(int)seq[crd[i][1] - 2]][(int)seq[crd[i][1] - 1]];
		if (cflag && seq[crd[i + 1][0] - 1] != 10) {
		  ++p_pcntS[phs];
		  ++p_pcnt[phs][(int)seq[crd[i][1] - 4]][(int)seq[crd[i][1] - 3]]
		    [(int)seq[crd[i][1] - 2]][(int)seq[crd[i][1] - 1]]
		    [(int)seq[crd[i + 1][0] - 1]];
		  if (cflag && crd[i + 1][1] - crd[i + 1][0] > 0 && seq[crd[i + 1][0]] != 10) {
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
	if (seq[crd[i][1] - 3] != 10) {
	  ++m_pcntS[phs];
	  ++m_pcnt[phs][(int)seq[crd[i][1] - 3]];
	  if (seq[crd[i][1] - 2] != 10) {
	    ++d_pcntS[phs];
	    ++d_pcnt[phs][(int)seq[crd[i][1] - 3]][(int)seq[crd[i][1] - 2]];
	    if (seq[crd[i][1] - 1] != 10) {
	      ++t_pcntS[phs];
	      ++t_pcnt[phs][(int)seq[crd[i][1] - 3]][(int)seq[crd[i][1] - 2]][(int)seq[crd[i][1] - 1]];
	      if (cflag && seq[crd[i + 1][0] - 1] != 10) {
		++r_pcntS[phs];
		++r_pcnt[phs][(int)seq[crd[i][1] - 3]][(int)seq[crd[i][1] - 2]]
		  [(int)seq[crd[i][1] - 1]][(int)seq[crd[i + 1][0] - 1]];
		if (cflag && crd[i + 1][1] - crd[i + 1][0] > 0 && seq[crd[i + 1][0]] != 10) {
		  ++p_pcntS[phs];
		  ++p_pcnt[phs][(int)seq[crd[i][1] - 3]][(int)seq[crd[i][1] - 2]]
		    [(int)seq[crd[i][1] - 1]][(int)seq[crd[i + 1][0] - 1]]
		    [(int)seq[crd[i + 1][0]]];
		  if (cflag && crd[i + 1][1] - crd[i + 1][0] > 1 && seq[crd[i + 1][0] + 1] != 10) {
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
	if (seq[crd[i][1] - 2] != 10) {
	  ++m_pcntS[phs];
	  ++m_pcnt[phs][(int)seq[crd[i][1] - 2]];
	  if (seq[crd[i][1] - 1] != 10) {
	    ++d_pcntS[phs];
	    ++d_pcnt[phs][(int)seq[crd[i][1] - 2]][(int)seq[crd[i][1] - 1]];
	    if (cflag && seq[crd[i + 1][0] - 1] != 10) {
	      ++t_pcntS[phs];
	      ++t_pcnt[phs][(int)seq[crd[i][1] - 2]][(int)seq[crd[i][1] - 1]]
		[(int)seq[crd[i + 1][0] - 1]];
	      if (cflag && crd[i + 1][1] - crd[i + 1][0] > 0 && seq[crd[i + 1][0]] != 10) {
		++r_pcntS[phs];
		++r_pcnt[phs][(int)seq[crd[i][1] - 2]][(int)seq[crd[i][1] - 1]]
		  [(int)seq[crd[i + 1][0] - 1]][(int)seq[crd[i + 1][0]]];
		if (cflag && crd[i + 1][1] - crd[i + 1][0] > 1 && seq[crd[i + 1][0] + 1] != 10) {
		  ++p_pcntS[phs];
		  ++p_pcnt[phs][(int)seq[crd[i][1] - 2]][(int)seq[crd[i][1] - 1]]
		    [(int)seq[crd[i + 1][0] - 1]][(int)seq[crd[i + 1][0]]]
		    [(int)seq[crd[i + 1][0] + 1]];
		  if (cflag && crd[i + 1][1] - crd[i + 1][0] > 2 && seq[crd[i + 1][0] + 2] != 10) {
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
      if (seq[crd[i][1] - 1] != 10) {
	++m_pcntS[phs];
	++m_pcnt[phs][(int)seq[crd[i][1] - 1]];
	if (cflag && seq[crd[i + 1][0] - 1] != 10) {
	  ++d_pcntS[phs];
	  ++d_pcnt[phs][(int)seq[crd[i][1] - 1]][(int)seq[crd[i + 1][0] - 1]];
	  if (cflag && crd[i + 1][1] - crd[i + 1][0] > 0 && seq[crd[i + 1][0]] != 10) {
	    ++t_pcntS[phs];
	    ++t_pcnt[phs][(int)seq[crd[i][1] - 1]][(int)seq[crd[i + 1][0] - 1]][(int)seq[crd[i + 1][0]]];
	    if (cflag && crd[i + 1][1] - crd[i + 1][0] > 1 && seq[crd[i + 1][0] + 1] != 10) {
	      ++r_pcntS[phs];
	      ++r_pcnt[phs][(int)seq[crd[i][1] - 1]][(int)seq[crd[i + 1][0] - 1]]
		[(int)seq[crd[i + 1][0]]][(int)seq[crd[i + 1][0] + 1]];
	      if (cflag && crd[i + 1][1] - crd[i + 1][0] > 2 && seq[crd[i + 1][0] + 2] != 10) {
		++p_pcntS[phs];
		++p_pcnt[phs][(int)seq[crd[i][1] - 1]][(int)seq[crd[i + 1][0] - 1]]
		  [(int)seq[crd[i + 1][0]]][(int)seq[crd[i + 1][0] + 1]]
		  [(int)seq[crd[i + 1][0] + 2]];
		if (cflag && crd[i + 1][1] - crd[i + 1][0] > 3 && seq[crd[i + 1][0] + 3] != 10) {
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
      if (crd[i][1] - crd[i][0] > 3 && seq[crd[i][1] - 5] != 10) {
	++m_pcntS[phs];
	++m_pcnt[phs][(int)seq[crd[i][1] - 5]];
	if (seq[crd[i][1] - 4] != 10) {
	  ++d_pcntS[phs];
	  ++d_pcnt[phs][(int)seq[crd[i][1] - 5]][(int)seq[crd[i][1] - 4]];
	  if (seq[crd[i][1] - 3] != 10) {
	    ++t_pcntS[phs];
	    ++t_pcnt[phs][(int)seq[crd[i][1] - 5]][(int)seq[crd[i][1] - 4]][(int)seq[crd[i][1] - 3]];
	    if (seq[crd[i][1] - 2] != 10) {
	      ++r_pcntS[phs];
	      ++r_pcnt[phs][(int)seq[crd[i][1] - 5]][(int)seq[crd[i][1] - 4]][(int)seq[crd[i][1] - 3]]
		[(int)seq[crd[i][1] - 2]];
	      if (seq[crd[i][1] - 1] != 10) {
		++p_pcntS[phs];
		++p_pcnt[phs][(int)seq[crd[i][1] - 5]][(int)seq[crd[i][1] - 4]]
		  [(int)seq[crd[i][1] - 3]][(int)seq[crd[i][1] - 2]][(int)seq[crd[i][1] - 1]];
	      }
	    }
	  }
	}
      }
      phs = (phs + 1) % 3;
      if (crd[i][1] - crd[i][0] > 2 && seq[crd[i][1] - 4] != 10) {
	++m_pcntS[phs];
	++m_pcnt[phs][(int)seq[crd[i][1] - 4]];
	if (seq[crd[i][1] - 3] != 10) {
	  ++d_pcntS[phs];
	  ++d_pcnt[phs][(int)seq[crd[i][1] - 4]][(int)seq[crd[i][1] - 3]];
	  if (seq[crd[i][1] - 2] != 10) {
	    ++t_pcntS[phs];
	    ++t_pcnt[phs][(int)seq[crd[i][1] - 4]][(int)seq[crd[i][1] - 3]][(int)seq[crd[i][1] - 2]];
	    if (seq[crd[i][1] - 1] != 10) {
	      ++r_pcntS[phs];
	      ++r_pcnt[phs][(int)seq[crd[i][1] - 4]][(int)seq[crd[i][1] - 3]][(int)seq[crd[i][1] - 2]]
		[(int)seq[crd[i][1] - 1]];
	    }
	  }
	}
      }
      phs = (phs + 1) % 3;
      if (crd[i][1] - crd[i][0] > 1 && seq[crd[i][1] - 3] != 10) {
	++m_pcntS[phs];
	++m_pcnt[phs][(int)seq[crd[i][1] - 3]];
	if (seq[crd[i][1] - 2] != 10) {
	  ++d_pcntS[phs];
	  ++d_pcnt[phs][(int)seq[crd[i][1] - 3]][(int)seq[crd[i][1] - 2]];
	  if (seq[crd[i][1] - 1] != 10) {
	    ++t_pcntS[phs];
	    ++t_pcnt[phs][(int)seq[crd[i][1] - 3]][(int)seq[crd[i][1] - 2]][(int)seq[crd[i][1] - 1]];
	  }
	}
      }
      phs = (phs + 1) % 3;
      if (crd[i][1] - crd[i][0] > 0 && seq[crd[i][1] - 2] != 10) {
	++m_pcntS[phs];
	++m_pcnt[phs][(int)seq[crd[i][1] - 2]];
	if (seq[crd[i][1] - 1] != 10) {
	  ++d_pcntS[phs];
	  ++d_pcnt[phs][(int)seq[crd[i][1] - 2]][(int)seq[crd[i][1] - 1]];
	}
      }
      phs = (phs + 1) % 3;
      if (seq[crd[i][1] - 1] != 10) {
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
  int slen=strlen(str);
  int nextFree=slen;
  int tmplen;

  for (i = ia; i <= ib + ni; ++i) {
    if (seq[i] == 24)
      ++ni;
    if ((i - ni - ia) % 10 == 0 && (i - ni - ia) % 60 != 0 && seq[i] != 24)
      str[nextFree++]=' ';
    if ((i - ni - ia) % 60 == 0 && seq[i] != 24) {
      if (nflag) {
	tmplen = sprintf(&(str[nextFree]), "\n%*d  ", ifwdth,
				(rflag)?(numbp - i + ni):(i - ni + 1));
      nextFree+=tmplen;
      }
      else {
	str[nextFree++]='\n';
      }
    }
    str[nextFree++]=ABC[(int)seq[i]];
  }
  str[nextFree++]='\n';
  str[nextFree++]='\0';

}				/* end prt_seq_segment_to_str() */



void prt_segment_trl(FILE *fp, int seq[], int numbp, int ia, int ib, int rflag,
		     int phsA, int cl, char ABC[])
	/* prints specified segment of seq[] in ABC translation */
{
  int i, rc, iaa = ia;

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
	fprintf(fp, "%c", ABC[codtoaa[(int)seq[i - 1]][(int)seq[i]][(int)seq[i + 1]]]);
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



static void free_gpa_list(struct gpalgnmnt *gpa)
{
  struct gpalgnmnt *tgpa = gpa;

  while (tgpa != NULL) {
    tgpa = gpa->next;
    free_gpa(gpa);
    gpa = tgpa;
  }

}						/* end free_gpa_list() */


int NALLOC_GCA = 0, NDEALLOC_GCA = 0;
int NALLOC_CDNA = 0, NDEALLOC_CDNA = 0;

int main(int argc, char *argv[])
{
  int i, j, k, numbp, anum = 0, frmterr, maxnest = 500, fromp, frompp, topp;
  int libop = 0, lnoffset = 4, gbflag = 0;
  int fposset = 0, tposset = 0, tmp;
  static char dfname[257];
#ifdef HTMLWS
  static char imgName[257];
#endif
  time_t tlc;
  int pmt[23], tsmat[23][23];
  char ch, ch2[2], aa[23], matname[257];
#ifdef MPI
  static char MPIV_outputfile[257];
  void MPI_derived_type_gca(struct gcalgnmnt_buf *gca,MPI_Datatype *MPI_gca);
  static struct gcalgnmnt_buf gca;
#endif

  sprm.pdg = 0.03F;
  sprm.ids = 2.0F;
  sprm.mms = -2.0F;
  sprm.nns = 0.0F;
  sprm.dls = -5.0F;
  sprm.unfrm = 0;
  sprm.min_intron_length = 30;
  sprm.min_exon_length = 5;
  sprm.min_nbr_endmatches = 2;
  sprm.tiny_exon = 20;
  sprm.short_exon = 200;
  sprm.long_intron = 300;
  sprm.poor_exon_score = 0.70;
  sprm.poor_donor_score = 0.5;
  sprm.poor_acptr_score = 0.5;
  sprm.join_length = 300;
  outfp = stdout;
  time(&tlc);

  if (argc < 2) {
    fprintf(stderr, USAGE1, argv[0]);
    fprintf(stderr, USAGE2);
    exit(0);
  }

#ifdef MPI
  MPI_Init(&argc,&argv);                    /* start MPI */
  MPI_Comm_rank(MPI_COMM_WORLD,&MPIV_rank); /* processor rank */
  MPI_Comm_size(MPI_COMM_WORLD,&MPIV_nprc); /* number of processors */
#endif
 
  for (i = 1; i < argc; ++i) {
    if (*argv[i] == '-')
      switch (*(argv[i] + 1)) {
      case 'h':
        htmlop = 1;
#ifdef HTMLWS
        SGMNTSZE = 60000;
#endif
        ++anum;
        break;
      case 'v':
	pstyle = 2;
	++anum;
	break;
      case 's':
#ifdef DAPBM
	if (strcmp("human", argv[i + 1]) == 0)
	  {Species =  0; MinMatchLen = 16; MinQualityCHAIN = 40;}
	else if (strcmp("mouse", argv[i + 1]) == 0)
	  {Species =  1; MinMatchLen = 16; MinQualityCHAIN = 40;}
	else if (strcmp("rat", argv[i + 1]) == 0)
	  {Species =  2; MinMatchLen = 16; MinQualityCHAIN = 40;}
	else if (strcmp("chicken", argv[i + 1]) == 0)
	  {Species =  3; MinMatchLen = 16; MinQualityCHAIN = 40;}
	else if (strcmp("Drosophila", argv[i + 1]) == 0)
	  {Species =  4; MinMatchLen = 16; MinQualityCHAIN = 40;}
	else if (strcmp("nematode", argv[i + 1]) == 0)
	  {Species =  5; MinMatchLen = 16; MinQualityCHAIN = 40;}
	else if (strcmp("yeast", argv[i + 1]) == 0)
	  {Species =  6; MinMatchLen = 16; MinQualityCHAIN = 40;}
	else if (strcmp("Aspergillus", argv[i + 1]) == 0)
	  {Species =  7; MinMatchLen = 16; MinQualityCHAIN = 40;}
	else if (strcmp("Arabidopsis", argv[i + 1]) == 0)
	  {Species =  8; MinMatchLen = 12; MinQualityCHAIN = 30;}
	else if (strcmp("maize", argv[i + 1]) == 0)
	  {Species =  9; MinMatchLen = 12; MinQualityCHAIN = 30;}
	else if (strcmp("rice", argv[i + 1]) == 0)
	  {Species = 10; MinMatchLen = 12; MinQualityCHAIN = 30;}
	else if (strcmp("Medicago", argv[i + 1]) == 0)
	  {Species = 11; MinMatchLen = 12; MinQualityCHAIN = 30;}
	else if (strcmp("Populus", argv[i + 1]) == 0)
	  {Species = 12; MinMatchLen = 12; MinQualityCHAIN = 30;}
	else if (strcmp("Daphnia", argv[i + 1]) == 0)
	  {Species = 13; MinMatchLen = 16; MinQualityCHAIN = 40;}
	else if (strcmp("generic", argv[i + 1]) == 0)
	  {Species = 99; MinMatchLen = 16; MinQualityCHAIN = 40;}
#else
	if (strcmp("maize", argv[i + 1]) == 0)
	  {Species =  0; MinMatchLen = 12; MinQualityCHAIN = 30;}
	else if (strcmp("Arabidopsis", argv[i + 1]) == 0)
	  {Species =  1; MinMatchLen = 12; MinQualityCHAIN = 30;}
	else if (strcmp("generic", argv[i + 1]) == 0)
	  {Species =  2; MinMatchLen = 16; MinQualityCHAIN = 40;}
#endif
	else {
	  fprintf(stderr, "\nINVALID species %s. Exit.\n", argv[i + 1]);
	  exit(0);
	}
	MinQualityHSP = MinMatchLen;
	anum += 2;
	break;
      case 'd':
	if (dbop == 1) {
	  fprintf(stderr,"\nPlease use either the -d OR the -D option. Exit.\n");
	  exit(0);
	}
	dbop = 1;
	enoffset = 4;
	do {
	  strcpy(estdbfnames[estdbnum], argv[i + 1 + estdbnum]);
	  if (++estdbnum > MAXDBNUM) {
	    fprintf(stderr, "\nMaximal allowed number of EST library files (%d) exceeded. Exit.\n",MAXDBNUM);
	    exit(0);
	  }
	} while (((i + 1 + estdbnum) < (argc - 1)) &&
		 (*argv[i + 1 + estdbnum] != '-'));
	anum += (estdbnum + 1);
	strcpy(estdbn, estdbfnames[0]);
	break;
      case 'D':
	if (dbop == 1) {
	  fprintf(stderr,"\nPlease use either the -d OR the -D option. Exit.\n");
	  exit(0);
	}
	dbop = 1;
	enoffset = 1;
	do {
	  strcpy(estdbfnames[estdbnum], argv[i + 1 + estdbnum]);
	  if (++estdbnum > MAXDBNUM) {
	    fprintf(stderr, "\nMaximal allowed number of EST library files (%d) exceeded. Exit.\n",MAXDBNUM);
	    exit(0);
	  }
	} while (((i + 1 + estdbnum) < (argc - 1)) &&
		 (*argv[i + 1 + estdbnum] != '-'));
	anum += (estdbnum + 1);
	strcpy(estdbn, estdbfnames[0]);
	break;
      case 'e':
	if (estop == 1) {
	  fprintf(stderr,"\nPlease use either the -e OR the -E option. Exit.\n");
	  exit(0);
	}
	estop = 1;
	enoffset = 4;
	strcpy(estdbn, argv[i + 1]);
	strcpy(estdbfnames[0], estdbn);
	estdbnum = 1;
	anum += 2;
	break;
      case 'E':
	if (estop == 1) {
	  fprintf(stderr,"\nPlease use either the -e OR the -E option. Exit.\n");
	  exit(0);
	}
	estop = 1;
	enoffset = 1;
	strcpy(estdbn, argv[i + 1]);
	strcpy(estdbfnames[0], estdbn);
	estdbnum = 1;
	anum += 2;
	break;
      case 'k':
	revestallowed = 0;
	++anum;
	break;
      case 'm':
	maxnest = atoi(argv[i + 1]);
	anum += 2;
	break;
      case 'M':
	SGMNTSZE = atoi(argv[i + 1]);
	if (SGMNTSZE < 3 * SGMNTSHFT) SGMNTSZE = 3 * SGMNTSHFT;
	anum += 2;
	break;
      case 'x':
	MinMatchLen = atoi(argv[i+1]);
	anum += 2;
	break;
      case 'y':
	MinQualityHSP = atoi(argv[i+1]);
	anum +=2;
	break;
      case 'z':
	MinQualityCHAIN = atoi(argv[i+1]);
	anum +=2;
	break;
      case 'w':
	MinESTcoverage = atof(argv[i+1]);
	anum +=2;
	break;
      case 'p':
	prmop = 1;
	strcpy(prmfile, argv[i + 1]);
	anum += 2;
	break;
      case 'q':
	qpop = 1;
	qnoffset = 4;
	strcpy(qpfname, argv[i + 1]);
	anum += 2;
	break;
      case 'Q':
	qpop = 1;
	qnoffset = 1;
	strcpy(qpfname, argv[i + 1]);
	anum += 2;
	break;
      case 'I':
	mflag = 1;
	strcpy(matname, argv[i + 1]);
	anum += 2;
	break;
      case 'a':
	FROMPOS = atoi(argv[i + 1]) - 1;
	fposset = 1;
	anum += 2;
	break;
      case 'b':
	TOPOS = atoi(argv[i + 1]) - 1;
	tposset = 1;
	anum += 2;
	break;
      case 'f':
	bflag = 0;
	++anum;
	break;
      case 'r':
	bflag = 2;
	++anum;
	break;
      case 'R':
	bflag = 1;
	++anum;
	break;
      case 'o':
	oflag = 1;
	strcpy(outfname, argv[i + 1]);
#ifdef MPI
	if (MPIV_rank == 0) {
#endif
	if ((outfp = fopen(outfname, "w")) == NULL) {
	  fprintf(stderr, "File %s cannot be opened.\n", outfname);
	  perror(outfname);
	  exit(-1);
	}
#ifdef MPI
	}
#endif
	anum += 2;
	break;
      case 'O':
	oflag = 2;
	strcpy(outfname, argv[i + 1]);
#ifdef MPI
	if (MPIV_rank == 0) {
#endif
	if ((outfp = fopen(outfname, "w")) == NULL) {
	  fprintf(stderr, "File %s cannot be opened.\n", outfname);
	  perror(outfname);
	  exit(-1);
	}
#ifdef MPI
	}
#endif
	anum += 2;
	break;
      case 'c':
	fcountA = atoi(argv[i + 1]);
	anum += 2;
	break;
      case 'C':
	fcountB = atoi(argv[i + 1]);
	anum += 2;
	break;
      case 'l':
	libop = 1;
	lnoffset = 4;
	strcpy(gdnafname, argv[i + 1]);
	anum += 2;
	break;
      case 'L':
	libop = 1;
	lnoffset = 1;
	strcpy(gdnafname, argv[i + 1]);
	anum += 2;
	break;
      case 'g':
	gbflag = 1;
	++anum;
	break;
      default:
	break;
      }
  }

#ifdef MPI
   if (MPIV_rank > 0) {
     sprintf(MPIV_outputfile,"%s_MPI_logfile_prc%d",outfname,MPIV_rank);
     if ( (outfp= fopen(MPIV_outputfile,"w")) == NULL ) {
       fprintf(stderr,"File %s cannot be opened.\n",MPIV_outputfile);
       perror(MPIV_outputfile); exit(-1);
     }
   }
#ifdef MPITIMING
   MPIV_begin = MPI_Wtime();
#endif
   MPI_derived_type_gca(&gca,&MPI_gca);
   SGMNTSZE += MPIV_nprc * INCRPP;
   if (SGMNTSZE > MAXSGMNTSZE) SGMNTSZE = MAXSGMNTSZE;
#endif

  if (Species < 0) {
    fprintf(stderr,
	    "\nYou must specify a species with the -s option to select the most appropriate\n\
               splice site models for your input.\n\n");
    fprintf(stderr, USAGE1, argv[0]);
    fprintf(stderr, USAGE2);
    exit(0);
  }
  if (!(libop + gbflag)) {
    fprintf(stderr,
	    "\nPlease specify an input file by either the -l or the -g option.\n\n");
    fprintf(stderr, USAGE1, argv[0]);
    fprintf(stderr, USAGE2);
    exit(0);
  }
  if (dbop == 0 && estop == 0 && qpop == 0) {
    fprintf(stderr,
	    "\nNothing to be done. Please specify at least one of the -dD, -eE, or -qQ options.\n\n");
    fprintf(stderr, USAGE1, argv[0]);
    fprintf(stderr, USAGE2);
    exit(0);
  }
  if (dbop == 1 && estop == 1) {
    fprintf(stderr, "\nPlease use either the -d OR the -e option.\n");
    exit(0);
  }

#ifdef HTMLWS
  if (htmlop) {
    if (oflag) {
      strcpy(imgName,outfname);
      strcat(imgName,"_img");
    }
    else {
      strcpy(imgName,"std_img");
    }
    if ((imageDataFh = fopen(imgName, "w")) == NULL) {
      fprintf(stderr, "File %s cannot be opened.\n", imgName);
      perror(imgName);
      exit(-1);
    }
  }
#endif

  fprintf(stderr, "\nNOW EXECUTING:  ");
  for (i = 0; i < argc; ++i)
    fprintf(stderr, " %s", argv[i]);
  fprintf(stderr, "\n\n");

  if (htmlop) {
    ADD_HEAD;
    fprintf(outfp, "<A NAME=\"TOP\"></A>\n");
  }

  fprintf(outfp, "GeneSeqer.   Version of February 11, 2019.\n");
  fprintf(outfp, "Date run: %s\n", ctime(&tlc));
#ifdef DAPBM
  if (Species < NMDLS)
    fprintf(outfp, "(Bayesian) Splice site model (species):\t%s\n",
	  name_model[Species]);
  else
    fprintf(outfp, "(Bayesian) Splice site model (species):\tgeneric\n");
#else
  if (Species < NMDLS)
    fprintf(outfp, "(Logitlinear) Splice site model (species):\t%s\n",
	  species_name[Species]);
  else
    fprintf(outfp, "(Logitlinear) Splice site model (species):\tgeneric\n");
#endif
  if (MinMatchLen < 12) MinMatchLen = 12;
  if (MinQualityHSP < MinMatchLen) MinQualityHSP = MinMatchLen;
  if (MinQualityCHAIN < MinQualityHSP) MinQualityCHAIN = MinQualityHSP;
  if (dbop || estop) {
    fprintf(outfp, "Fast search parameters:  (-x) wsize %d, (-y) minqHSP %d, (-z) minqHSPc %d.\n",
		MinMatchLen, MinQualityHSP, MinQualityCHAIN);
    fprintf(outfp, "                         (-w) minESTc %4.2f.\n", MinESTcoverage);
  }

  if (oflag == 1) {
    if (htmlop)
      fprintf(stdout, "<A NAME=\"TOP\"></A>\n");
    fprintf(stdout, "GeneSeqer.   Version of February 11, 2019.\n");
    fprintf(stdout, "Date run: %s\n", ctime(&tlc));
#ifdef DAPBM
    if (Species < NMDLS)
      fprintf(stdout, "(Bayesian) Splice site model (species):\t%s\n",
		name_model[Species]);
    else
      fprintf(stdout, "(Bayesian) Splice site model (species):\tgeneric\n");
#else
    if (Species < NMDLS)
      fprintf(stdout, "(Logitlinear) Splice site model (species):\t%s\n",
		species_name[Species]);
    else
      fprintf(stdout, "(Logitlinear) Splice site model (species):\tgeneric\n");
#endif
    if (dbop || estop) {
      fprintf(stdout, "Fast search parameters:  (-x) wsize %d, (-y) minqHSP %d, (-z) minqHSPc %d.\n",
		MinMatchLen, MinQualityHSP, MinQualityCHAIN);
      fprintf(stdout, "                         (-w) minESTc %4.2f.\n", MinESTcoverage);
    }
  }

  rflag = 0;

  if (estop) {
    MAKENAMES(estdbn);
    EstCreateDatIndFiles(estdbn, dataName, indexName);
    DataSize = EstReadDatFile(dataName, &Tstring);
    if (DataSize == 0) {
	printf("\nSorry. The specified EST database is empty.\n\n");
	delfile(dataName);
	delfile(indexName);
	if (!qpop){
          if (htmlop) ADD_TAIL;
          exit(-1);
        }
    }
    else {
	SuffixArrayCreate(DataSize);
	SuffixArraySort();
	SuffixArraySave(suffixName, lcptreeName);
	SuffixArrayDestroy();
	EstDestroy();
	ResultCreate(indexName);
    }
  }

  if (prmop) {
    if ((prmfp = fopen(prmfile, "r")) == NULL) {
      fprintf(stderr, "File %s cannot be opened.\n", prmfile);
      perror(prmfile);
      if (htmlop) ADD_TAIL;
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
    fprintf(outfp, "\n\tMIN_NBR_ENDMATCHES	%3d", sprm.min_nbr_endmatches);
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
      if (htmlop) ADD_TAIL;
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
    }
    else {
      if ((matfp = fopen(matname, "r")) == NULL) {
        fprintf(stderr, "File %s cannot be opened.\n", matname);
        perror(matname);
        if (htmlop) ADD_TAIL;
        exit(-1);
      }
      fprintf(stderr, "Read scoring matrix file %s.\n", matname);

      i = 0;
      while ((ch = (char) (fgetc(matfp))) != '\n')
        if (ch > 62 && ch < 91)
	  aa[i++] = ch;
      if (i != 23) {
        fprintf(outfp, "\nERROR: Matrix %s has less than 23 column labels.\n", matname);
        if (htmlop) ADD_TAIL;
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
        if (htmlop) ADD_TAIL;
        exit(0);
      }
    }/* end  of else */

    if (pstyle == 2 || mflag) {
      fprintf(outfp, "Amino acid substitution scoring matrix: %s \n", matname);
      if (oflag == 1)
	fprintf(stdout, "Amino acid substitution scoring matrix: %s \n", matname);
    }
/*                                                                            */
/******************************************************************************/
  }

  if (libop) {
    if (openRawTextFile(INPUT_FILE, gdnafname, 1) != 0) {
      fprintf(stderr, "File %s cannot be opened.\n", gdnafname);
      if (htmlop) ADD_TAIL;
      exit(-1);
    }
    fprintf(stderr, "Library file %s has been opened for reading.\n", gdnafname);
    fprintf(outfp, "Input file  : %s \n", gdnafname);
    frmterr = 1;
    while ((numbp = getlns(INPUT_FILE, 1, lnoffset, sfname, gdna)) != 0) {
      if (numbp == -1)
	continue;
      if ((gdna = (char *) calloc(numbp, sizeof(char))) == NULL)
	fatal_error("Error: memory allocation failed. Exit.\n");
      numbp = getlns(INPUT_FILE, 0, lnoffset, sfname, gdna);
      frmterr = 0;
      if (++fcount % 100 == 0)
	fprintf(stderr,
		" ... now processing nucleotide sequence %4d (%s)\n", fcount, sfname);
      if (fcount < fcountA  ||  fcount > fcountB) {
	free(gdna);
	continue;
      }
      if (numbp < 1000000) ifwdth = 6; else ifwdth = 10;
      if ((gdnaR = (char *) calloc(numbp, sizeof(char))) == NULL)
	fatal_error("Error: memory allocation failed. Exit.\n");
      complement_seq(gdna, numbp, gdnaR);
      if (CHECKFLAG) {
	fprintf(outfp, "\n");
	for (i = 0; i < 80; ++i)
	  fprintf(outfp, "*");
        if (htmlop)
	  fprintf(outfp, "\nSequence %4d (File: %s)\n", fcount, getDnaLink(1,sfname));
        else
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
      strcpy(dfname, sfname);
      if (htmlop && (dbop || (estop && (DataSize > 0)))) {
	for (fromp= frompp= FROMPOS; fromp <= TOPOS; fromp += (SGMNTSZE - SGMNTSHFT)) {
	  topp = (frompp + SGMNTSZE - 1 > TOPOS) ? TOPOS : frompp + SGMNTSZE - 1;
	  frompA = fromp;
	  if (fromp > FROMPOS) frompA += SGMNTSHFT;
	  top = (fromp + SGMNTSZE - 1 > TOPOS) ? TOPOS : fromp + SGMNTSZE - 1;
	  if (fromp > FROMPOS)
	    fprintf(outfp, "\n<A HREF=\"#BOTTOM-PGL-%s-%d-%d\">Scroll down to \"Segment %d-%d\"</A>",dfname,frompp,topp,frompA+1,top+1);
	  frompp = frompA;
	  if (top == TOPOS)
	    break;
	}
	fprintf(outfp,"\n");
      }
      for (fromp = FROMPOS; fromp <= TOPOS; fromp += (SGMNTSZE - SGMNTSHFT)) {
	frompA = fromp;
	if (fromp > FROMPOS) frompA += SGMNTSHFT;
	top = (fromp + SGMNTSZE - 1 > TOPOS) ? TOPOS : fromp + SGMNTSZE - 1;
	strcpy(sfname, dfname);
	rflag = 0;
#ifdef MPI
        MPI_Barrier(MPI_COMM_WORLD);
#endif
	doit(numbp, fromp, top, rflag, maxnest);
	rflag = 1;
	strcpy(sfname, dfname);
	if (CHECKFLAG && fromp == FROMPOS) {
	  fprintf(outfp, "\n");
	  for (i = 0; i < 80; ++i)
	    fprintf(outfp, "*");
	  if (htmlop)
            fprintf(outfp, "\nSequence %4d (File: %s)\n", fcount, getDnaLink(1,sfname));
	  else
            fprintf(outfp, "\nSequence %4d (File: %s)\n", fcount, sfname);
	  fprintf(outfp, "Reverse strand.\n");
	  prt_seq_segment(outfp, gdnaR, numbp, 0, numbp - 1, rflag, NAUC, 1);
	}
#ifdef MPI
        MPI_Barrier(MPI_COMM_WORLD);
#endif
	doit(numbp, fromp, top, rflag, maxnest);
	if (top == TOPOS)
	  break;
      }
      free(gdna);
      free(gdnaR);
      fflush(outfp);
    }
    closeRawTextFile(INPUT_FILE);
    if (frmterr) {
      fprintf(outfp,
	      "\nNo FASTA-formatted input sequence data found in file %s.\n",
	      gdnafname);
      if (oflag == 1)
        fprintf(stdout,
	      "\nNo FASTA-formatted input sequence data found in file %s.\n",
	      gdnafname);
    }
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
      strcpy(gdnafname, argv[anum]);
      if (anum + 1 == argc && gbfcount == 1) {
	fprintf(outfp, "GenBank file: %s\n", gdnafname);
	if (oflag == 1)
	  fprintf(stdout, "GenBank file: %s\n", gdnafname);
      }
      else {
	fprintf(outfp, "\n\n");
	for (i = 0; i < 80; ++i)
	  fprintf(outfp, "*");
	fprintf(outfp, "\nGenBank file: %s (# %4d)\n", gdnafname, gbfcount);
	if (oflag == 1) {
	  fprintf(stdout, "\n\n");
	  for (i = 0; i < 80; ++i)
	    fprintf(stdout, "*");
	  fprintf(stdout, "\nGenBank file: %s (# %4d)\n", gdnafname, gbfcount);
	}
      }
      while ((numbp = getgbs(INPUT_FILE, 1, sfname, gdna)) != 0) {
	if (numbp == -1)
	  continue;
	if ((gdna = (char *) calloc(numbp, sizeof(char))) == NULL) {
	  fatal_error("Error: memory allocation failed. Exit.\n");
	}
	numbp = getgbs(INPUT_FILE, 0, sfname, gdna);
	frmterr = 0;
	if (++fcount % 100 == 0) {
	  fprintf(stderr, " ... now processing nucleotide sequence %4d (%s)\n",
		  fcount, sfname);
	}
        if (fcount < fcountA  ||  fcount > fcountB) {
	  free(gdna);
	  continue;
	}
        if (numbp < 1000000) ifwdth = 6; else ifwdth = 10;
	time(&tlc);
	fprintf(stderr, "\nSeq #%3d = %s     %s", fcount, sfname, ctime(&tlc));
	if ((gdnaR = (char *) calloc(numbp, sizeof(char))) == NULL)
	  fatal_error("Error: memory allocation failed. Exit.\n");
	complement_seq(gdna, numbp, gdnaR);
	if (CHECKFLAG) {
	  fprintf(outfp, "\n");
	  for (i = 0; i < 80; ++i)
	    fprintf(outfp, "*");
	  if (htmlop)
            fprintf(outfp, "\nSequence %4d (File: %s)\n", fcount, getDnaLink(1,sfname));
	  else
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
	strcpy(dfname, sfname);
	if (htmlop && (dbop || (estop && (DataSize > 0)))) {
	  for (fromp= frompp= FROMPOS; fromp <= TOPOS; fromp += (SGMNTSZE - SGMNTSHFT)) {
	    topp = (frompp + SGMNTSZE - 1 > TOPOS) ? TOPOS : frompp + SGMNTSZE - 1;
	    frompA = fromp;
	    if (fromp > FROMPOS) frompA += SGMNTSHFT;
	    top = (fromp + SGMNTSZE - 1 > TOPOS) ? TOPOS : fromp + SGMNTSZE - 1;
	    if (fromp > FROMPOS)
	      fprintf(outfp, "\n<A HREF=\"#BOTTOM-PGL-%s-%d-%d\">Scroll down to \"Segment %d-%d\"</A>",dfname,frompp,topp,frompA+1,top+1);
	    frompp = frompA;
	    if (top == TOPOS)
	      break;
	  }
	  fprintf(outfp,"\n");
	}
	for (fromp = FROMPOS; fromp <= TOPOS; fromp += (SGMNTSZE - SGMNTSHFT)) {
	  frompA = fromp;
	  if (fromp > FROMPOS) frompA += SGMNTSHFT;
	  top = (fromp + SGMNTSZE - 1 > TOPOS) ? TOPOS : fromp + SGMNTSZE - 1;
	  strcpy(sfname, dfname);
	  rflag = 0;
#ifdef MPI
          MPI_Barrier(MPI_COMM_WORLD);
#endif
	  doit(numbp, fromp, top, rflag, maxnest);
	  rflag = 1;
	  strcpy(sfname, dfname);
	  if (CHECKFLAG && fromp == FROMPOS) {
	    fprintf(outfp, "\n");
	    for (i = 0; i < 80; ++i)
	      fprintf(outfp, "*");
	    if (htmlop)
	      fprintf(outfp, "\nSequence %4d (File: %s)\n", fcount, getDnaLink(1,sfname));
	    else
	      fprintf(outfp, "\nSequence %4d (File: %s)\n", fcount, sfname);
            fprintf(outfp, "Reverse strand.\n");
	    prt_seq_segment(outfp, gdnaR, numbp, 0, numbp - 1, rflag, NAUC, 1);
	  }
#ifdef MPI
          MPI_Barrier(MPI_COMM_WORLD);
#endif
	  doit(numbp, fromp, top, rflag, maxnest);
	  if (top == TOPOS)
	    break;
	}
	free(gdna);
	free(gdnaR);
	fflush(outfp);
      }
      closeRawTextFile(INPUT_FILE);
      if (frmterr) {
	fprintf(outfp, "\nNo GenBank-formatted input sequence data found in file %s.\n", gdnafname);
	if (oflag == 1)
	  fprintf(stdout, "\nNo GenBank-formatted input sequence data found in file %s.\n", gdnafname);

      }
    }


  if (estop && (DataSize > 0)) {
    delfile(indexName);
    delfile(dataName);
    delfile(suffixName);
    delfile(lcptreeName);
  }

  if (qpop)
    closeRawTextFile(QP_FILE);

  if (htmlop) {
#ifdef HTMLWS
    fclose(imageDataFh);
#endif
    ADD_TAIL;
  }

  if (oflag)
    fclose(outfp);

  for(i=0;i<8;i++)
    free(logV[i]);
  free(global_algnmnt);

#ifdef MPI
#ifdef MPITIMING
  if (MPIV_rank == 0) {
    time(&tlc);
    MPIV_finish = MPI_Wtime();
    printf("\nThe Master Processor 0 (of %4d) has everything wrapped up after a total of %e minutes at: %s\n",
		MPIV_nprc, (MPIV_finish - MPIV_begin) / 60, ctime(&tlc));
    fflush(stdout);
  }
#endif
  MPI_Finalize();
#endif

  return 0;

}				/* end main() */



void doit(int numbp, int ia, int ib, int rflag, int maxnest)
{
  time_t tlc;

  static struct gcalgnmnt *gca;
  static struct gcalgnmnt *gcaR;
  static struct gcalgnmnt *gcaheadp;
  static struct gpalgnmnt *gpa;
  static struct gpalgnmnt *gpaheadp;

  static int gecnt;
  static int lecnt;
  static float *pd;
  static float *pa;
  static float *pdR;
  static float *paR;

  auto thePair *hits_result;

  auto int i;
  auto int estlg;
  auto int numaa;
  static int qpcnt;
  static int phcnt;
  auto int iaR;
  auto int ibR;
  auto int gia;
  auto int gib;
  auto int iaA;
  auto int ibRA;
  auto int glgth;
  auto char *cdna;
  auto char *protein;
  auto char dfnameB[257]={0};
  auto char dfname[257]={0};
  auto char gsgmntn[513]={0};
  auto char efname[257]={0};
  auto float avgpdpa = 0.0;
  auto float MINAVGPDPA;
  auto float val = 0.0;
  auto float ExdF;
  auto int hitsCount;
  auto int hitsNum;
  auto int left;
  auto int right;
  auto int rightShift;
  auto int dbcount;

#ifdef MPI
  static int ecount;
  auto struct gcalgnmnt_buf gca_buf;
#endif

  cdna = NULL;
  protein = NULL;
  if (rflag == 0) {
    gcaheadp = NULL;
    lecnt = 0;
    gpaheadp = NULL;
  }
#ifdef DAPBM
  if (Species == 99)
    sprm.unfrm = 1;
  MINAVGPDPA = 0.5F;
#else
  if (Species == 2)
    sprm.unfrm = 1;
  MINAVGPDPA = 0.2F;
#endif
  iaR = numbp - 1 - ib;
  ibR = numbp - 1 - ia;
  iaA = ia;
  ibRA = ibR;
  if (ia > FROMPOS) {
    iaA += SGMNTSHFT;
    ibRA -= SGMNTSHFT;
  }
  else if (rflag == 0) {
    gecnt = 0;
#ifdef MPI
    ecount = 0;
#endif
  }

#ifdef MPI
  if (ecount > maxnest)	return;	/* ... maximal number of attempted EST alignments reached */
#else
  if (gecnt > maxnest)	return;	/* ... maximal number of attempted EST alignments reached */
#endif

  strcpy(dfnameB, sfname);
  strcpy(dfname, sfname);
  sprintf(gsgmntn,"%s-%d-%d",dfname,ia,ib);
  if (rflag)
    strcat(dfname, "-");
  else
    strcat(dfname, "+");

  if (rflag == 0) {
    fprintf(outfp, "\n");
    for (i = 0; i < 80; ++i)
      fprintf(outfp, "_");
    if (htmlop) {
      ADD_NAME(outfp,gsgmntn);
      fprintf(outfp, "\nSequence %4d:   %s, ", fcount, getDnaLink(1,sfname));
    }
    else
      fprintf(outfp, "\nSequence %4d:   %s, ", fcount, sfname);
    fprintf(outfp, "from %d to %d", iaA + 1, ib + 1);
    if (bflag == 1)
      fprintf(outfp, ", both strands analyzed");
    if (bflag == 2)
      fprintf(outfp, ", reverse strand");
    fprintf(outfp, ".\n");
    if (ib - iaA + 1 == 1) {
      fprintf(outfp, "\nWARNING: Empty sequence range specified.");
      fprintf(outfp, "\n  A likely cause is that the -a, -b arguments were given as greater than the sequence length (%d bp).", numbp);
      fprintf(outfp, "\n  Exiting program.\n\n");
      fprintf(stderr, "\nWARNING: Empty sequence range specified.");
      fprintf(stderr, "\n  A likely cause is that the -a, -b arguments were given as greater than the sequence length (%d bp).", numbp);
      fprintf(stderr, "\n  Exiting program.\n\n");
      exit(0);
    }
    time(&tlc);
    fprintf(outfp, "... started  at: %s\n", ctime(&tlc));
    if (htmlop && (dbop || (estop && (DataSize > 0))))
      fprintf(outfp, "\n<A HREF=\"#HEAD-PGL-%s\">Scroll down to \"Predicted gene locations\"</A>\n",gsgmntn);
  }

  if (oflag == 1) {
    if (rflag == 0) {
      fprintf(stdout, "\n");
      for (i = 0; i < 80; ++i)
	fprintf(stdout, "_");
      if (htmlop) {
        ADD_NAME(stdout,gsgmntn);
        fprintf(stdout, "\nSequence %4d:   %s, ", fcount, getDnaLink(1,sfname));
      }
      else
        fprintf(stdout, "\nSequence %4d:   %s, ", fcount, sfname);
      fprintf(stdout, "from %d to %d", iaA + 1, ib + 1);
      if (bflag == 1)
	fprintf(stdout, ", both strands analyzed");
      if (bflag == 2)
	fprintf(stdout, ", reverse strand");
      fprintf(stdout, ".\n");
    }
  }

  if (pstyle == 2 && (dbop || (estop && (DataSize > 0)))) {
    if (bflag < 2  &&  rflag == 0) {
      fprintf(outfp, "\n\nForward strand:");
      prt_seq_segment(outfp, gdna, numbp, iaA, ib, rflag, NAUC, 1);
      base_composition(gdna, iaA, ib, &tbcv);
      prt_base_composition(outfp, &tbcv);
    }
    if (bflag == 2  &&  rflag == 1) {
      fprintf(outfp, "\n\nReverse strand:");
      prt_seq_segment(outfp, gdnaR, numbp, iaR, ibRA, rflag, NAUC, 1);
      base_composition(gdnaR, iaR, ibRA, &tbcv);
      prt_base_composition(outfp, &tbcv);
    }
  }

  if (rflag == 0) {
    if ((pd = (float *) calloc((ib - ia + 2) , sizeof(float))) == NULL) {
      fatal_error("Error: memory allocation failed. Exit.\n");
    }
    if ((pa = (float *) calloc((ib - ia + 2) , sizeof(float))) == NULL) {
      fatal_error("Error: memory allocation failed. Exit.\n");
    }
#ifdef DAPBM
    det_daPbm7(&(gdna[ia]), ib - ia + 1, 0, ib - ia + 1, pd, pa, sprm.unfrm);
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
    det_daPbm7(&(gdnaR[iaR]), ib - ia + 1, 0, ib - ia + 1, pdR, paR, sprm.unfrm);
#else
    det_daPll(&(gdnaR[iaR]), ib - ia + 1, 0, ib - ia + 1, pdR, paR, sprm.unfrm);
#endif
    ComputeLogValue(pd,pa,pdR,paR,ib-ia+2);
  }

  if (dbop || (estop && (DataSize > 0))) {
#ifdef MPITIMING
      MPIV_firststart = MPI_Wtime();
#endif
    for (dbcount = 0; dbcount < estdbnum; dbcount++) {
      if (dbop) {
	MAKENAMES(estdbfnames[dbcount]);
	strcpy(estdbn, estdbfnames[dbcount]);
	ResultCreate(indexName);
      }
      if (openRawTextFile(EST_FILE, estdbn, 1) != 0) {
	fprintf(stderr, "File %s cannot be opened.\n", estdbn);
	if (htmlop) ADD_TAIL;
	exit(-1);
      }
      fprintf(stderr, "EST library file %s has been opened for reading.\n",
	      estdbn);
      fprintf(outfp, "\nEST library file:\t%s;", estdbn);
      if (rflag == 0)
	fprintf(outfp,"\tmatching gDNA +strand ...");
      else
	fprintf(outfp,"\tmatching gDNA -strand ...");
      fflush(outfp);
      if (oflag == 1) {
	fprintf(stdout, "\nEST library file:\t%s;", estdbn);
        if (rflag == 0)
	  fprintf(stdout,"\tmatching gDNA +strand ...\n\n");
        else
	  fprintf(stdout,"\tmatching gDNA -strand ...\n\n");
	fflush(stdout);
      }
      if (rflag == 0) {
	hitsNum = ResultSearch2(&(gdna[ia]), ib - ia + 1, MinMatchLen);
	rightShift = ia;
      }
      else {
	hitsNum = ResultSearch2(&(gdnaR[iaR]), ibR - iaR + 1, MinMatchLen);
	rightShift = iaR;
      }

      ResultFirst();
#ifdef MPITIMING
      MPIV_finish = MPI_Wtime();
      printf ("\nProcessor %2d used %e seconds to find matches to database %s for rflag %d",
		MPIV_rank, MPIV_finish - MPIV_start, estdbn, rflag);
      fflush(stdout);
      MPIV_start = MPI_Wtime();
#endif
      for (hitsCount = 0; hitsCount < hitsNum; hitsCount++) {

	hits_result = ResultNext();
	setPositionRawTextFile(EST_FILE, hits_result->value.OriginalPos - 1);
	if ((estlg = getlns(EST_FILE, 1, enoffset, sfname, cdna)) == -1) {
	  fprintf(stderr,"\nERROR: EST %s APPEARS TO BE OF LENGTH 0. PLEASE CHECK EST FILE.\n",sfname);
	  exit(0);
	}

	if (((hits_result->value).End_Est_Match -
	     (hits_result->value).Start_Est_Match + 1) < MinWidthOfEST ||
	    (hits_result->value.End_GenomicDNA_Match -
	     hits_result->value.Start_GenomicDNA_Match + 1) < MinWidthOfGDNA ||
	    hits_result->value.quality < MinQualityCHAIN ||
	    ((hits_result->value).End_Est_Match -
	     (hits_result->value).Start_Est_Match + 1) < MinESTcoverage*estlg) {
	  continue;
	}

#ifdef MPI
	++gecnt;
	if (++ecount > maxnest) {
#else
	if (++gecnt > maxnest) {
#endif
	  fprintf(outfp,"\nWARNING: Maximal matching EST count (%d) reached.\n", maxnest);
	  fprintf(outfp,
		  "  Only the first %d matches will be displayed for sequence %s.\n",
		  maxnest, dfnameB);
	  fprintf(outfp,
		  "  You may change this limit by setting the -m argument.\n\n");
	  if (oflag == 1) {
	    fprintf(stdout,"\nWARNING: Maximal matching EST count (%d) reached.\n",
		    maxnest);
	    fprintf(stdout,
		    "  Only the first %d matches will be displayed for sequence %s.\n",
		    maxnest, dfnameB);
	    fprintf(stdout,
		  "  You may change this limit by setting the -m argument.\n\n");
	  }
	  fprintf(stderr,"\nWARNING: Maximal matching EST count (%d) reached.\n", maxnest);
	  fprintf(stderr,
		  "  Only the first %d matches will be displayed for sequence %s.\n",
		  maxnest, dfnameB);
	  fprintf(stderr,
		  "  You may change this limit by setting the -m argument.\n\n");
	  break;	/* ... out of hitsCount loop */
	}

#ifdef MPI
#ifdef MPITIMING
	MPIV_begec = MPI_Wtime();
#endif
	if (MPIV_rank == ecount%(MPIV_nprc-1)+1) {		/* client processor */
#endif

	if ((cdna = (char *) calloc(estlg , sizeof(char))) == NULL) {
	  fatal_error("Error: memory allocation failed. Exit.\n");
	}
	if (CHECKFLAG) fprintf(stdout,"\nALLOC_CDNA1 %ld # %d",(long)cdna,++NALLOC_CDNA);
	estlg = getlns(EST_FILE, 0, enoffset, sfname, cdna);
	strcpy(efname, sfname);
	strcat(efname, "+");
	++lecnt;

#ifdef MPI
	MPIV_tag= 1;
#endif

	ExdF = (double)(MAXGLGTH - (hits_result->value.End_GenomicDNA_Match -
		hits_result->value.Start_GenomicDNA_Match) - 2*ExdLength -1) /
		(double)(estlg - (hits_result->value.End_Est_Match -
		hits_result->value.Start_Est_Match));
	if (ExdF <= 0) ExdF = 1.0;
	else           ExdF = MIN((double)ExdFactor,ExdF);
	left = (int) (hits_result->value.Start_GenomicDNA_Match - ExdF *
		hits_result->value.Start_Est_Match - ExdLength);
	right = (int) (hits_result->value.End_GenomicDNA_Match + ExdF *
		(estlg - hits_result->value.End_Est_Match) + ExdLength);
	left += rightShift;
	right += rightShift;

	if (rflag == 0) {
	  gia = (left > ia) ? left : ia;
	  gib = (right < ib) ? right : ib;
	  if ((ib < TOPOS && gib == ib) || (ia > FROMPOS && gib <= ia + SGMNTSHFT)) {
	    if (CHECKFLAG) fprintf(stdout,"\nNDEALLOC_CDNA1 %ld # %d",(long)cdna,++NDEALLOC_CDNA);
	    free(cdna);
	    --lecnt;
	    --gecnt;
#ifdef MPI
	    MPIV_tag= 0;
#ifdef MPITIMING
MPIV_bsend = MPI_Wtime();
#endif
	    MPI_Send(&gca_buf,1,MPI_gca,0,MPIV_tag,MPI_COMM_WORLD);
#ifdef MPITIMING
MPIV_asend = MPI_Wtime();
printf("\nTIME-send01 %9.6f Processor %2d alignment %5d (MPI_TAG 0, db %s, rflag %d)\n",
	MPIV_asend - MPIV_bsend, MPIV_rank, ecount, estdbn, rflag);
fflush(stdout);
#endif
#endif
	    continue;
	  }

	  glgth = gib - gia + 1;
	  if (glgth > MAXGLGTH) {
	    fprintf(outfp,"\nWARNING:\n\
A potential gene segment from %d to %d on the + strand of the gDNA\n\
(%s) exceeds the maximal gDNA length set (%d).  The matching EST is\n\
%s (length: %d).  The predicted gene may contain long introns.  The\n\
program will attempt a spliced alignment in the central %d window of the\n\
potential gene segment.  This may miss a significant alignment or display\n\
only part of a significant alignment.\n",
	gia, gib, dfname, MAXGLGTH - 3, efname, estlg, MAXGLGTH - 3);
	    gia += (glgth - MAXGLGTH + 4) / 2 + 2;
	    gib -= (glgth - MAXGLGTH + 4) / 2 - 2;
	  }

	  if (bflag == 2) {
	  if (revestallowed) {
	    efname[strlen(efname) - 1] = '-';
	    reverse_seq(cdna, estlg);
	    dfname[strlen(dfname) - 1] = '-';
	    if ((gcaR = (struct gcalgnmnt *) calloc(1 , sizeof(struct gcalgnmnt))) == NULL)
	      fatal_error("Error: memory allocation failed. Exit.\n");
	    if (CHECKFLAG) fprintf(stdout,"\nALLOC_GCA1 %ld # %d",(long)gcaR,++NALLOC_GCA);
	    strcpy(gcaR->gname, dfname);
	    strcpy(gcaR->gsgmntn, gsgmntn);
	    strcpy(gcaR->cname, efname);
	    strcpy(gcaR->estdbn, estdbn);
	    gcaR->offset = hits_result->value.OriginalPos - 1;
	    gcaR->gia = numbp - 1 - gib;
	    gcaR->gib = numbp - 1 - gia;
	    gcaR->clgth = estlg;

	    if (sahmtD(outfp, gdna, gdnaR, iaR, ibR, numbp, 1, pd, pa, pdR, paR, cdna, gcaR, sprm, gecnt)) {
	      if (insert_gca(&gcaheadp, gcaR) == 0) {
		free_gca(gcaR);
		--lecnt;
		--gecnt;
#ifdef MPI
		MPIV_tag= 0;
#ifdef MPITIMING
MPIV_bsend = MPI_Wtime();
#endif
		MPI_Send(&gca_buf,1,MPI_gca,0,MPIV_tag,MPI_COMM_WORLD);
#ifdef MPITIMING
MPIV_asend = MPI_Wtime();
printf("\nTIME-send02 %9.6f Processor %2d alignment %5d (MPI_TAG 0, db %s, rflag %d)\n",
	MPIV_asend - MPIV_bsend, MPIV_rank, ecount, estdbn, rflag);
fflush(stdout);
#endif
#endif
	      }
	      else {
		if (oflag == 1) {
		  fprintf(stdout, "%s", gcaR->algnmnt);
		  fflush(stdout);
		}
#ifdef MPI
		gca2buf(gcaR,&gca_buf);
#ifdef MPITIMING
MPIV_bsend = MPI_Wtime();
#endif
		MPI_Send(&gca_buf,1,MPI_gca,0,MPIV_tag,MPI_COMM_WORLD);
#ifdef MPITIMING
MPIV_asend = MPI_Wtime();
printf("\nTIME-send03 %9.6f Processor %2d alignment %5d (EST %d %s %4d, db %s, rflag %d)\n",
	MPIV_asend - MPIV_bsend, MPIV_rank, ecount, gca_buf.calln, gca_buf.cname, gca_buf.clgth, estdbn, rflag);
fflush(stdout);
#endif
#endif
	      }
	    }
	  } /* end revestallowed */
	  }
	  else {
	    if ((gca = (struct gcalgnmnt *) calloc(1 , sizeof(struct gcalgnmnt))) == NULL)
	      fatal_error("Error: memory allocation failed. Exit.\n");
	    if (CHECKFLAG) fprintf(stdout,"\nALLOC_GCA2 %ld # %d",(long)gca,++NALLOC_GCA);
	    strcpy(gca->gname, dfname);
	    strcpy(gca->gsgmntn, gsgmntn);
	    strcpy(gca->cname, efname);
	    strcpy(gca->estdbn, estdbn);
	    gca->offset = hits_result->value.OriginalPos - 1;
	    gca->gia = gia;
	    gca->gib = gib;
	    gca->clgth = estlg;

	    if (sahmtD(outfp, gdna, gdnaR, ia, ib, numbp, rflag, pd, pa, pdR, paR, cdna, gca, sprm, gecnt)) {
	      if (bflag == 0  ||  revestallowed == 0) {
	        if (insert_gca(&gcaheadp, gca) == 0) {
		  free_gca(gca);
		  --lecnt;
		  --gecnt;
#ifdef MPI
		  MPIV_tag= 0;
#ifdef MPITIMING
MPIV_bsend = MPI_Wtime();
#endif
		  MPI_Send(&gca_buf,1,MPI_gca,0,MPIV_tag,MPI_COMM_WORLD);
#ifdef MPITIMING
MPIV_asend = MPI_Wtime();
printf("\nTIME-send04 %9.6f Processor %2d alignment %5d (MPI_TAG 0, db %s, rflag %d)\n",
	MPIV_asend - MPIV_bsend, MPIV_rank, ecount, estdbn, rflag);
fflush(stdout);
#endif
#endif
	        }
	        else {
		  if (oflag == 1) {
		    fprintf(stdout, "%s", gca->algnmnt);
		    fflush(stdout);
		  }
#ifdef MPI
		  gca2buf(gca,&gca_buf);
#ifdef MPITIMING
MPIV_bsend = MPI_Wtime();
#endif
		  MPI_Send(&gca_buf,1,MPI_gca,0,MPIV_tag,MPI_COMM_WORLD);
#ifdef MPITIMING
MPIV_asend = MPI_Wtime();
printf("\nTIME-send05 %9.6f Processor %2d alignment %5d (EST %d %s %4d, db %s, rflag %d)\n",
	MPIV_asend - MPIV_bsend, MPIV_rank, ecount, gca_buf.calln, gca_buf.cname, gca_buf.clgth, estdbn, rflag);
fflush(stdout);
#endif
#endif
	        }
	      }
	      else {
	        if (gca->exn > 1) {
		  avgpdpa = 0.0;
		  for (i = 0; i < gca->exn - 1; ++i) {
		    avgpdpa += gca->itrscr[i][0];
		    avgpdpa += gca->itrscr[i][1];
		  }
		  avgpdpa /= (2 * (gca->exn - 1));
	        }
	        if ((gca->exn == 1 && gca->ppa[0] > gca->ppa[1])  ||
		    (gca->exn == 2 && avgpdpa < 1.5 * MINAVGPDPA) ||
		    (gca->exn  > 2 && avgpdpa <       MINAVGPDPA)   ) {
		  if (gca->exn > 1)
		    val = avgpdpa;
		  efname[strlen(efname) - 1] = '-';
		  reverse_seq(cdna, estlg);
		  dfname[strlen(dfname) - 1] = '-';
		  if ((gcaR = (struct gcalgnmnt *) calloc(1 , sizeof(struct gcalgnmnt))) == NULL)
		    fatal_error("Error: memory allocation failed. Exit.\n");
		  if (CHECKFLAG) fprintf(stdout,"\nALLOC_GCA3 %ld # %d",(long)gcaR,++NALLOC_GCA);
		  strcpy(gcaR->gname, dfname);
		  strcpy(gcaR->gsgmntn, gsgmntn);
		  strcpy(gcaR->cname, efname);
		  strcpy(gcaR->estdbn, estdbn);
		  gcaR->offset = hits_result->value.OriginalPos - 1;
		  gcaR->gia = numbp - 1 - gib;
		  gcaR->gib = numbp - 1 - gia;
		  gcaR->clgth = estlg;

		  if (sahmtD(outfp, gdna, gdnaR, iaR, ibR, numbp, 1, pd, pa, pdR, paR, cdna, gcaR, sprm, gecnt)) {
		    if (gcaR->exn > 1) {
		      avgpdpa = 0.0;
		      for (i = 0; i < gcaR->exn - 1; ++i) {
		        avgpdpa += gcaR->itrscr[i][0];
		        avgpdpa += gcaR->itrscr[i][1];
		      }
		      avgpdpa /= (2 * (gcaR->exn - 1));
		    }
		    if (gca->exn == 1 || (gcaR->exn > 1  &&  avgpdpa > val)) {
		      if (insert_gca(&gcaheadp, gcaR) == 0) {
		        free_gca(gcaR);
		        --lecnt;
		        --gecnt;
#ifdef MPI
			MPIV_tag= 0;
#ifdef MPITIMING
MPIV_bsend = MPI_Wtime();
#endif
	 		MPI_Send(&gca_buf,1,MPI_gca,0,MPIV_tag,MPI_COMM_WORLD);
#ifdef MPITIMING
MPIV_asend = MPI_Wtime();
printf("\nTIME-send06 %9.6f Processor %2d alignment %5d (MPI_TAG 0, db %s, rflag %d)\n",
	MPIV_asend - MPIV_bsend, MPIV_rank, ecount, estdbn, rflag);
fflush(stdout);
#endif
#endif
		      }
		      else {
			if (oflag == 1) {
		          fprintf(stdout, "%s", gcaR->algnmnt);
		          fflush(stdout);
			}
#ifdef MPI
			gca2buf(gcaR,&gca_buf);
#ifdef MPITIMING
MPIV_bsend = MPI_Wtime();
#endif
	 		MPI_Send(&gca_buf,1,MPI_gca,0,MPIV_tag,MPI_COMM_WORLD);
#ifdef MPITIMING
MPIV_asend = MPI_Wtime();
printf("\nTIME-send07 %9.6f Processor %2d alignment %5d (EST %d %s %4d, db %s, rflag %d)\n",
	MPIV_asend - MPIV_bsend, MPIV_rank, ecount, gca_buf.calln, gca_buf.cname, gca_buf.clgth, estdbn, rflag);
fflush(stdout);
#endif
#endif
		      }
		      free_gca(gca);
		    }
		    else {
		      if (insert_gca(&gcaheadp, gca) == 0) {
		        free_gca(gca);
		        --lecnt;
		        --gecnt;
#ifdef MPI
			MPIV_tag= 0;
#ifdef MPITIMING
MPIV_bsend = MPI_Wtime();
#endif
	 		MPI_Send(&gca_buf,1,MPI_gca,0,MPIV_tag,MPI_COMM_WORLD);
#ifdef MPITIMING
MPIV_asend = MPI_Wtime();
printf("\nTIME-send08 %9.6f Processor %2d alignment %5d (MPI_TAG 0, db %s, rflag %d)\n",
	MPIV_asend - MPIV_bsend, MPIV_rank, ecount, estdbn, rflag);
fflush(stdout);
#endif
#endif
		      }
		      else {
			if (oflag == 1) {
			  fprintf(stdout, "%s", gca->algnmnt);
			  fflush(stdout);
			}
#ifdef MPI
			gca2buf(gca,&gca_buf);
#ifdef MPITIMING
MPIV_bsend = MPI_Wtime();
#endif
	 		MPI_Send(&gca_buf,1,MPI_gca,0,MPIV_tag,MPI_COMM_WORLD);
#ifdef MPITIMING
MPIV_asend = MPI_Wtime();
printf("\nTIME-send09 %9.6f Processor %2d alignment %5d (EST %d %s %4d, db %s, rflag %d)\n",
	MPIV_asend - MPIV_bsend, MPIV_rank, ecount, gca_buf.calln, gca_buf.cname, gca_buf.clgth, estdbn, rflag);
fflush(stdout);
#endif
#endif
		      }
		      free_gca(gcaR);
		    }
		  }
		  else {
		    if (insert_gca(&gcaheadp, gca) == 0) {
		      free_gca(gca);
		      --lecnt;
		      --gecnt;
#ifdef MPI
		      MPIV_tag= 0;
#ifdef MPITIMING
MPIV_bsend = MPI_Wtime();
#endif
	 	      MPI_Send(&gca_buf,1,MPI_gca,0,MPIV_tag,MPI_COMM_WORLD);
#ifdef MPITIMING
MPIV_asend = MPI_Wtime();
printf("\nTIME-send10 %9.6f Processor %2d alignment %5d (MPI_TAG 0, db %s, rflag %d)\n",
	MPIV_asend - MPIV_bsend, MPIV_rank, ecount, estdbn, rflag);
fflush(stdout);
#endif
#endif
		    }
		    else {
		      if (oflag == 1) {
			fprintf(stdout, "%s", gca->algnmnt);
			fflush(stdout);
		      }
#ifdef MPI
		      gca2buf(gca,&gca_buf);
#ifdef MPITIMING
MPIV_bsend = MPI_Wtime();
#endif
	 	      MPI_Send(&gca_buf,1,MPI_gca,0,MPIV_tag,MPI_COMM_WORLD);
#ifdef MPITIMING
MPIV_asend = MPI_Wtime();
printf("\nTIME-send11 %9.6f Processor %2d alignment %5d (EST %d %s %4d, db %s, rflag %d)\n",
	MPIV_asend - MPIV_bsend, MPIV_rank, ecount, gca_buf.calln, gca_buf.cname, gca_buf.clgth, estdbn, rflag);
fflush(stdout);
#endif
#endif
		    }
		    free_gca(gcaR);
		  }
		  dfname[strlen(dfname) - 1] = '+';
	        }
	        else {
		  if (insert_gca(&gcaheadp, gca) == 0) {
		    free_gca(gca);
		    --lecnt;
		    --gecnt;
#ifdef MPI
		    MPIV_tag= 0;
#ifdef MPITIMING
MPIV_bsend = MPI_Wtime();
#endif
		    MPI_Send(&gca_buf,1,MPI_gca,0,MPIV_tag,MPI_COMM_WORLD);
#ifdef MPITIMING
MPIV_asend = MPI_Wtime();
printf("\nTIME-send12 %9.6f Processor %2d alignment %5d (MPI_TAG 0, db %s, rflag %d)\n",
	MPIV_asend - MPIV_bsend, MPIV_rank, ecount, estdbn, rflag);
fflush(stdout);
#endif
#endif
		  }
		  else {
		    if (oflag == 1) {
		      fprintf(stdout, "%s", gca->algnmnt);
		      fflush(stdout);
		    }
#ifdef MPI
		    gca2buf(gca,&gca_buf);
#ifdef MPITIMING
MPIV_bsend = MPI_Wtime();
#endif
	 	    MPI_Send(&gca_buf,1,MPI_gca,0,MPIV_tag,MPI_COMM_WORLD);
#ifdef MPITIMING
MPIV_asend = MPI_Wtime();
printf("\nTIME-send13 %9.6f Processor %2d alignment %5d (EST %d %s %4d, db %s, rflag %d)\n",
	MPIV_asend - MPIV_bsend, MPIV_rank, ecount, gca_buf.calln, gca_buf.cname, gca_buf.clgth, estdbn, rflag);
fflush(stdout);
#endif
#endif
		  }
	        }
	      }
	    }
	    else {
	      free_gca(gca);
	      --lecnt;
	      --gecnt;
#ifdef MPI
	      MPIV_tag= 0;
#ifdef MPITIMING
MPIV_bsend = MPI_Wtime();
#endif
	      MPI_Send(&gca_buf,1,MPI_gca,0,MPIV_tag,MPI_COMM_WORLD);
#ifdef MPITIMING
MPIV_asend = MPI_Wtime();
printf("\nTIME-send14 %9.6f Processor %2d alignment %5d (MPI_TAG 0, db %s, rflag %d)\n",
	MPIV_asend - MPIV_bsend, MPIV_rank, ecount, estdbn, rflag);
fflush(stdout);
#endif
#endif
	    }
	  }
	  fflush(outfp);
	}
	else {			/* rflag==1 */
	  gia = (left > iaR) ? left : iaR;
	  gib = (right < ibR) ? right : ibR;
	  if ((iaR > numbp - 1 - TOPOS && gia == iaR) ||
	      (ibR < numbp - 1 - FROMPOS && gia >= ibR - SGMNTSHFT)) {
	    if (CHECKFLAG) fprintf(stdout,"\nNDEALLOC_CDNA2 %ld # %d",(long)cdna,++NDEALLOC_CDNA);
	    free(cdna);
	    --lecnt;
	    --gecnt;
#ifdef MPI
	    MPIV_tag= 0;
#ifdef MPITIMING
MPIV_bsend = MPI_Wtime();
#endif
	    MPI_Send(&gca_buf,1,MPI_gca,0,MPIV_tag,MPI_COMM_WORLD);
#ifdef MPITIMING
MPIV_asend = MPI_Wtime();
printf("\nTIME-send15 %9.6f Processor %2d alignment %5d (MPI_TAG 0, db %s, rflag %d)\n",
	MPIV_asend - MPIV_bsend, MPIV_rank, ecount, estdbn, rflag);
fflush(stdout);
#endif
#endif
	    continue;
	  }

	  glgth = gib - gia + 1;
	  if (glgth > MAXGLGTH) {
	    fprintf(outfp,"\nWARNING:\n\
A potential gene segment from %d to %d on the - strand of the gDNA\n\
(%s) exceeds the maximal gDNA length set (%d).  The matching EST is\n\
%s (length: %d).  The predicted gene may contain long introns.  The\n\
program will attempt a spliced alignment in the central %d window of the\n\
potential gene segment.  This may miss a significant alignment or display\n\
only part of a significant alignment.\n",
	numbp - gib + 1, numbp - gia + 1, dfname, MAXGLGTH - 3, efname, estlg,
	MAXGLGTH - 3);
	    gia += (glgth - MAXGLGTH + 4) / 2 + 2;
	    gib -= (glgth - MAXGLGTH + 4) / 2 - 2;
	  }

	  if (bflag == 0) {
	  if (revestallowed) {
	    efname[strlen(efname) - 1] = '-';
	    reverse_seq(cdna, estlg);
	    dfname[strlen(dfname) - 1] = '+';
	    if ((gcaR = (struct gcalgnmnt *) calloc(1 , sizeof(struct gcalgnmnt))) == NULL)
	      fatal_error("Error: memory allocation failed. Exit.\n");
	    if (CHECKFLAG) fprintf(stdout,"\nALLOC_GCA4 %ld # %d",(long)gcaR,++NALLOC_GCA);
	    strcpy(gcaR->gname, dfname);
	    strcpy(gcaR->gsgmntn, gsgmntn);
	    strcpy(gcaR->cname, efname);
	    strcpy(gcaR->estdbn, estdbn);
	    gcaR->offset = hits_result->value.OriginalPos - 1;
	    gcaR->gia = numbp - 1 - gib;
	    gcaR->gib = numbp - 1 - gia;
	    gcaR->clgth = estlg;

	    if (sahmtD(outfp, gdna, gdnaR, ia, ib, numbp, 0, pd, pa, pdR, paR, cdna, gcaR, sprm, gecnt)) {
	      if (insert_gca(&gcaheadp, gcaR) == 0) {
		free_gca(gcaR);
		--lecnt;
		--gecnt;
#ifdef MPI
		MPIV_tag= 0;
#ifdef MPITIMING
MPIV_bsend = MPI_Wtime();
#endif
		MPI_Send(&gca_buf,1,MPI_gca,0,MPIV_tag,MPI_COMM_WORLD);
#ifdef MPITIMING
MPIV_asend = MPI_Wtime();
printf("\nTIME-send16 %9.6f Processor %2d alignment %5d (MPI_TAG 0, db %s, rflag %d)\n",
	MPIV_asend - MPIV_bsend, MPIV_rank, ecount, estdbn, rflag);
fflush(stdout);
#endif
#endif
	      }
	      else {
		if (oflag == 1) {
		  fprintf(stdout, "%s", gcaR->algnmnt);
		  fflush(stdout);
		}
#ifdef MPI
		gca2buf(gcaR,&gca_buf);
#ifdef MPITIMING
MPIV_bsend = MPI_Wtime();
#endif
		MPI_Send(&gca_buf,1,MPI_gca,0,MPIV_tag,MPI_COMM_WORLD);
#ifdef MPITIMING
MPIV_asend = MPI_Wtime();
printf("\nTIME-send17 %9.6f Processor %2d alignment %5d (EST %d %s %4d, db %s, rflag %d)\n",
	MPIV_asend - MPIV_bsend, MPIV_rank, ecount, gca_buf.calln, gca_buf.cname, gca_buf.clgth, estdbn, rflag);
fflush(stdout);
#endif
#endif
	      }
	    }
	  } /* end revestallowed */
	  }
	  else {
	    if ((gca = (struct gcalgnmnt *) calloc(1 , sizeof(struct gcalgnmnt))) == NULL)
	      fatal_error("Error: memory allocation failed. Exit.\n");
	    if (CHECKFLAG) fprintf(stdout,"\nALLOC_GCA5 %ld # %d",(long)gca,++NALLOC_GCA);
	    strcpy(gca->gname, dfname);
	    strcpy(gca->gsgmntn, gsgmntn);
	    strcpy(gca->cname, efname);
	    strcpy(gca->estdbn, estdbn);
	    gca->offset = hits_result->value.OriginalPos - 1;
	    gca->gia = gia;
	    gca->gib = gib;
	    gca->clgth = estlg;

	    if (sahmtD(outfp, gdna, gdnaR, iaR, ibR, numbp, rflag, pd, pa, pdR, paR, cdna, gca, sprm, gecnt)) {
	      if (bflag == 2  ||  revestallowed == 0) {
	        if (insert_gca(&gcaheadp, gca) == 0) {
		  free_gca(gca);
		  --lecnt;
		  --gecnt;
#ifdef MPI
		  MPIV_tag= 0;
#ifdef MPITIMING
MPIV_bsend = MPI_Wtime();
#endif
		  MPI_Send(&gca_buf,1,MPI_gca,0,MPIV_tag,MPI_COMM_WORLD);
#ifdef MPITIMING
MPIV_asend = MPI_Wtime();
printf("\nTIME-send18 %9.6f Processor %2d alignment %5d (MPI_TAG 0, db %s, rflag %d)\n",
	MPIV_asend - MPIV_bsend, MPIV_rank, ecount, estdbn, rflag);
fflush(stdout);
#endif
#endif
	        }
	        else {
		  if (oflag == 1) {
		    fprintf(stdout, "%s", gca->algnmnt);
		    fflush(stdout);
		  }
#ifdef MPI
		  gca2buf(gca,&gca_buf);
#ifdef MPITIMING
MPIV_bsend = MPI_Wtime();
#endif
		  MPI_Send(&gca_buf,1,MPI_gca,0,MPIV_tag,MPI_COMM_WORLD);
#ifdef MPITIMING
MPIV_asend = MPI_Wtime();
printf("\nTIME-send19 %9.6f Processor %2d alignment %5d (EST %d %s %4d, db %s, rflag %d)\n",
	MPIV_asend - MPIV_bsend, MPIV_rank, ecount, gca_buf.calln, gca_buf.cname, gca_buf.clgth, estdbn, rflag);
fflush(stdout);
#endif
#endif
	        }
	      }
	      else {
	        if (gca->exn > 1) {
		  avgpdpa = 0.0;
		  for (i = 0; i < gca->exn - 1; ++i) {
		    avgpdpa += gca->itrscr[i][0];
		    avgpdpa += gca->itrscr[i][1];
		  }
		  avgpdpa /= (2 * (gca->exn - 1));
	        }
	        if ((gca->exn == 1 && gca->ppa[0] > gca->ppa[1]) ||
		    (gca->exn == 2 && avgpdpa < 1.5 * MINAVGPDPA) ||
		    (gca->exn  > 2 && avgpdpa <       MINAVGPDPA)   ) {
		  if (gca->exn > 1)
		    val = avgpdpa;
		  efname[strlen(efname) - 1] = '-';
		  reverse_seq(cdna, estlg);
		  dfname[strlen(dfname) - 1] = '+';
		  if ((gcaR = (struct gcalgnmnt *) calloc(1 , sizeof(struct gcalgnmnt))) == NULL)
		    fatal_error("Error: memory allocation failed. Exit.\n");
		  if (CHECKFLAG) fprintf(stdout,"\nALLOC_GCA6 %ld # %d",(long)gcaR,++NALLOC_GCA);
		  strcpy(gcaR->gname, dfname);
		  strcpy(gcaR->gsgmntn, gsgmntn);
		  strcpy(gcaR->cname, efname);
		  strcpy(gcaR->estdbn, estdbn);
		  gcaR->offset = hits_result->value.OriginalPos - 1;
		  gcaR->gia = numbp - 1 - gib;
		  gcaR->gib = numbp - 1 - gia;
		  gcaR->clgth = estlg;

		  if (sahmtD(outfp, gdna, gdnaR, ia, ib, numbp, 0, pd, pa, pdR, paR, cdna, gcaR, sprm, gecnt)) {
		    if (gcaR->exn > 1) {
		      avgpdpa = 0.0;
		      for (i = 0; i < gcaR->exn - 1; ++i) {
		        avgpdpa += gcaR->itrscr[i][0];
		        avgpdpa += gcaR->itrscr[i][1];
		      }
		      avgpdpa /= (2 * (gcaR->exn - 1));
		    }
		    if (gca->exn == 1 || (gcaR->exn > 1  &&  avgpdpa > val)) {
		      if (insert_gca(&gcaheadp, gcaR) == 0) {
		        free_gca(gcaR);
		        --lecnt;
		        --gecnt;
#ifdef MPI
			MPIV_tag= 0;
#ifdef MPITIMING
MPIV_bsend = MPI_Wtime();
#endif
			MPI_Send(&gca_buf,1,MPI_gca,0,MPIV_tag,MPI_COMM_WORLD);
#ifdef MPITIMING
MPIV_asend = MPI_Wtime();
printf("\nTIME-send20 %9.6f Processor %2d alignment %5d (MPI_TAG 0, db %s, rflag %d)\n",
	MPIV_asend - MPIV_bsend, MPIV_rank, ecount, estdbn, rflag);
fflush(stdout);
#endif
#endif
		      }
		      else {
			if (oflag == 1) {
			  fprintf(stdout, "%s", gcaR->algnmnt);
			  fflush(stdout);
			}
#ifdef MPI
			gca2buf(gcaR,&gca_buf);
#ifdef MPITIMING
MPIV_bsend = MPI_Wtime();
#endif
			MPI_Send(&gca_buf,1,MPI_gca,0,MPIV_tag,MPI_COMM_WORLD);
#ifdef MPITIMING
MPIV_asend = MPI_Wtime();
printf("\nTIME-send21 %9.6f Processor %2d alignment %5d (EST %d %s %4d, db %s, rflag %d)\n",
	MPIV_asend - MPIV_bsend, MPIV_rank, ecount, gca_buf.calln, gca_buf.cname, gca_buf.clgth, estdbn, rflag);
fflush(stdout);
#endif
#endif
		      }
		      free_gca(gca);
		    }
		    else {
		      if (insert_gca(&gcaheadp, gca) == 0) {
		        free_gca(gca);
		        --lecnt;
		        --gecnt;
#ifdef MPI
			MPIV_tag= 0;
#ifdef MPITIMING
MPIV_bsend = MPI_Wtime();
#endif
			MPI_Send(&gca_buf,1,MPI_gca,0,MPIV_tag,MPI_COMM_WORLD);
#ifdef MPITIMING
MPIV_asend = MPI_Wtime();
printf("\nTIME-send22 %9.6f Processor %2d alignment %5d (MPI_TAG 0, db %s, rflag %d)\n",
	MPIV_asend - MPIV_bsend, MPIV_rank, ecount, estdbn, rflag);
fflush(stdout);
#endif
#endif
		      }
		      else {
			if (oflag == 1) {
			  fprintf(stdout, "%s", gca->algnmnt);
			  fflush(stdout);
			}
#ifdef MPI
			gca2buf(gca,&gca_buf);
#ifdef MPITIMING
MPIV_bsend = MPI_Wtime();
#endif
			MPI_Send(&gca_buf,1,MPI_gca,0,MPIV_tag,MPI_COMM_WORLD);
#ifdef MPITIMING
MPIV_asend = MPI_Wtime();
printf("\nTIME-send23 %9.6f Processor %2d alignment %5d (EST %d %s %4d, db %s, rflag %d)\n",
	MPIV_asend - MPIV_bsend, MPIV_rank, ecount, gca_buf.calln, gca_buf.cname, gca_buf.clgth, estdbn, rflag);
fflush(stdout);
#endif
#endif
		      }
		      free_gca(gcaR);
		    }
		  }
		  else {
		    if (insert_gca(&gcaheadp, gca) == 0) {
		      free_gca(gca);
		      --lecnt;
		      --gecnt;
#ifdef MPI
		      MPIV_tag= 0;
#ifdef MPITIMING
MPIV_bsend = MPI_Wtime();
#endif
		      MPI_Send(&gca_buf,1,MPI_gca,0,MPIV_tag,MPI_COMM_WORLD);
#ifdef MPITIMING
MPIV_asend = MPI_Wtime();
printf("\nTIME-send24 %9.6f Processor %2d alignment %5d (MPI_TAG 0, db %s, rflag %d)\n",
	MPIV_asend - MPIV_bsend, MPIV_rank, ecount, estdbn, rflag);
fflush(stdout);
#endif
#endif
		    }
		    else {
		      if (oflag == 1) {
			fprintf(stdout, "%s", gca->algnmnt);
			fflush(stdout);
		      }
#ifdef MPI
		      gca2buf(gca,&gca_buf);
#ifdef MPITIMING
MPIV_bsend = MPI_Wtime();
#endif
		      MPI_Send(&gca_buf,1,MPI_gca,0,MPIV_tag,MPI_COMM_WORLD);
#ifdef MPITIMING
MPIV_asend = MPI_Wtime();
printf("\nTIME-send25 %9.6f Processor %2d alignment %5d (EST %d %s %4d, db %s, rflag %d)\n",
	MPIV_asend - MPIV_bsend, MPIV_rank, ecount, gca_buf.calln, gca_buf.cname, gca_buf.clgth, estdbn, rflag);
fflush(stdout);
#endif
#endif
		    }
		    free_gca(gcaR);
		  }
		  dfname[strlen(dfname) - 1] = '-';
	        }
	        else {
		  if (insert_gca(&gcaheadp, gca) == 0) {
		    free_gca(gca);
		    --lecnt;
		    --gecnt;
#ifdef MPI
		    MPIV_tag= 0;
#ifdef MPITIMING
MPIV_bsend = MPI_Wtime();
#endif
		    MPI_Send(&gca_buf,1,MPI_gca,0,MPIV_tag,MPI_COMM_WORLD);
#ifdef MPITIMING
MPIV_asend = MPI_Wtime();
printf("\nTIME-send26 %9.6f Processor %2d alignment %5d (MPI_TAG 0, db %s, rflag %d)\n",
	MPIV_asend - MPIV_bsend, MPIV_rank, ecount, estdbn, rflag);
fflush(stdout);
#endif
#endif
		  }
		  else {
		    if (oflag == 1) {
		      fprintf(stdout, "%s", gca->algnmnt);
		      fflush(stdout);
		    }
#ifdef MPI
		    gca2buf(gca,&gca_buf);
#ifdef MPITIMING
MPIV_bsend = MPI_Wtime();
#endif
		    MPI_Send(&gca_buf,1,MPI_gca,0,MPIV_tag,MPI_COMM_WORLD);
#ifdef MPITIMING
MPIV_asend = MPI_Wtime();
printf("\nTIME-send27 %9.6f Processor %2d alignment %5d (EST %d %s %4d, db %s, rflag %d)\n",
	MPIV_asend - MPIV_bsend, MPIV_rank, ecount, gca_buf.calln, gca_buf.cname, gca_buf.clgth, estdbn, rflag);
fflush(stdout);
#endif
#endif
		  }
	        }
	      }
	    }
	    else {
	      free_gca(gca);
	      --lecnt;
	      --gecnt;
#ifdef MPI
	      MPIV_tag= 0;
#ifdef MPITIMING
MPIV_bsend = MPI_Wtime();
#endif
	      MPI_Send(&gca_buf,1,MPI_gca,0,MPIV_tag,MPI_COMM_WORLD);
#ifdef MPITIMING
MPIV_asend = MPI_Wtime();
printf("\nTIME-send28 %9.6f Processor %2d alignment %5d (MPI_TAG 0, db %s, rflag %d)\n",
	MPIV_asend - MPIV_bsend, MPIV_rank, ecount, estdbn, rflag);
fflush(stdout);
#endif
#endif
	    }
	  }
	  fflush(outfp);
	}
	if (CHECKFLAG) fprintf(stdout,"\nNDEALLOC_CDNA3 %ld # %d",(long)cdna,++NDEALLOC_CDNA);
	free(cdna);

#ifdef MPI
#ifdef MPITIMING
MPIV_finish = MPI_Wtime();
if (MPIV_tag == 1) {
  printf("\nTIME-2send  %9.6f Processor %2d alignment %5d (EST %d %s %4d, db %s, rflag %d)\n",
	MPIV_finish - MPIV_begec, MPIV_rank, ecount, gca_buf.calln, gca_buf.cname, gca_buf.clgth, estdbn, rflag);
} else {
  printf("\nTIME-2send  %9.6f Processor %2d alignment %5d (MPI_TAG 0, db %s, rflag %d)\n",
	MPIV_finish - MPIV_begec, MPIV_rank, ecount, estdbn, rflag);
}
fflush(stdout);
#endif
	}
	else if (MPIV_rank == 0) {				/* master processor */
	  MPIV_source= ecount%(MPIV_nprc-1)+1;
	  MPI_Recv(&gca_buf,1,MPI_gca,MPIV_source,MPI_ANY_TAG,MPI_COMM_WORLD,&MPIV_status);
#ifdef MPITIMING
MPIV_finish = MPI_Wtime();
if ((MPIV_status).MPI_TAG == 1) {
  printf("\nIt took %e seconds until the Master Processor received alignment %5d (EST %d %s %4d, db %s, rflag %d) from Processor %2d\n",
	MPIV_finish - MPIV_begin, ecount, gca_buf.calln, gca_buf.cname, gca_buf.clgth, estdbn, rflag, MPIV_source);
  printf("\nTIME-2recv  %9.6f Processor %2d alignment %5d (EST %d %s %4d, db %s, rflag %d)\n",
	MPIV_finish - MPIV_begec, MPIV_source, ecount, gca_buf.calln, gca_buf.cname, gca_buf.clgth, estdbn, rflag);
} else {
  printf("\nIt took %e seconds until the Master Processor received alignment %5d (MPI_TAG 0, db %s, rflag %d) from Processor %2d\n",
	MPIV_finish - MPIV_begin, ecount, estdbn, rflag, MPIV_source);
  printf("\nTIME-2recv  %9.6f Processor %2d alignment %5d (MPI_TAG 0, db %s, rflag %d)\n",
	MPIV_finish - MPIV_begec, MPIV_rank, ecount, estdbn, rflag);
}
fflush(stdout);
#endif
	  if ((MPIV_status).MPI_TAG == 1) {
	    if ( (gca= (struct gcalgnmnt *) malloc(sizeof(struct gcalgnmnt))) == NULL)
	      fatal_error("Error: memory allocation failed. Exit.\n");
	    buf2gca(&gca_buf,gca);
	    if (insert_gca(&gcaheadp,gca) != 0) ++lecnt;
	  }
	}
#endif

      }				/* end of for(hitsCount ... */

      ResultDestroy();
      closeRawTextFile(EST_FILE);
#ifdef MPITIMING
      MPIV_finish = MPI_Wtime();
      printf ("\nProcessor %2d used %e seconds to generate the alignments to database %s for rflag %d\n",
		MPIV_rank, MPIV_finish - MPIV_start, estdbn, rflag);
      fflush(stdout);
      MPIV_start = MPI_Wtime();
#endif
    }				/* end of for(dbcount ... */

#ifdef MPI
    if (MPIV_rank==0) {
    if (rflag == 1  ||  ecount > maxnest) {
      if (ecount == 0  &&  ecount <= maxnest) {
#else
    if (rflag == 1  ||  gecnt > maxnest) {
      if (lecnt == 0  &&  gecnt <= maxnest) {
#endif
	if (htmlop)
	  fprintf(outfp, "<A NAME=\"HEAD-PGL-%s\"></A>\n",gsgmntn);
	fprintf(outfp, "\n\nNo significant EST matches were found.\n");
	if (oflag == 1)
	  fprintf(stdout, "\n\nNo significant EST matches were found.\n");
      }
      else {
#ifndef HTMLWS
	if (htmlop && oflag == 1) {
	  fprintf(stdout, "<A NAME=\"HEAD-PGL-%s\"></A>\n",gsgmntn);
	  fprintf(stdout, "\nPlease see file %s for predicted gene locations based on\nconsensus EST evidence.\n",outfname);
	}
#endif
#ifdef MPITIMING
	time(&tlc);
	MPIV_finish = MPI_Wtime();
	printf("\nAfter %e minutes, the Master Processor 0 enters det_pgl at: %s",
		(MPIV_finish - MPIV_begin) / 60, ctime(&tlc));
	fflush(stdout);
	MPIV_start = MPI_Wtime();
#endif
	det_pgl(outfp, &gcaheadp, dfname, gsgmntn, gdna, gdnaR, ia, ib, numbp, pd, pa, pdR, paR, sprm, bflag, revestallowed);
#ifdef MPITIMING
	time(&tlc);
	MPIV_finish = MPI_Wtime();
	printf("\nAfter %e minutes, the Master Processor 0 exits  det_pgl at: %s\n",
		(MPIV_finish - MPIV_start) / 60, ctime(&tlc));
	fflush(stdout);
#endif
      }
      if (htmlop) {
        fprintf(outfp, "<A NAME=\"BOTTOM-PGL-%s\"></A>\n",gsgmntn);
	if (oflag == 1)
          fprintf(stdout, "<A NAME=\"BOTTOM-PGL-%s\"></A>\n",gsgmntn);
      }
      fflush(outfp);
      free_gca_list(gcaheadp);
    }
#ifdef MPI
    }
    else {
      if (rflag == 1) {
	fflush(outfp);
	free_gca_list(gcaheadp);
      }
#ifdef MPITIMING
      time(&tlc);
      MPIV_finish = MPI_Wtime();
      printf("\nProcessor %2d used %e seconds for sequence %s (rflag %d)",
		MPIV_rank, MPIV_finish - MPIV_firststart, gsgmntn, rflag);
      if (rflag == 1)
	printf(" and finished after a total of %e minutes at: %s", (MPIV_finish - MPIV_begin) / 60, ctime(&tlc));
      fflush(stdout);
#endif
    }
#endif

  }				/* end if (dbop || estop) */


  if (qpop) {
    if (bflag == 1  ||  (bflag == 0  &&  rflag == 0)  ||  (bflag == 2  &&  rflag == 1)) {
      setPositionRawTextFile(QP_FILE, 0);
      qpcnt = 0;
      if (bflag != 1  ||  rflag == 0) phcnt = 0;
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
	if (protein[numaa-1] != 23) {
          protein[numaa] = 23;
          ++numaa;
        }

	if (rflag == 0) {
	  gia = ia;
	  gib = ib;
	  if ((gpa = (struct gpalgnmnt *) calloc(1 , sizeof(struct gpalgnmnt))) == NULL)
	    fatal_error("Error: memory allocation failed. Exit.\n");
	  strcpy(gpa->gname, dfname);
	  strcpy(gpa->gsgmntn, gsgmntn);
	  strcpy(gpa->pname, sfname);
	  strcpy(gpa->qpdbn, qpfname);
	  gpa->gia = gia;
	  gpa->gib = gib;
	  gpa->plgth = numaa;
	  logPia=ia;
	  if (sahmtP(outfp, &(gdna[gia]), dfname, gib - gia + 1, gia,
                     numbp, rflag, &(pd[gia - ia]), &(pa[gia - ia]),
                     protein, sfname, numaa, smat, gpa, sprm, qpcnt)) {
	    if (oflag == 1)
	      prt_gpa_list(stdout, gpa, gsgmntn);
	    insert_gpa(&gpaheadp, gpa);
	    ++phcnt;
	  }
	  else
	    free_gpa(gpa);
	}
	else {
	  gia = iaR;
	  gib = ibR;
	  if ((gpa = (struct gpalgnmnt *) calloc(1 , sizeof(struct gpalgnmnt))) == NULL)
	    fatal_error("Error: memory allocation failed. Exit.\n");
	  strcpy(gpa->gname, dfname);
	  strcpy(gpa->gsgmntn, gsgmntn);
	  strcpy(gpa->pname, sfname);
	  strcpy(gpa->qpdbn, qpfname);
	  gpa->gia = gia;
	  gpa->gib = gib;
	  gpa->plgth = numaa;
	  logPia=iaR;
	  if (sahmtP(outfp, &(gdnaR[gia]), dfname, gib - gia + 1, gia,
                     numbp, rflag, &(pdR[gia - iaR]), &(paR[gia - iaR]),
                     protein, sfname, numaa, smat, gpa, sprm, qpcnt)) {
	    if (oflag == 1)
	      prt_gpa_list(stdout, gpa, gsgmntn);
	    insert_gpa(&gpaheadp, gpa);
	    ++phcnt;
	  }
	  else
	    free_gpa(gpa);
	}
        free(protein);
      }
    }

    if (rflag == 1) {
      if (phcnt == 0) {
	fprintf(outfp, "\n\nNo significant protein matches were found.\n");
	if (oflag == 1)
	  fprintf(stdout, "\n\nNo significant protein matches were found.\n");
      }
      else
	prt_gpa_list(outfp, gpaheadp, gsgmntn);
      fflush(outfp);
      free_gpa_list(gpaheadp);
    }
  }				/* end if (qpop .. */

  if (rflag == 1) {
    free(pd);
    free(pa);
    free(pdR);
    free(paR);
    time(&tlc);
    fprintf(outfp, "... finished at: %s\n", ctime(&tlc));
  }

}				/* end doit() */



void fatal_error(char *buf)
{

  fprintf(stderr, "%s\n", buf);
  exit(1);

}				/* end fatal_error() */



#ifdef MPI
void buf2gca(struct gcalgnmnt_buf *buf, struct gcalgnmnt *gca)
{
  int len;

  memcpy(gca,buf,sizeof(struct gcalgnmnt)-3*sizeof(int)); /* leave three pointers out */
  len=strlen(buf->algnmnt);
  if ( (gca->algnmnt= (char *) malloc((len+1)*sizeof(char))) == NULL)
    fatal_error("Error: memory allocation failed. Exit.\n");
  strcpy(gca->algnmnt,buf->algnmnt);
  gca->link=NULL;
  gca->next=NULL;		
}



void gca2buf(struct gcalgnmnt *gca, struct gcalgnmnt_buf *buf)
{

  memcpy(buf,gca,sizeof(struct gcalgnmnt)-3*sizeof(int)); /* leave three pointers out */
  if (gca->algnmnt != NULL) {
    strcpy(buf->algnmnt,gca->algnmnt);
  }
  buf->link=NULL;
  buf->next=NULL;		
}



void MPI_derived_type_gca(struct gcalgnmnt_buf *gca,MPI_Datatype *MPI_gca)
{
  int bl[23];			/* block lengths: number of elements in each block */
  MPI_Aint dp[23];		/* displacemtns of each element from start of new type */
  MPI_Datatype tl[23];		/* data type */
  MPI_Aint sa, a;		/* start address and address */
   
  /* elm 0: char gname[257] */
  bl[0]= 257;	
  tl[0]= MPI_CHAR;
  dp[0]= 0;				
  MPI_Get_address(gca->gname,&sa);	
   
  /* elm 1: char gsgmntn[257] */ 
  bl[1]= 257;
  tl[1]= MPI_CHAR;
  MPI_Get_address(gca->gsgmntn,&a);
  dp[1]= a-sa;
   
  /* elm2: char cname[257] */
  bl[2]= 257;
  tl[2]= MPI_CHAR;
  MPI_Get_address(gca->cname,&a);
  dp[2]= a-sa;
   
  /* elm3: char estdbn[257] */
  bl[3]= 257;
  tl[3]= MPI_CHAR;
  MPI_Get_address(gca->estdbn,&a);
  dp[3]= a-sa;
   
  /* elm4: int calln */
  bl[4]= 1;
  tl[4]= MPI_INT;
  MPI_Get_address(&(gca->calln),&a);
  dp[4]= a-sa;
   
  /* elm5: int offset */
  bl[5]= 1;
  tl[5]= MPI_INT;
  MPI_Get_address(&(gca->offset),&a);
  dp[5]= a-sa;
   
  /* elm6: int gia */
  bl[6]= 1;
  tl[6]= MPI_INT;
  MPI_Get_address(&(gca->gia),&a);
  dp[6]= a-sa;
   
  /* elm7: int gib */
  bl[7]= 1;
  tl[7]= MPI_INT;
  MPI_Get_address(&(gca->gib),&a);
  dp[7]= a-sa;
   
  /* elm8: int clgth */
  bl[8]= 1;
  tl[8]= MPI_INT;
  MPI_Get_address(&(gca->clgth),&a);
  dp[8]= a-sa;
   
  /* elm9: int exn */
  bl[9]= 1;
  tl[9]= MPI_INT;
  MPI_Get_address(&(gca->exn),&a);
  dp[9]= a-sa;
   
  /* elm10: int gcds[MAXNEXNS][2] */
  bl[10]= MAXNEXNS*2;
  tl[10]= MPI_INT;
  MPI_Get_address(&(gca->gcds[0][0]),&a);
  dp[10]= a-sa;
   
  /* elm11: int ccds[MAXNEXNS][2] */
  bl[11]= MAXNEXNS*2;
  tl[11]= MPI_INT;
  MPI_Get_address(&(gca->ccds[0][0]),&a);
  dp[11]= a-sa;   
   
  /* elm12: int ppa[2] */
  bl[12]= 2;
  tl[12]= MPI_INT;
  MPI_Get_address(&(gca->ppa),&a);
  dp[12]= a-sa;
   
  /* elm13: float exnscr[MAXNEXNS] */
  bl[13]= MAXNEXNS;
  tl[13]= MPI_FLOAT;
  MPI_Get_address(&(gca->exnscr),&a);
  dp[13]= a-sa;
   
  /* elm14: float itscr[MAXNEXNS][2] */
  bl[14]= MAXNEXNS*2;
  tl[14]= MPI_FLOAT;
  MPI_Get_address(&(gca->itrscr),&a);
  dp[14]= a-sa;
   
  /* elm15: float score */
  bl[15]= 1;
  tl[15]= MPI_FLOAT;
  MPI_Get_address(&(gca->score),&a);
  dp[15]= a-sa;
   
  /* elm16: int pcds[MAXNPPS][MAXNEXNS][2];*/
  bl[16]= 2*MAXNPPS*MAXNEXNS;
  tl[16]= MPI_INT;
  MPI_Get_address(&(gca->pcds),&a);
  dp[16]= a-sa;
   
  /* elm17: int npsg[MAXNPPS] */
  bl[17]= MAXNPPS	;
  tl[17]= MPI_INT;
  MPI_Get_address(&(gca->npsg),&a);
  dp[17]= a-sa;
   
  /*elm18: int npps; */
  bl[18]= 1;
  tl[18]= MPI_INT;
  MPI_Get_address(&(gca->npps),&a);
  dp[18]= a-sa;
   
  /*elm19: int ags[MAXNAGS] */
  bl[19]= MAXNAGS;
  tl[19]= MPI_INT;
  MPI_Get_address(&(gca->ags),&a);
  dp[19]= a-sa;
   
  /* elm20: char algnmnt[6*MAXGLGTH+2*MAXCLGTH+2*MAXNEXNS*80]; */
  bl[20]= 6*MAXGLGTH+2*MAXCLGTH+2*MAXNEXNS*80;
  tl[20]= MPI_CHAR;
  MPI_Get_address(gca->algnmnt,&a);
  dp[20]= a-sa;
   
  /* elm21: struct gcalgnmnt *link; */
  bl[21]= 1;
  tl[21]= MPI_INT;
  MPI_Get_address(&(gca->link),&a);
  dp[21]= a-sa;
   
  /* elm22: struct gcalgnmnt *next; */
  bl[22]= 1;
  tl[22]= MPI_INT;
  MPI_Get_address(&(gca->next),&a);
  dp[22]= a-sa;
   
  MPI_Type_create_struct(23,bl,dp,tl,MPI_gca);
  MPI_Type_commit(MPI_gca);
   
} /* end MPI_derived_type_gca() */
#endif
