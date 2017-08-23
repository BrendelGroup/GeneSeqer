/*sequence.h                                         Version of June 2, 2000. */

#ifndef _sequence_h_
#define _sequence_h_

int gm_pcnt[3][4], gm_pcntS[3], gm_ccnt[4], gm_ccntS = 0;
int gd_pcnt[3][4][4], gd_pcntS[3], gd_ccnt[4][4], gd_ccntS = 0;
int gt_pcnt[3][4][4][4], gt_pcntS[3], gt_ccnt[4][4][4], gt_ccntS = 0;
int gr_pcnt[3][4][4][4][4], gr_pcntS[3], gr_ccnt[4][4][4][4], gr_ccntS = 0;
int gp_pcnt[3][4][4][4][4][4], gp_pcntS[3], gp_ccnt[4][4][4][4][4], gp_ccntS =
0;
int gh_pcnt[3][4][4][4][4][4][4], gh_pcntS[3], gh_ccnt[4][4][4][4][4][4],
  gh_ccntS = 0;

void reverse_seq (char *seq, int numbp);
void complement_seq (char *seq, int numbp, char *seqR);
void reverse_ftc (int ftc[MAXNBAF][2], int nbaf, int numbp);
int base_composition (char *seq, int ia, int ib, struct bcvct *bcv);
void prt_base_composition (FILE * fp, struct bcvct *bcv);
void word_composition (FILE * fp, char *seq, int crd[MAXNBAF][2], int ia,
		       int ib, int pflag, int cflag, int tflag,
		       float MINgcPCNT, float MAXgcPCNT);
void prt_word_composition (FILE * fp, int m_cnt[4], int nm, int d_cnt[4][4],
			   int nd, int t_cnt[4][4][4], int nt,
			   int r_cnt[4][4][4][4], int nr,
			   int p_cnt[4][4][4][4][4], int np,
			   int h_cnt[4][4][4][4][4][4], int nh);
void prt_seq_segment (FILE * fp, char *seq, int numbp, int ia, int ib,
		      int rflag, char ABC[], int nflag);
void prt_seq_segment_to_str (char *str, char *seq, int numbp, int ia, int ib,
			     int rflag, char ABC[], int nflag);
void prt_segment_trl (FILE * fp, int seq[], int numbp, int ia, int ib,
		      int rflag, int phsA, int cl, char ABC[]);
int trl_seq_segment (int dseq[], int ia, int ib, int phsA, int pseq[],
		     int imrkr[], int nintrns);
int translate (int dseq[], int ia, int ib, int pseq[], int numaa);

#endif
