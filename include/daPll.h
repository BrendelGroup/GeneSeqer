/*daPll.h;                 Last update: February 15, 2014, by Volker Brendel. */
/*Logitlinear splice site model parameters and subroutines for GeneSeqer.c    */

#define NSPECIES 2
#define NMDLS  3
#define NLEVELS  4
#define DONSIZE  7
#define ACTSIZE 15


/* CURRENT SETTINGS:

   species 0 = Zea mays
   1 = Arabidopsis thaliana

   model 0 = logitlinear, without subclassification
   1 = logitlinear, main subclass (GGT for donors, CAG for acceptors)
   2 = logitlinear, minor subclass (HGT for donors, DAG for acceptors)

   levels 0 = display all GU (AG) sites provided at least 50 bases to the
   left and right are available in the input sequence
   1 = threshold corresponding to 100\% sensitivity on the training set
   2 = threshold corresponding to  95\% sensitivity on the training set
   3 = threshold corresponding to maximal Tau value on the training set

   notation: ct = constant term (alpha)
   fU = weight for U-contrast (delta)
   fS = weight for S-contrast (mu)
   wt = weights for signal positions (-5 to +4 for donors, -15 to
   +4 for acceptors; note: the GU and AG consensus sites are
   assumed present and are not scored)
 */

char species_name[NSPECIES + 1][25] =
{"Maize", "Arabidopsis thaliana", "generic"};
char model_name[NMDLS][50] =
{"no subclassification",
 "subclassification G/H-GU, C/D-AG"
};

double dtv[NSPECIES][NMDLS][NLEVELS] =
{
  {
    {-0.0001, 0.0013, 0.0440, 0.2823},
    {-0.0001, 0.0189, 0.0701, 0.5411},
    {-0.0001, 0.0050, 0.0081, 0.2189}
  },
  {
    {-0.0001, 0.0017, 0.0426, 0.2515},
    {-0.0001, 0.0036, 0.0750, 0.2332},
    {-0.0001, 0.0028, 0.0200, 0.2940}
  }
};



double dct[NSPECIES][NMDLS] =
{
  {-12.063100, 3.3567, 0.1692},
  {-11.236000, 3.4594, 0.3677}
};
double dfU[NSPECIES][NMDLS] =
{
  {-10.217900, -8.2873, -17.6069},
  {-9.297700, -9.7055, -8.8654}
};
double dfS[NSPECIES][NMDLS] =
{
  {9.283300, 10.1065, 6.5313},
  {12.869700, 13.6899, 13.0875}
};

double dwt[NSPECIES][NMDLS][DONSIZE][4] =
{
  {
    {
      {0.000000, 0.959300, 1.116500, 0.312500},
      {0.000000, -0.009900, 1.926400, -0.793700},
      {0.000000, -0.193500, -0.046200, 3.698300},
      {0.000000, 0.591300, 4.037900, 1.538800},
      {0.000000, 0.919300, 2.176100, -0.108000},
      {0.000000, -0.322900, 0.198100, 2.741100},
      {0.000000, -0.785700, -1.468900, -1.852900}
    },
    {
      {-1.636300, -0.046800, 0.000000, -1.171200},
      {-2.157300, -2.533000, 0.000000, -2.911400},
      {0.000000, 0.000000, 0.000000, 0.000000},
      {-3.759100, -3.376700, 0.000000, -2.390200},
      {-1.504200, -0.658800, 0.000000, -1.648300},
      {-2.499500, -2.915700, -2.428400, 0.000000},
      {0.000000, -0.824100, -1.255200, -1.782600}
    },
    {
      {0.638500, -0.388900, 0.000000, 0.364500},
      {-0.906500, -0.921800, 0.000000, -2.723400},
      {0.000000, -0.239800, -0.446100, 0.000000},
      {-11.613000, -3.932200, 0.000000, -2.667000},
      {-4.893800, -3.915800, 0.000000, -11.429200},
      {-4.035800, -4.424800, -4.858900, 0.000000},
      {0.000000, -1.201400, -2.648300, -2.569600}
    }
  },
  {
    {
      {0.000000, 1.171000, 0.873900, -0.069300},
      {0.000000, 0.214100, 2.382700, -0.092900},
      {0.000000, -0.140000, -0.302000, 3.476700},
      {0.000000, 0.071600, 3.170500, 0.224000},
      {0.000000, 0.707200, 2.094300, -0.478500},
      {0.000000, -0.312200, 0.031900, 3.005500},
      {0.000000, -0.822800, -1.311300, -1.748100}
    },
    {
      {-0.796200, 0.319200, 0.000000, -1.199300},
      {-2.742000, -3.211500, 0.000000, -3.046300},
      {0.000000, 0.000000, 0.000000, 0.000000},
      {-2.867100, -3.059100, 0.000000, -2.988700},
      {-1.777200, -1.082400, 0.000000, -2.108900},
      {-2.720100, -2.824400, -2.529600, 0.000000},
      {0.000000, -0.820200, -1.353000, -1.797000}
    },
    {
      {-0.906200, 0.203000, 0.000000, -0.793100},
      {-1.258600, -0.707300, 0.000000, -1.291400},
      {0.000000, -0.105400, -0.172700, 0.000000},
      {-5.207600, -3.176000, 0.000000, -2.674100},
      {-2.710100, -2.130600, 0.000000, -4.523800},
      {-3.727800, -5.561500, -3.979700, 0.000000},
      {0.000000, -0.896500, -1.166000, -1.539000}
    }
  }
};



double atv[NSPECIES][NMDLS][NLEVELS] =
{
  {
    {-0.0001, 0.0018, 0.0358, 0.1701},
    {-0.0001, 0.0018, 0.0488, 0.2697},
    {-0.0001, 0.0004, 0.0501, 0.3577}
  },
  {
    {-0.0001, 0.0010, 0.0536, 0.2773},
    {-0.0001, 0.0013, 0.0973, 0.3577},
    {-0.0001, 0.0006, 0.0456, 0.3012}
  }
};

double act[NSPECIES][NMDLS] =
{
  {-3.107600, 3.5991, -2.6044},
  {-1.259600, 5.1969, 2.7254}
};
double afU[NSPECIES][NMDLS] =
{
  {15.006500, 12.8622, 28.6355},
  {10.345600, 11.3194, 10.7442}
};
double afS[NSPECIES][NMDLS] =
{
  {-9.000700, -10.8430, -11.4700},
  {-12.507500, -11.4307, -13.8548}
};

double awt[NSPECIES][NMDLS][ACTSIZE][4] =
{
  {
    {
      {0.000000, -0.401200, -0.787100, -1.148900},
      {0.000000, -0.093500, -0.553300, -0.499800},
      {0.000000, -0.093600, -1.171000, -0.073100},
      {0.000000, -0.128900, -0.860700, 0.036500},
      {0.000000, -0.320600, -0.770000, -0.603200},
      {0.000000, -0.128900, -0.437600, -0.630000},
      {0.000000, 0.301400, -0.650800, -0.499300},
      {0.000000, -0.230900, -0.663000, -0.547900},
      {0.000000, -0.891400, -1.017900, -1.110300},
      {0.000000, -0.814300, -0.818600, -0.804200},
      {0.000000, -0.625300, -1.376100, -1.176300},
      {0.000000, -0.072400, 0.627800, 1.394200},
      {0.000000, 1.807700, -1.634300, -10.279900},
      {0.000000, 0.206600, 1.408500, 2.955200},
      {0.000000, -1.147100, -1.832900, -1.409700}
    },
    {
      {0.000000, -0.477800, -0.586100, -1.183800},
      {0.000000, -0.086500, -0.473500, -0.284300},
      {0.000000, -0.140300, -1.445100, -0.333200},
      {0.000000, 0.031300, -1.114000, -0.152200},
      {0.000000, -0.340100, -1.228700, -0.421300},
      {0.000000, -0.084700, -0.378300, -0.468900},
      {0.000000, 0.256200, -0.486000, -0.477900},
      {0.000000, 0.006300, -0.740700, -0.430200},
      {0.000000, -0.836800, -1.152000, -1.351800},
      {0.000000, -0.947600, -0.746200, -1.066900},
      {0.000000, -0.807800, -1.484400, -1.363800},
      {-2.098400, -1.887300, -0.726700, 0.000000},
      {0.000000, 0.000000, 0.000000, 0.000000},
      {-2.792500, -2.870700, -1.257100, 0.000000},
      {0.000000, -1.541300, -1.757200, -1.563600}
    },
    {
      {0.000000, -0.376800, -1.744100, -1.496700},
      {0.000000, -0.281100, -1.303700, -1.654400},
      {0.000000, 0.345100, -1.633600, 1.678100},
      {0.000000, 0.311600, 0.111900, 1.450800},
      {0.000000, -0.371200, -0.452400, -1.963300},
      {0.000000, 0.697100, -0.349100, -0.969600},
      {0.000000, 1.906800, -1.030900, -0.074600},
      {0.000000, -0.652300, 0.095700, -0.910400},
      {0.000000, -1.458300, -0.688800, -0.620000},
      {0.000000, -0.018700, -0.948900, -0.352900},
      {0.000000, -0.423800, -2.678400, -1.138400},
      {0.460000, -0.891100, -0.546200, 0.000000},
      {0.000000, 0.000000, -1.929500, -13.017400},
      {-5.885500, -3.197000, -3.157500, 0.000000},
      {0.000000, 0.195600, -2.610400, -2.266400}
    }
  },
  {
    {
      {0.000000, -0.480500, -0.494500, -0.369600},
      {0.000000, -0.210300, -0.769700, -0.725600},
      {0.000000, -0.945900, -0.994300, -0.755400},
      {0.000000, -1.065100, -0.890600, -0.996700},
      {0.000000, -0.934600, -0.686100, -0.731600},
      {0.000000, -0.643400, -0.596500, -0.851200},
      {0.000000, -0.917700, -0.871400, -0.787700},
      {0.000000, -1.086800, -0.812900, -0.760000},
      {0.000000, -1.198000, -0.707000, -0.966400},
      {0.000000, -1.046500, -0.878100, -0.700500},
      {0.000000, -1.221400, -0.734300, -1.250100},
      {0.000000, -0.178700, 1.057100, 1.534000},
      {0.000000, 1.819300, -2.468800, -5.535600},
      {0.000000, -0.320200, 0.764300, 2.604200},
      {0.000000, -1.013100, -1.306300, -1.033900}
    },
    {
      {0.000000, -0.428200, -0.251000, -0.445000},
      {0.000000, -0.375000, -0.725100, -0.748300},
      {0.000000, -0.830200, -0.874900, -1.013400},
      {0.000000, -1.334700, -1.019000, -0.999200},
      {0.000000, -1.532800, -0.781400, -0.783800},
      {0.000000, -0.970100, -0.645200, -1.316200},
      {0.000000, -0.940200, -1.240800, -0.993800},
      {0.000000, -1.320900, -0.989100, -1.105400},
      {0.000000, -1.426500, -0.361000, -1.001600},
      {0.000000, -1.386300, -1.062900, -0.983900},
      {0.000000, -1.278700, -0.683100, -1.199300},
      {-1.488600, -1.908800, -0.330400, 0.000000},
      {0.000000, 0.000000, 0.000000, 0.000000},
      {-2.383400, -2.873300, -1.680400, 0.000000},
      {0.000000, -1.046000, -1.090900, -1.167200}
    },
    {
      {0.000000, -0.643800, -1.086400, -0.401800},
      {0.000000, -0.022600, -0.945700, -0.789800},
      {0.000000, -1.244000, -1.315100, -0.484200},
      {0.000000, -0.785400, -1.075400, -1.074100},
      {0.000000, -0.117100, -0.440000, -0.599500},
      {0.000000, -0.258300, -0.589900, -0.251400},
      {0.000000, -1.010700, -0.412500, -0.424900},
      {0.000000, -0.924200, -0.568500, -0.432400},
      {0.000000, -0.835600, -1.411900, -0.843700},
      {0.000000, -0.929100, -0.709500, -0.428200},
      {0.000000, -1.344200, -0.829600, -1.586800},
      {-1.581700, -1.695000, -0.629300, 0.000000},
      {0.000000, 0.000000, -2.610000, -5.900600},
      {-3.355200, -3.179000, -2.433800, 0.000000},
      {0.000000, -1.133200, -2.069500, -1.000000}
    }
  }
};

double dqntl[2][3][4] =
{
  {
    {0.760, 0.520, 0.240, 0.075},
    {0.174, 0.054, 0.010, 0.002},
    {0.600, 0.300, 0.085, 0.000001}
  },
  {
    {0.720, 0.449, 0.180, 0.073},
    {0.116, 0.024, 0.005, 0.001},
    {0.497, 0.224, 0.060, 0.009}}
};
double aqntl[2][3][4] =
{
  {
    {0.710, 0.500, 0.250, 0.110},
    {0.148, 0.063, 0.015, 0.003},
    {0.570, 0.192, 0.071, 0.000001}
  },
  {
    {0.718, 0.553, 0.240, 0.085},
    {0.186, 0.072, 0.015, 0.002},
    {0.458, 0.170, 0.049, 0.006}}
};

int Species = -1, Wlgth = 50, Model = 1, Level = 1;

struct ssite_stats
  {
    int Y, loc, bq;
    double pval, ct, Lscr, delU, cU, delS, cS, logit;
  }
dsite, asite;



int qlty_branchpoint (char seq[], int numbp, int pos, int wa, int wb)
{
  int i;

  {
    if (pos - wa - 1 >= 0 && pos - wb - 2 < numbp)
      {
	for (i = pos - wa + 1; i <= pos - wb - 2; ++i)
	  if (seq[i] == 1 && seq[i + 1] == 0 && seq[i + 3] == 2)
	    return (2);
	/* CTNA pattern */
	for (i = pos - wa + 1; i <= pos - wb - 2; ++i)
	  if (seq[i] == 0 && seq[i + 1] == 0 && seq[i + 3] == 2)
	    return (1);
	/* TTNA pattern */
      }
    return (0);
  }

}				/* end qlty_branchpoint() */



double is_donor (char seq[], int numbp, int pos, int wl, int species,
	         int model, struct ssite_stats *psite)
{
  int i, j;
  struct bcvct bcvl, bcvr;
  double delU, delS, Lscr = 0.0, logit, pval = -0.5;
  if (seq[pos] != 3 || seq[pos + 1] != 0)
    return (-1.0);
  else
    {
      if (pos - wl >= 0 && pos + wl + 1 < numbp)
	{
	  base_composition (seq, pos - wl, pos - 1, &bcvl);
	  base_composition (seq, pos + 2, pos + wl + 1, &bcvr);
	  delU = bcvl.bfrq[0] - bcvr.bfrq[0];
	  delS =
	    (bcvl.bfrq[1] + bcvl.bfrq[3]) - (bcvr.bfrq[1] + bcvr.bfrq[3]);
	  for (i = pos - 3; i <= pos - 1; ++i)
	    {
	      if (seq[i] != 10)
		Lscr += dwt[species][model][i - pos + 3][(int)seq[i]];
	      else
		for (j = 0; j < 4; ++j)
		  Lscr += (0.25 * dwt[species][model][i - pos + 3][j]);
	    }
	  for (i = pos + 2; i <= pos + 5; ++i)
	    {
	      if (seq[i] != 10)
		Lscr += dwt[species][model][i - pos + 1][(int)seq[i]];
	      else
		for (j = 0; j < 4; ++j)
		  Lscr += (0.25 * dwt[species][model][i - pos + 1][j]);
	    }
	  logit = dct[species][model] + dfU[species][model] * delU +
	    dfS[species][model] * delS + Lscr;
	  pval = exp (logit) / (exp (logit) + 1.00);
	  pval = exp (logit) / (exp (logit) + 1.00);
	  psite->ct = dct[species][model];
	  psite->delU = delU;
	  psite->cU = dfU[species][model] * delU;
	  psite->delS = delS;
	  psite->cS = dfS[species][model] * delS;
	  psite->Lscr = Lscr;
	  psite->logit = logit;
	  psite->pval = pval;
	  psite->bq = -1;
	}
      return (pval);
    }

}				/* end is_donor() */



double is_acptr (char seq[], int numbp, int pos, int wl, int species, int model,
	         struct ssite_stats *psite)
{
  int i, j;
  int qlty_branchpoint (char seq[], int numbp, int pos, int wa, int wb);
  struct bcvct bcvl, bcvr;
  double delU, delS, Lscr = 0.0;
  double logit, pval=0.0;

  if (seq[pos - 1] != 2 || seq[pos] != 3)
    return (-1.0);
  else
    {
      if (pos - wl - 1 >= 0 && pos + wl < numbp)
	{
	  base_composition (seq, pos - wl - 1, pos - 2, &bcvl);
	  base_composition (seq, pos + 1, pos + wl, &bcvr);
	  delU = bcvl.bfrq[0] - bcvr.bfrq[0];
	  delS =
	    (bcvl.bfrq[1] + bcvl.bfrq[3]) - (bcvr.bfrq[1] + bcvr.bfrq[3]);
	  for (i = pos - 14; i <= pos - 2; ++i)
	    {
	      if (seq[i] != 10)
		Lscr += awt[species][model][i - pos + 14][(int)seq[i]];
	      else
		for (j = 0; j < 4; ++j)
		  Lscr += (0.25 * awt[species][model][i - pos + 14][j]);
	    }
	  for (i = pos + 1; i <= pos + 2; ++i)
	    {
	      if (seq[i] != 10)
		Lscr += awt[species][model][i - pos + 12][(int)seq[i]];
	      else
		for (j = 0; j < 4; ++j)
		  Lscr += (0.25 * awt[species][model][i - pos + 12][j]);
	    }
	  logit = act[species][model] + afU[species][model] * delU +
	    afS[species][model] * delS + Lscr;

	  pval = exp (logit) / (exp (logit) + 1.00);
	  pval = exp (logit) / (exp (logit) + 1.00);
	  psite->ct = act[species][model];
	  psite->delU = delU;
	  psite->cU = afU[species][model] * delU;
	  psite->delS = delS;
	  psite->cS = afS[species][model] * delS;
	  psite->Lscr = Lscr;
	  psite->logit = logit;
	  psite->pval = pval;
	  psite->bq = qlty_branchpoint (seq, numbp, pos, 60, 21);
	}
      return (pval);
    }

}				/* end is_acptr() */



double donorP(char *seq, int numbp, int pos, int wl, int species, int model)
{
  int i, j;
  struct bcvct bcvl, bcvr;
  double delU, delS, Lscr = 0.0, logit, pval = -0.5;
  if (seq[pos] != 3 || seq[pos + 1] != 0)
    return (-1.0);
  else {
    if (pos - wl >= 0 && pos + wl + 1 < numbp) {
      base_composition(seq, pos - wl, pos - 1, &bcvl);
      base_composition(seq, pos + 2, pos + wl + 1, &bcvr);
      delU = bcvl.bfrq[0] - bcvr.bfrq[0];
      delS = (bcvl.bfrq[1] + bcvl.bfrq[3]) - (bcvr.bfrq[1] + bcvr.bfrq[3]);
      for (i = pos - 3; i <= pos - 1; ++i) {
	if (seq[i] != 10)
	  Lscr += dwt[species][model][i - pos + 3][(int)seq[i]];
	else
	  for (j = 0; j < 4; ++j)
	    Lscr += (0.25 * dwt[species][model][i - pos + 3][j]);
      }
      for (i = pos + 2; i <= pos + 5; ++i) {
	if (seq[i] != 10)
	  Lscr += dwt[species][model][i - pos + 1][(int)seq[i]];
	else
	  for (j = 0; j < 4; ++j)
	    Lscr += (0.25 * dwt[species][model][i - pos + 1][j]);
      }
      logit = dct[species][model] + dfU[species][model] * delU +
	dfS[species][model] * delS + Lscr;
      pval = exp(logit) / (exp(logit) + 1.00);
      pval = exp(logit) / (exp(logit) + 1.00);
    }
    return (pval);
  }

}				/* end donorP() */



double acptrP(char *seq, int numbp, int pos, int wl, int species, int model)
{
  int i, j;
  struct bcvct bcvl, bcvr;
  double delU, delS, Lscr = 0.0;
  double logit, pval;
  pval = 0.0;
  if (seq[pos - 1] != 2 || seq[pos] != 3)
    return (-1.0);
  else {
    if (pos - wl - 1 >= 0 && pos + wl < numbp) {
      base_composition(seq, pos - wl - 1, pos - 2, &bcvl);
      base_composition(seq, pos + 1, pos + wl, &bcvr);
      delU = bcvl.bfrq[0] - bcvr.bfrq[0];
      delS = (bcvl.bfrq[1] + bcvl.bfrq[3]) - (bcvr.bfrq[1] + bcvr.bfrq[3]);
      for (i = pos - 14; i <= pos - 2; ++i) {
	if (seq[i] != 10)
	  Lscr += awt[species][model][i - pos + 14][(int)seq[i]];
	else
	  for (j = 0; j < 4; ++j)
	    Lscr += (0.25 * awt[species][model][i - pos + 14][j]);
      }
      for (i = pos + 1; i <= pos + 2; ++i) {
	if (seq[i] != 10)
	  Lscr += awt[species][model][i - pos + 12][(int)seq[i]];
	else
	  for (j = 0; j < 4; ++j)
	    Lscr += (0.25 * awt[species][model][i - pos + 12][j]);
      }
      logit = act[species][model] + afU[species][model] * delU +
	afS[species][model] * delS + Lscr;

      pval = exp(logit) / (exp(logit) + 1.00);
      pval = exp(logit) / (exp(logit) + 1.00);
    }
    return (pval);
  }

}				/* end acptrP() */



void det_daPll (char *seq, int numbp, int ia, int ib, float *pd, float *pa, int unifrm)
{
  int i, smodel;
  double donorP(char *seq, int numbp, int pos, int wl, int species, int model);
  double dscr, ascr;
  double acptrP(char *seq, int numbp, int pos, int wl, int species, int model);
  double dmin, amin;

  smodel = 0;
  pd[0] = pa[0] = 0.000001F;
  for (i = 1; i < numbp - 1; ++i) {
    if (seq[i] == 3) {
      if (seq[i + 1] == 0) {
	if (unifrm)
	  pd[i] = 0.05F;
	else
	  pd[i] = 0.001F;
      }
      else if (seq[i + 1] == 1)
	pd[i] = 0.0001F;
      else
	pd[i] = 0.000001F;
      if (seq[i - 1] == 2) {
	if (unifrm)
	  pa[i] = 0.05F;
	else
	  pa[i] = 0.001F;
      }
      else
	pa[i] = 0.000001F;
    }
    else {
      pd[i] = pa[i] = 0.000001F;
    }
  }
  pd[numbp - 1] = pa[numbp - 1] = pd[numbp] = pa[numbp] = 0.000001F;
  if (unifrm)
    return;


  for (i = ia; i <= ib; ++i) {
    if (i - Wlgth >= 0 && i + Wlgth + 1 < numbp) {
      if (Model == 0)
	smodel = 0;
      else if (Model == 1) {
	if (seq[i - 1] == 3)
	  smodel = 1;
	else
	  smodel = 2;
      }
      dmin = dtv[Species][smodel][Level];
      if ((dscr = donorP(seq, numbp, i, Wlgth, Species, smodel)) >= dmin)
	pd[i] = (float) (dscr);
    }
    if (i - Wlgth - 1 >= 0 && i + Wlgth < numbp) {
      if (Model == 0)
	smodel = 0;
      else if (Model == 1) {
	if (seq[i - 2] == 1)
	  smodel = 1;
	else
	  smodel = 2;
      }
      amin = atv[Species][smodel][Level];
      if ((ascr = acptrP(seq, numbp, i, Wlgth, Species, smodel)) >= amin)
	pa[i] = (float) (ascr);
    }

  }

}				/* end det_daPll() */
