#ifndef _HTML_H_
#define _HTML_H_
#include <stdio.h>


static char link[1024]="";
extern FILE *imageDataFh;

#ifndef HTMLWSF
#ifdef HTMLWS
static int hwsflag = 1;
#else
static int hwsflag = 0;
#endif 
#define HTMLWSF
#endif

#define ADD_HEAD {fprintf(outfp,"<!DOCTYPE HTML PUBLIC \"-//IETF//DTD HTML//EN\">\n<HTML><HEAD><TITLE>GeneSeqer Output</TITLE></HEAD><BODY bgcolor=white text=black link=blue><FONT FACE=\"Courier\"><PRE>\n"); if(!hwsflag && oflag) fprintf(stdout,"<!DOCTYPE HTML PUBLIC \"-//IETF//DTD HTML//EN\">\n<HTML><HEAD><TITLE>GeneSeqer Output</TITLE></HEAD><BODY bgcolor=white text=black link=blue><FONT FACE=\"Courier\"><PRE>\n");}

#define ADD_TAIL {fprintf(outfp,"</PRE></BODY></HTML>\n"); if(!hwsflag && oflag) fprintf(stdout,"</PRE></BODY></HTML>\n");}

#define ADD_NAME(fp,name) {fprintf(fp,"\n<A NAME=\"%s\"></A>\n",name);}

#define ADD_PGL_LINK(fp,name) {fprintf(fp, "\n  <A HREF=\"#%s\">PGS</A> (",name);}


#define BLASTP_HEAD "  http://www.ncbi.nlm.nih.gov/blast/Blast.cgi?CMD=Web&LAYOUT=TwoWindows&AUTO_FORMAT=Semiauto&ALIGNMENTS=50&ALIGNMENT_VIEW=Pairwise&CDD_SEARCH=on&CLIENT=web&COMPOSITION_BASED_STATISTICS=on&DATABASE=nr&DESCRIPTIONS=100&ENTREZ_QUERY=%28none%29&EXPECT=10&FILTER=L&FORMAT_OBJECT=Alignment&FORMAT_TYPE=HTML&I_THRESH=0.005&MATRIX_NAME=BLOSUM62&NCBI_GI=on&PAGE=Proteins&PROGRAM=blastp&SERVICE=plain&SET_DEFAULTS.x=41&SET_DEFAULTS.y=5&SHOW_OVERVIEW=on&END_OF_HTTPGET=Yes&SHOW_LINKOUT=yes&GET_SEQUENCE=yes&QUERY="

#define BLASTP_TAIL "&END_OF_HTTPGET=Yes"


#ifdef FPROTOTYPE
static char *getName(long address)
#else
static char *getName(address) 
    long address;
#endif
{
    extern char link[];
    sprintf(link,"\n<A NAME=\"%ld\"></A>\n",address);
    return link;
}
    

#ifdef FPROTOTYPE
static char *getDnaLink(int db, char *gid)
#else
static char *getDnaLink(db, gid) 
    int db;
    char *gid;
#endif
{
    extern char link[];
    char gidb[257];
    
    strcpy(gidb,gid);
    if (gidb[strlen(gidb) - 1] == '+'  ||  gidb[strlen(gidb) - 1] == '-') {
      gidb[strlen(gidb) - 1] = '\0';
    }
    if (db == 1) {	/* GenBank */
      sprintf(link,"<A HREF=\"http://www.ncbi.nlm.nih.gov/entrez/query.fcgi?cmd=Search&db=Nucleotide&term=%s&doptcmdl=GenBank\" TARGET=\"NUCLEOTIDE SEQUENCE SEARCH\">%s</A>",gidb,gid);
    } else if (db == 2) {	/* PlantGDB */
      sprintf(link,"<A HREF=\"http://www.plantgdb.org/search/display/data.php?Seq_ID=%s\" TARGET=\"NUCLEOTIDE SEQUENCE SEARCH\">%s</A>",gidb,gid);
    }
    else {	/* default: GenBank */
      sprintf(link,"<A HREF=\"http://www.ncbi.nlm.nih.gov/entrez/query.fcgi?cmd=Search&db=Nucleotide&term=%s&doptcmdl=GenBank\" TARGET=\"NUCLEOTIDE SEQUENCE SEARCH\">%s</A>",gidb,gid);
    }
    return link;
} 


#ifdef FPROTOTYPE
static char *getProteinLink(char *pid)
#else
static char *getProteinLink(pid) 
    char *pid;
#endif
{
    extern char link[];
    sprintf(link,"<A HREF=\"http://www.ncbi.nlm.nih.gov/entrez/query.fcgi?cmd=Search&db=Protein&term=%s&doptcmdl=GenPept\" TARGET=\"PROTEIN SEQUENCE SEARCH\">%s</A>",pid,pid);
    return link;
} 

#endif
