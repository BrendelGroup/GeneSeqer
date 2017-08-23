#!/usr/bin/env python
'''Allows loading gthxml PGS data into xGDB (http://sourceforge.net/projects/xgdb)

This module allows one to load GeneSeqer or GenomeThreader spliced
alignment results, for nucleotide- and/or protein-type query sequences,
into an xGDB backend.  It parses results stored in the gthxml format, as
defined at http://www.public.iastate.edu/~mespar1/gthxml/GenomeThreader.rng.txt .
(Note that GenomeThreader can generate this output directly using the
-xmlout option, whereas GeneSeqer/GeneSeqer2 plain text output can be
retrofitted to this standard via the GSQ2XML.pl tool, available
at http://www.public.iastate.edu/~mespar1/gthxml/GSQ2XML.pl .)

Our local xGDB administrators have historically loaded only the
transcript-based hqPGS results from GeneSeqer output; I refer to these data
as ``autocorrelated spliced alignments".  There is (currently) no such equivalent
in the GenomeThreader world, though.  Therefore, the user should specify
if his/her input data is based on GeneSeqer results .AND. if
they wish to use hqPGS data rather than raw spliced alignments
(set --geneseqer); the default for this is GenomeThreader-style.
Note that hqPGS data are not available for protein spliced alignments.

Usage:
    python xGDBload_PgsFromgthxml.py
        -h --help           # prints this message and exits
        -x --xreffile file  # cross reference table file for genomic template indices
                            # per-line format is either 
                            #   '(\S+)\s+(\S+)\s+chr(\d+)' .OR.
                            #   '(\S+)\s+(\S+)\s+(\d+)'
                            # depending on whether a chromosome was used as the genomic
                            # reference type or not, respectively. The third token will
                            # be used as a unique numeric key for a given genomic template,
                            # while the identifier for that template, as used in the
                            # spliced alignment gthxml output, must be indicated in
                            # either field one or two of the record.
        -c --cxreffile file # cross reference file for unique identifiers, e.g., GI
                            # numbers, of full-length cDNA sequences. Optional. These are
                            # not required to be numeric, but cannot contain whitespace.
        -p --pxreffile file # cross reference file for unique identifiers, e.g., GI
                            # numbers, of protein sequences.  These are not required to
                            # be numeric, but cannot contain whitespace.  If the gthxml
                            # data contains protein spliced alignment data, this option
                            # is mandatory.
        -i --infile file    # Input file identifier. Input stream can be gthxml data
                            # piped by stdin (default), a file on the resident file
                            # system, or the URL for a (strictly!) gzip'ed gthxml file
                            # from the Web
        -o --outfile file   # Output file identifier. Can be stdout (default) or a
                            # local file, which will be overwritten if it exists
        -t --table name     # Handle for xGDB tables to be populated (default is 'good_pgs')
        -u --cdnatable name # If query sequences contain full-length cDNAs, you can
                            # specify a table name for storing cDNA-based results, separate
                            # from those of ESTs and/or proteins. Use in conjunction with
                            # the -c option, if desired (default table name is
                            # 'cdna_good_pgs', enacted only if -c is user-defined).
        -v --prottable name # If query sequences contain proteins, you can specify a
                            # table name for storing protein-based results, separate
                            # from those of ESTs/cDNAs. Use in conjunction with
                            # the -p option, if desired (default table name is
                            # 'prot_good_pgs', enacted only if -p is user-defined).
        --chromosome        # Were alignments made with fully assembled pseudomolecules?
        --geneseqer         # Load hqPGS data (autocorrelated spliced alignments) from
                            # GeneSeqer results, rather than PGS?  Note that hqPGS data
                            # are not available for protein alignments.
        --nuclloscore float # minimum overall alignment score of EST/cDNA alignments (default 0.8)
        --nuclhiscore float # maximum overall alignment score of EST/cDNA alignments (default 2.0)
        --protloscore float # minimum overall alignment score of protein alignments (default 0.3)
        --prothiscore float # maximum overall alignment score of protein alignments (default 2.0)
        --nucllocov float   # minimum coverage of max(template,EST/cDNA) coverage (default 0.8)
        --nuclhicov float   # maximum coverage of max(template,EST/cDNA) coverage (default 2.0)
        --protlocov float   # minimum coverage of max(template,protein) coverage (default 0.5)
        --prothicov float   # maximum coverage of max(template,protein) coverage (default 2.0)
'''

__author__ = 'Michael E Sparks (mespar1@iastate.edu)'
__version__ = '1.2.1'
__copyright__ = 'Copyright (c) 2006 Michael E Sparks, all rights reserved'
__license__ = 'GNU General Public License (GPL)'

# Version history
# 1.2.1 - Records 'PGS_line' data regardless, overwriting
#         with hqPGS data, if desired
#         Modified some comments
#
# 1.2 - Added support for parsing protein alignment data
#       cDNA xref file now only requires a single column
#       Modified some comments
#
# 1.1 - Minor modifications to comments
#       Enclosed gi column in single quotes
#
# 1.0 - Initial release
#       Parsed EST/cDNA-based PGS and hqPGS results

import getopt, re, sys
from xml.sax.handler import ContentHandler

ABS_MIN_SCORE = 0.00
ABS_MAX_SCORE = 2.00 # I don't like this, but it works...mes
ABS_MIN_COV = 0.00
ABS_MAX_COV = 2.00 # I don't like this, but it works...mes

class xGDBgthxmlPGSError(Exception):
    '''Exception subclass for addressing fatal errors in xGDBgthxmlPGS objects'''
    pass

class xGDBgthxmlPGS(ContentHandler):
    '''A handler to populate xGDB with spliced alignment results in gthxml format.

    Allows one to parse a gthxml v1.0 file from a local filesystem or over the
    Web and print executable SQL statements to populate appropriate xGDB backend
    tables.
    '''

    def __init__(self,argv):
        '''Class "constructor" method'''

        # Set default values prior to parsing argument list
        # Note that these may trigger exceptions later in the class
        # methods, but this is the intended effect!  In short, don't
        # monkey with this part of the code.
        self.filexref = '/dev/null'   # cross ref table file for genomic template indices
        self.cfilexref = '/dev/null'  # cross ref table file for cDNA identification
                                      # strings, such as GI numbers (optional)
        self.cfilexrefp = False       # Did user specify a cDNA cross ref file?
        self.pfilexref = '/dev/null'  # cross ref table file for prot identification
                                      # strings, such as GI numbers (optional)
        self.pfilexrefp = False       # Did user specify a protein cross ref file?
        self.table = 'good_pgs'       # handle for the xGDB tables we want to populate
        self.ctable = 'cdna_good_pgs' # handle for xGDB tables we want to populate with cDNA info
        self.ptable = 'prot_good_pgs' # handle for xGDB tables we want to populate with prot info
        self.chrp = False             # Were alignments made with fully assembled pseudomolecules?
        self.gsqp = False             # Input is hqPGS from GeneSeqer, or no?  (Not possible if
                                      # source program was GenomeThreader, nor for protein
                                      # alignments!)
        self.infile = sys.stdin       # Input file identifier. Input stream can be gthxml data
                                      # piped by stdin, a file on the resident file system, or a
                                      # (strictly!) gzip'ed gthxml file from the Web
        self.instreamp = False        # Have we opened the input stream yet?
        self.webp = False             # If self.infile contains an 'http://' prefix, this will be
                                      # be set to True in the code below.  User will not directly
                                      # modify this data member
        self.outfile = sys.stdout     # Output file identifier. Can be stdout or a local file,
                                      # which will be overwritten if it exists
        self.outstreamp = False       # Have we opened the output stream yet?
        self.nuclloscore = 0.80          # minimum overall alignment score of EST/cDNA alignments
        self.nuclhiscore = ABS_MAX_SCORE # maximum overall alignment score of EST/cDNA alignments
        self.protloscore = 0.30          # minimum overall alignment score of protein alignments
        self.prothiscore = ABS_MAX_SCORE # maximum overall alignment score of protein alignments
        self.nucllocov = 0.80         # minimum coverage of max(template,EST/cDNA) coverage
        self.nuclhicov = ABS_MAX_COV  # maximum coverage of max(template,EST/cDNA) coverage
        self.protlocov = 0.50         # minimum coverage of max(template,protein) coverage
        self.prothicov = ABS_MAX_COV  # maximum coverage of max(template,protein) coverage
        self.scoreDNE = '-1.0'        # Stub value for protein donor/acceptor scores,
                                      # which do not exist.
        self.phase_symbols = ('a', 'b', 'c') # Used to correlate gap phases to symbols
                                             # in protein alignments

        # Parse option list
        try:
            opts, args = getopt.getopt(argv, 'hx:c:p:t:u:v:i:o:',
                                            ['help',
                                             'xreffile=',
                                             'cxreffile=',
                                             'pxreffile=',
                                             'table=',
                                             'cdnatable=',
                                             'prottable=',
                                             'infile=',
                                             'outfile=',
                                             'chromosome',
                                             'geneseqer',
                                             'nuclloscore=', 'nuclhiscore=',
                                             'protloscore=', 'prothiscore=',
                                             'nucllocov=', 'nuclhicov=',
                                             'protlocov=', 'prothicov='])
        except getopt.GetoptError:
            print __doc__
            sys.exit(2)

        for o, a in opts:
            if o in ('-h', '--help'):
                print __doc__
                sys.exit(None)
            elif o in ('-x', '--xreffile'):
                self.filexref = a
            elif o in ('-c', '--cxreffile'):
                self.cfilexref = a
            elif o in ('-p', '--pxreffile'):
                self.pfilexref = a
            elif o in ('-t', '--table'):
                self.table = a
            elif o in ('-u', '--cdnatable'):
                self.ctable = a
            elif o in ('-v', '--prottable'):
                self.ptable = a
            elif o in ('-i', '--infile'):
                self.infile = a
            elif o in ('-o', '--outfile'):
                self.outfile = a
            elif o == '--chromosome':
                self.chrp = True
            elif o == '--geneseqer':
                self.gsqp = True
            elif o == '--nuclloscore':
                self.nuclloscore = float(a)
            elif o == '--nuclhiscore':
                self.nuclhiscore = float(a)
            elif o == '--protloscore':
                self.protloscore = float(a)
            elif o == '--prothiscore':
                self.prothiscore = float(a)
            elif o == '--nucllocov':
                self.nucllocov = float(a)
            elif o == '--nuclhicov':
                self.nuclhicov = float(a)
            elif o == '--protlocov':
                self.protlocov = float(a)
            elif o == '--prothicov':
                self.prothicov = float(a)

        # Verify nucleotide/protein score and coverage values are reasonable.
        self.__validateparms()

        # Populate cross-reference dictionaries
        self.__xrefLookup(self.filexref,self.chrp)
        self.__cpxrefLookup(self.cfilexref,'cDNA')
        self.__cpxrefLookup(self.pfilexref,'prot')

        # Attempt to open input stream
        if self.infile is not sys.stdin:
            if re.search(r'^http://',self.infile):
                # This flag will come in handy when we need to close compressed_gthxml_stream
                self.webp = True
 
                # I only support parsing gzip'ed gthxml over the web,
                # to enforce lower overhead on network traffic.
                import gzip, StringIO, urllib

                websocket = urllib.urlopen(self.infile)
                compressed_gthxml = websocket.read()
                websocket.close()

                self.compressed_gthxml_stream = StringIO.StringIO(compressed_gthxml)
                self.instream = gzip.GzipFile(fileobj=self.compressed_gthxml_stream)
            else:
                self.instream = open(self.infile,"r")
        else:
            self.instream = self.infile
        self.instreamp = True

        # Attempt to open output stream
        if self.outfile is not sys.stdout:
            self.outstream = open(self.outfile,"w")
        else:
            self.outstream = self.outfile
        self.outstreamp = True

        # initialize various data members
        self.__reset()

    def __validateparms(self):
        '''Verify nucleotide/protein score and coverage values are reasonable.'''

        # Verify nucleotide score values were legitimate
        if self.nuclloscore > self.nuclhiscore:
            sys.stderr.write('Warning: nuclloscore > nuclhiscore, swapping instead!??\n')
            garbage = self.nuclloscore
            self.nuclloscore = self.nuclhiscore
            self.nuclhiscore = garbage
        if self.nuclloscore < ABS_MIN_SCORE:
            sys.stderr.write('Warning: Bad --nuclloscore, using %s instead\n' \
                             % str(ABS_MIN_SCORE))
            self.nuclloscore = ABS_MIN_SCORE
        if self.nuclhiscore > ABS_MAX_SCORE:
            sys.stderr.write('Warning: Bad --nuclhiscore, using %s instead\n' \
                             % str(ABS_MAX_SCORE))
            self.nuclhiscore = ABS_MAX_SCORE

        # Verify protein score values were legitimate
        if self.protloscore > self.prothiscore:
            sys.stderr.write('Warning: protloscore > prothiscore, swapping instead!??\n')
            garbage = self.protloscore
            self.protloscore = self.prothiscore
            self.prothiscore = garbage
        if self.protloscore < ABS_MIN_SCORE:
            sys.stderr.write('Warning: Bad --protloscore, using %s instead\n' \
                             % str(ABS_MIN_SCORE))
            self.protloscore = ABS_MIN_SCORE
        if self.prothiscore > ABS_MAX_SCORE:
            sys.stderr.write('Warning: Bad --prothiscore, using %s instead\n' \
                             % str(ABS_MAX_SCORE))
            self.prothiscore = ABS_MAX_SCORE

        # Verify nucleotide coverage values were legitimate
        if self.nucllocov > self.nuclhicov:
            sys.stderr.write('Warning: nucllocov > nuclhicov, swapping instead!??\n')
            garbage = self.nucllocov
            self.nucllocov = self.nuclhicov
            self.nuclhicov = garbage
        if self.nucllocov < ABS_MIN_COV:
            sys.stderr.write('Warning: Bad --nucllocov, using %s instead\n' % str(ABS_MIN_COV))
            self.nucllocov = ABS_MIN_COV
        if self.nuclhicov > ABS_MAX_COV:
            sys.stderr.write('Warning: Bad --nuclhicov, using %s instead\n' % str(ABS_MAX_COV))
            self.nuclhicov = ABS_MAX_COV

        # Verify protein coverage values were legitimate
        if self.protlocov > self.prothicov:
            sys.stderr.write('Warning: protlocov > prothicov, swapping instead!??\n')
            garbage = self.protlocov
            self.protlocov = self.prothicov
            self.prothicov = garbage
        if self.protlocov < ABS_MIN_COV:
            sys.stderr.write('Warning: Bad --protlocov, using %s instead\n' % str(ABS_MIN_COV))
            self.protlocov = ABS_MIN_COV
        if self.prothicov > ABS_MAX_COV:
            sys.stderr.write('Warning: Bad --prothicov, using %s instead\n' % str(ABS_MAX_COV))
            self.prothicov = ABS_MAX_COV

    def __xrefLookup(self,indfilename,chrp):
        '''Populate dictionary for storing numeric indices for genomic templates

        This function consults a lookup table (indfilename, a plain text file in the
        resident file system) and populates dictionaries for storing numeric indices
        based on Locus or Accession identifiers.

        The agreed-upon format for this text file is
            Alias {\t} Alias {\t} ID {\n}
        Where this is generally Locus \t Accession \t GI.
        The third field will be used as the indexed key.
        At least one of fields one and two must refer to the genomic
        template *as it is indexed in the gthxml output!*
        '''

        if indfilename == '/dev/null':
            raise xGDBgthxmlPGSError, \
              'Error: Must specify a cross reference file for genomic templates!'

        # We must distinguish alternate formats for the xref index file
        # based on whether a chromosome was used as the genomic reference
        # type or not.
        if chrp:
            regex = re.compile('^(\S+)\s+(\S+)\s+chr(\S+)$')
        else:
            regex = re.compile('^(\S+)\s+(\S+)\s+(\S+)$')

        self.xref_dictLOC2GI = {}
        self.xref_dictACC2GI = {}

        xref = open(indfilename,'r')
        for xrefline in xref:
            match = regex.match(xrefline)
            if match and not re.search(r'\D',match.group(3)):
                # populate the xref dictionaries
                self.xref_dictLOC2GI[match.group(1)] = int(match.group(3))
                self.xref_dictACC2GI[match.group(2)] = int(match.group(3))
            else:
                xref.close()
                raise xGDBgthxmlPGSError, \
                  'Error: Blatantly non-conforming crud in xref file!'
        xref.close()

    def __cpxrefLookup(self,cpindfilename,qtype):
        '''Populate cDNA or protein dictionary

        Populate dictionary for storing identifiers of full-length cDNA or
        protein sequences in the input file.  This data structure is useful
        for appropriate segregation of results based on mixtures of query
        types, e.g., an EST+cDNA mixture, an EST+cDNA+protein mixture, etc.

        The agreed-upon format for this text file is
            '^(\S+)$'
        where each record corresponds to the identification string
        (possibly GI, but not necessarily) of a cDNA or protein.
        '''

        if cpindfilename == '/dev/null':
            # This is a legitimate condition; self.cfilexrefp or
            # self.pfilexrefp are set to False by default.
            pass
        else:
            if qtype not in ('cDNA', 'prot'):
                raise xGDBgthxmlPGSError, \
                  'Error: Internal misuse of __cpxrefLookup method!'
            elif qtype == 'cDNA': 
                self.cfilexrefp = True
                self.xref_cdnas = {}
            else: # qtype == 'prot'
                self.pfilexrefp = True
                self.xref_prots = {}

            cpxref = open(cpindfilename,'r')
            cpregex = re.compile('^(\S+)$')

            for cpxrefline in cpxref:
                match = cpregex.match(cpxrefline)
                if match:
                    if qtype == 'cDNA':
                        # populate the xref_cdnas dictionary
                        self.xref_cdnas[match.group(1)] = True
                    else: # qtype == 'prot'
                        # populate the xref_prots dictionary
                        self.xref_prots[match.group(1)] = True
                else:
                    cpxref.close()
                    if qtype == 'cDNA':
                        raise xGDBgthxmlPGSError, \
                          'Error: Blatantly non-conforming crud in cdna xref file!'
                    else: # qtype == 'prot'
                        raise xGDBgthxmlPGSError, \
                          'Error: Blatantly non-conforming crud in prot xref file!'
            cpxref.close()

    def __reset(self):
        # Bools for triggering SAX events
        self.in_alnmodule = False
        self.in_exintinfo = False
        self.in_matchline = False
        self.in_totalnscore = False
        self.in_cumlen = False
        self.in_PGS = False
        self.in_alignment = False
        self.in_genome_strand = False
        self.in_mrna_strand = False
        self.in_prot_strand = False
        self.in_hqPGS = False

        # Does gene structure contain >= 1 intron?
        self.intronp = False

        # strings used to process alignments
        self.genome_strand = ''
        self.mrna_strand = ''
        self.prot_strand = ''

        # Stores query type: cDNA or Protein
        self.querytype = ''

        # Lists for storing exon and intron data
        self.exon_list = []
        self.intron_list = []

        # 15 dictionary fields for main
        # xGDB alignment table
        #
        # Okay, so you might notice that there are
        # actually 16 fields.  This is because xGDB
        # arbitrarily uses a field named "chr" in
        # tables storing alignments involving
        # assembled pseudomolecules, and "gseg_gi"
        # in those involving BAC/GSS assembly/etc.-type
        # genomic templates.  They are the same data
        # type in the backend, just named different.
        self.mainspaln = { 'uid'       : 'NULL', \
                           'gi'        : '', \
                           'E_O'       : '', \
                           'sim'       : '', \
                           'mlength'   : '', \
                           'cov'       : '', \
                           'gseg_gi'   : '', \
                           'chr'       : '', \
                           'G_O'       : '', \
                           'l_pos'     : '', \
                           'r_pos'     : '', \
                           'pgs'       : '', \
                           'pgs_lpos'  : '', \
                           'pgs_rpos'  : '', \
                           'gseg_gaps' : '', \
                           'pgs_gaps'  : ''  \
                         }

        # 7 dictionary fields for exon xGDB alignment table
        self.exonspaln = { 'pgs_uid'    : 'LAST_INSERT_ID()', \
                           'num'        : '', \
                           'gseg_start' : '', \
                           'gseg_stop'  : '', \
                           'pgs_start'  : '', \
                           'pgs_stop'   : '', \
                           'score'      : ''  \
                         }

        # 8 dictionary fields for intron xGDB alignment table
        self.intronspaln = { 'pgs_uid'    : 'LAST_INSERT_ID()', \
                             'num'        : '', \
                             'gseg_start' : '', \
                             'gseg_stop'  : '', \
                             'Dscore'     : '', \
                             'Dsim'       : '', \
                             'Ascore'     : '', \
                             'Asim'       : ''  \
                           }

    # Override startElement method from ContentHandler
    def startElement(self,name,attrs):
        if name == 'alignment_module':
            self.in_alnmodule = True
        if self.in_alnmodule:
            if name == 'exon-intron_info':
                self.in_exintinfo = True
            if self.in_exintinfo:
                if name == 'exon':
                    self.currserial = int(attrs['e_serial'])
                    self.exon_list.append([]) # Note: self.exon_list is set to []
                                              #       in the __reset method
                if name == 'gDNA_exon_boundary':
                    self.exon_list[self.currserial - 1].extend(( \
                      attrs['g_start'], \
                      attrs['g_stop']))
                if name == 'reference_exon_boundary':
                    self.querytype = attrs['r_type']
                    self.exon_list[self.currserial - 1].extend(( \
                      attrs['r_start'], \
                      attrs['r_stop'], \
                      attrs['r_score']))
                if name == 'intron':
                    self.intronp = True
                    self.currserial = int(attrs['i_serial'])
                    self.intron_list.append([]) # Note: self.intron_list is set to []
                                                #       in the __reset method
                if name == 'gDNA_intron_boundary':
                    self.intron_list[self.currserial - 1].extend(( \
                      attrs['i_start'], \
                      attrs['i_stop']))
                if name == 'donor':
                    # One exon must come before any introns, so this test is safe
                    if self.querytype == 'cDNA':
                        self.intron_list[self.currserial - 1].extend(( \
                          attrs['d_prob'], \
                          attrs['d_score']))
                    else:
                        # The Dsim field, as used in the xGDB table structure, is
                        # not defined for protein alignments--I set it to a default,
                        # ``impossible" value of -1.0 (see self.scoreDNE).
                        self.intron_list[self.currserial - 1].extend(( \
                          attrs['d_prob'], \
                          self.scoreDNE))
                if name == 'acceptor':
                    # One exon must come before any introns, so this test is safe
                    if self.querytype == 'cDNA':
                        self.intron_list[self.currserial - 1].extend(( \
                          attrs['a_prob'], \
                          attrs['a_score']))
                    else:
                        # The Asim field, as used in the xGDB table structure, is
                        # not defined for protein alignments--I set it to a default,
                        # ``impossible" value of -1.0 (see self.scoreDNE).
                        self.intron_list[self.currserial - 1].extend(( \
                          attrs['a_prob'], \
                          self.scoreDNE))
            if name == 'MATCH_line':
                self.mainspaln['gi'] = attrs['ref_id']
                self.in_matchline = True
            if self.in_matchline:
                if name == 'total_alignment_score':
                    self.in_totalnscore = True
                if name == 'cumulative_length_of_scored_exons':
                    self.in_cumlen = True
                if name == 'coverage':
                    self.mainspaln['cov'] = attrs['percentage']
            if name == 'PGS_line':
                self.in_PGS = True
            if self.in_PGS:
                if name == 'gDNA':
                    if self.chrp: # access pseudomolecule xGDB table
                        if self.xref_dictLOC2GI.has_key(attrs['gen_id']):
                            self.mainspaln['chr'] = self.xref_dictLOC2GI[attrs['gen_id']]
                        elif self.xref_dictACC2GI.has_key(attrs['gen_id']):
                            self.mainspaln['chr'] = self.xref_dictACC2GI[attrs['gen_id']]
                        else:
                            raise xGDBgthxmlPGSError, \
                              'Error: xreffile does not index a given gDNA!'
                    else: # other genomic template types
                        if self.xref_dictLOC2GI.has_key(attrs['gen_id']):
                            self.mainspaln['gseg_gi'] = self.xref_dictLOC2GI[attrs['gen_id']]
                        elif self.xref_dictACC2GI.has_key(attrs['gen_id']):
                            self.mainspaln['gseg_gi'] = self.xref_dictACC2GI[attrs['gen_id']]
                        else:
                            raise xGDBgthxmlPGSError, \
                              'Error: xreffile does not index a given gDNA!'
                    self.mainspaln['G_O'] = attrs['gen_strand']
                if name == 'rDNA':
                    self.mainspaln['eseg_gi'] = attrs['rDNA_id']
                    self.mainspaln['E_O'] = attrs['rDNA_strand']
                if name == 'rProt':
                    self.mainspaln['eseg_gi'] = attrs['rProt_id']
                    self.mainspaln['E_O'] = '+' # Reasonable default! :)
                if name == 'exon':
                    self.mainspaln['pgs'] = \
                      self.mainspaln['pgs'] + \
                      attrs['e_start'] + \
                      '  ' + \
                      attrs['e_stop'] + \
                      ','
            if name == 'alignment':
                self.in_alignment = True
            if self.in_alignment:
                if name == 'genome_strand':
                    self.in_genome_strand = True
                if name == 'mrna_strand':
                    self.in_mrna_strand = True
                if name == 'queryProt':
                    self.in_prot_strand = True
            if name == 'hqPGS_line' and self.gsqp: # autocorrelated spliced alignments
                # Processing these elements will overwrite the raw PGS data stored
                # in the mainspaln dictionary.  Support for hqPGS data from protein
                # alignments (not available in GeneSeqer or GenomeThreader, as of
                # version 1.2.1 of this script) is implemented to support this feature
                # if it should become available.
                self.in_hqPGS = True
            if self.in_hqPGS:
                if name == 'hqgDNA':
                    if self.chrp: # access pseudomolecule xGDB table
                        if self.xref_dictLOC2GI.has_key(attrs['gen_id']):
                            self.mainspaln['chr'] = self.xref_dictLOC2GI[attrs['gen_id']]
                        elif self.xref_dictACC2GI.has_key(attrs['gen_id']):
                            self.mainspaln['chr'] = self.xref_dictACC2GI[attrs['gen_id']]
                        else:
                            raise xGDBgthxmlPGSError, \
                              'Error: xreffile does not index a given gDNA!'
                    else: # other genomic template types
                        if self.xref_dictLOC2GI.has_key(attrs['gen_id']):
                            self.mainspaln['gseg_gi'] = self.xref_dictLOC2GI[attrs['gen_id']]
                        elif self.xref_dictACC2GI.has_key(attrs['gen_id']):
                            self.mainspaln['gseg_gi'] = self.xref_dictACC2GI[attrs['gen_id']]
                        else:
                            raise xGDBgthxmlPGSError, \
                              'Error: xreffile does not index a given gDNA!'
                    self.mainspaln['G_O'] = attrs['gen_strand']
                if name == 'hqrDNA':
                    self.mainspaln['eseg_gi'] = attrs['rDNA_id']
                    self.mainspaln['E_O'] = attrs['rDNA_strand']
                if name == 'hqrProt':
                    self.mainspaln['eseg_gi'] = attrs['rProt_id']
                    self.mainspaln['E_O'] = '+' # Reasonable default! :)
                if name == 'hqexon':
                    self.mainspaln['pgs'] = \
                      self.mainspaln['pgs'] + \
                      attrs['e_start'] + \
                      '  ' + \
                      attrs['e_stop'] + \
                      ','

    # Override characters method from ContentHandler
    # Note that the parser will break content up into chunks of 1024 bytes,
    # so it is safest (well, it's completely necessary!) to iterate over each
    # chunk (hence, the append (+) operator).
    def characters(self,characters):
        if self.in_totalnscore:
            self.mainspaln['sim'] = self.mainspaln['sim'] + str(characters)
        if self.in_cumlen:
            self.mainspaln['mlength'] = self.mainspaln['mlength'] + str(characters)
        if self.in_genome_strand:
            self.genome_strand = self.genome_strand + str(characters)
        if self.in_mrna_strand:
            self.mrna_strand = self.mrna_strand + str(characters)
        if self.in_prot_strand:
            self.prot_strand = self.prot_strand + str(characters)

    # Override endElement method from ContentHandler
    def endElement(self,name):
        if self.in_alnmodule:
            if name == 'exon-intron_info': # set a value or two
                self.mainspaln['l_pos'] = min( \
                  int(self.exon_list[0][0]), \
                  int(self.exon_list[len(self.exon_list) - 1][1]))
                self.mainspaln['r_pos'] = max( \
                  int(self.exon_list[0][0]), \
                  int(self.exon_list[len(self.exon_list) - 1][1]))
                self.mainspaln['pgs_lpos'] = min( \
                  int(self.exon_list[0][2]), \
                  int(self.exon_list[len(self.exon_list) - 1][3]))
                self.mainspaln['pgs_rpos'] = max( \
                  int(self.exon_list[0][2]), \
                  int(self.exon_list[len(self.exon_list) - 1][3]))
                self.in_exintinfo = False
            if name == 'total_alignment_score':
                self.in_totalnscore = False
            if name == 'cumulative_length_of_scored_exons':
                self.in_cumlen = False
            if name == 'MATCH_line':
                self.in_matchline = False
            if name == 'PGS_line' and self.in_PGS:
                self.mainspaln['pgs'] = self.mainspaln['pgs'][:-1] # chop trailing ','
                self.in_PGS = False
            if name == 'genome_strand':
                self.in_genome_strand = False
            if name == 'mrna_strand':
                self.in_mrna_strand = False
            if name == 'queryProt':
                self.in_prot_strand = False
            if name == 'alignment':
                self.in_alignment = False
            if name == 'hqPGS_line' and self.in_hqPGS:
                self.mainspaln['pgs'] = self.mainspaln['pgs'][:-1] # chop trailing ','
                self.in_hqPGS = False
            if name == 'predicted_gene_structure' and \
               self.in_alnmodule: # Triggers match recording
                # First, let's verify we have the query type set reasonably;
                # we will need this information in the methods called below.
                if self.querytype not in ('cDNA', 'Protein'):
                    raise xGDBgthxmlPGSError, \
                      'Error: Invalid query type! (You did not validate your gthxml?!!)'
                self.__compress_alignment()
                self.__recordmatch()
                self.__reset()
                self.in_alnmodule = True # not necessarily finished...
            if name == 'alignment_module': # Prevent processing data in PGL_module of gthxml
                self.in_alnmodule = False

    def __compress_alignment(self):
        '''Compress the alignment information using the xGDB routine.'''

        if self.querytype == 'cDNA' and \
           (len(self.genome_strand) != len(self.mrna_strand)):
            raise xGDBgthxmlPGSError, \
              'Error: Genomic and EST/cDNA alignment strings incompatible!'
        elif self.querytype == 'Protein' and \
             (len(self.genome_strand) != len(self.prot_strand)):
            raise xGDBgthxmlPGSError, \
              'Error: Genomic and Protein alignment strings incompatible!'
        else:
            pass

        # Process genomic template first
        self.genomic_Watson = \
          (int(self.exon_list[0][0]) < \
           int(self.exon_list[len(self.exon_list) - 1][1])) and True or False
        in_gap = False
        genomicposition = int(self.exon_list[0][0])
        for i in range(len(self.genome_strand)):
            if self.genome_strand[i] == '-':
                if in_gap:
                    len_gap += 1
                else:
                    in_gap = True
                    len_gap = 1
                    start_gap = genomicposition
            elif re.search(r'\w',self.genome_strand[i]):
                if in_gap:
                    if self.genomic_Watson:
                        self.mainspaln['gseg_gaps'] = \
                          self.mainspaln['gseg_gaps'] + \
                          ':' + str(genomicposition - 1) + \
                          ':' + str(len_gap)
                    else:
                        self.mainspaln['gseg_gaps'] = \
                          self.mainspaln['gseg_gaps'] + \
                          ':' + str(genomicposition) + \
                          ':' + str(len_gap)

                in_gap = False

                # update genomicposition
                if self.genomic_Watson:
                    genomicposition += 1
                else:
                    genomicposition -= 1

            else:
                raise xGDBgthxmlPGSError, \
                  'Error: Encountered invalid symbol in genomic template.'

        self.mainspaln['gseg_gaps'] = self.mainspaln['gseg_gaps'][1:] # chop leading ':'

        # Now, process the query sequence

        # Code to handle cDNA alignments
        if self.querytype == 'cDNA':
            in_acc_gap = in_gap = in_intron = False
            queryposition = int(self.exon_list[0][2])
            self.query_Watson = \
              (int(self.exon_list[0][2]) < \
               int(self.exon_list[len(self.exon_list) - 1][3])) and True or False
            for i in range(len(self.mrna_strand)):
                if self.mrna_strand[i] == '-':
                    if in_intron: # gap adjacent to acceptor
                        if in_acc_gap:
                            len_acc_gap += 1
                        else:
                            in_acc_gap = True
                            len_acc_gap = 1
                            start_acc_gap = queryposition + 1
                    else:
                        if in_gap:
                            len_gap += 1
                        else:
                            in_gap = True
                            len_gap = 1
                            start_gap = queryposition
                elif self.mrna_strand[i] == '.':
                    if in_intron:
                        len_intron += 1
                    else:
                        in_intron = True
                        len_intron = 1
                        start_intron = queryposition
                elif re.search(r'\w',self.mrna_strand[i]):
                    if in_gap:
                        if self.query_Watson:
                            self.mainspaln['pgs_gaps'] = \
                              self.mainspaln['pgs_gaps'] + \
                              ':' + str(start_gap - 1) + \
                              ':' + str(len_gap)
                        else:
                            self.mainspaln['pgs_gaps'] = \
                              self.mainspaln['pgs_gaps'] + \
                              ':' + str(start_gap) + \
                              ':' + str(len_gap)
                    if in_acc_gap:
                        if self.query_Watson:
                            self.mainspaln['pgs_gaps'] = \
                              self.mainspaln['pgs_gaps'] + \
                              ':' + str(start_acc_gap - 1) + \
                              ':' + str(-(len_acc_gap))
                        else:
                            self.mainspaln['pgs_gaps'] = \
                              self.mainspaln['pgs_gaps'] + \
                              ':' + str(start_acc_gap) + \
                              ':' + str(-(len_acc_gap))
                    if in_intron:
                        if self.query_Watson:
                            self.mainspaln['pgs_gaps'] = \
                              self.mainspaln['pgs_gaps'] + \
                              ':' + str(-(start_intron - 1)) + \
                              ':' + str(len_intron)
                        else:
                            self.mainspaln['pgs_gaps'] = \
                              self.mainspaln['pgs_gaps'] + \
                              ':' + str(-(start_intron)) + \
                              ':' + str(len_intron)
    
                    in_acc_gap = in_gap = in_intron = False
    
                    # update queryposition
                    if self.query_Watson:
                        queryposition += 1
                    else:
                        queryposition -= 1
                else:
                    raise xGDBgthxmlPGSError, \
                      'Error: Encountered invalid symbol in query sequence: %s' % \
                      self.mrna_strand[i]

        # Code to handle protein alignments
        # Protein alignments are similar to transcript ones, but must
        # take phase relative to the genomic template into account.
        # Also, proteins are exclusively in forward orientation.
        #
        # gap/intron phase conventions:
        #     a ->  C O D - C O D
        #     b ->  C - O D C O D
        #     c ->  C O - D C O D
        #
        # Gap/intron start positions are given in terms of the query
        # sequence (amino acids), whereas lengths of the same are
        # given in terms of the genomic template (nucleotides).
        elif self.querytype == 'Protein':
            in_acc_gap = in_gap = in_intron = was_in_intron = False
            queryposition = int(self.exon_list[0][2])
            start_qprot = True  # Are we beginning the alignment parse?
            dist_from_Psite = 1 # Just the distance from most recent
                                # AA symbol in prot_strand, i.e.,
                                #      01201
                                #     GCAACC etc.
                                #      A  T
                                # 'Psite' comes from ribosomal EPA nomenclature.
                                # This is used to assist gap phase classification;
                                # it is always given in terms of where point is
                                # at a given moment, never where it's going or
                                # where it's been.
            safety_space_ct = 0 # Protect against gap symbols that may
                                # fail to print in acceptor gaps

            for i in range(len(self.prot_strand)):
                if self.prot_strand[i] == ' ':
                    if in_gap or in_acc_gap:
                        pass
                    elif in_intron:
                        # Identify the intron's phase
                        gap_phase = self.phase_symbols[ (dist_from_Psite + 1) % 3 ]
                        self.mainspaln['pgs_gaps'] = \
                          self.mainspaln['pgs_gaps'] + \
                          ':' + str(-(start_intron - 1)) + gap_phase + \
                          ':' + str(len_intron)
                        in_intron = False
                        was_in_intron = True
                        safety_space_ct = 1
                    elif was_in_intron:
                        safety_space_ct += 1
                    else:
                        dist_from_Psite += 1

                elif self.prot_strand[i] == '-':
                    # I assume that '.-' never occurs in the query alignment string.(?)
                    if was_in_intron and not in_acc_gap: # gap adjacent to acceptor
                        in_acc_gap = True
                        len_acc_gap = 3
                        start_acc_gap = queryposition + 1
                        was_in_intron = False
                    elif in_acc_gap:
                        # I assume that all gaps in the protein query incur
                        # gap lengths of multiples of 3.(?)
                        len_acc_gap += 3
                    else:
                        # I assume that all gaps in the protein query incur
                        # gap lengths of multiples of 3.(?)
                        if in_gap:
                            len_gap += 3
                        else:
                            in_gap = True
                            len_gap = 3
                            start_gap = queryposition

                elif self.prot_strand[i] == '.':
                    if start_qprot:
                        # The example below, corresponding to this case,
                        # is a special one, where the first codon position
                        # is in an upstream exon of sufficiently poor quality
                        # as not to be included in the final alignment.
                        #
                        # If the residue A in the query protein below were
                        # at position 114, this intron would be represented
                        # as '-113b:1'.
                        #
                        #    GCAACCCTTG     genome
                        #     A  T  L       translated genome
                        #     |  |  .  
                        #    .A  T  F       query protein
                        #
                        # Note that I assume also that such ``introductory
                        # introns" do not occur in phases a or c, but only b.(?)
                        in_intron = True
                        len_intron = 1
                        start_intron = queryposition
                        dist_from_Psite = 2 # Special trick to get the phase (b) correct.
                    elif in_intron:
                        len_intron += 1
                    else:
                        in_intron = True
                        len_intron = 1
                        start_intron = queryposition
                        if not in_gap:
                            dist_from_Psite += 1
                            dist_from_Psite %= 3

                elif re.search(r'[\w\*]',self.prot_strand[i]): # '*' is Stop codon
                    if in_gap:
                        # Identify the gap's phase
                        gap_phase = self.phase_symbols[ (dist_from_Psite + 1) % 3 ]
                        self.mainspaln['pgs_gaps'] = \
                          self.mainspaln['pgs_gaps'] + \
                          ':' + str(start_gap - 1) + gap_phase + \
                          ':' + str(len_gap)
                    if in_acc_gap:
                        # Identify the gap's phase
                        gap_phase = self.phase_symbols[ (dist_from_Psite + 1) % 3 ]
                        self.mainspaln['pgs_gaps'] = \
                          self.mainspaln['pgs_gaps'] + \
                          ':' + str(start_acc_gap - 1) + gap_phase + \
                          ':' + str(-(len_acc_gap))
                    if in_intron:
                        # Identify the intron's phase
                        gap_phase = self.phase_symbols[ (dist_from_Psite + 1) % 3 ]
                        self.mainspaln['pgs_gaps'] = \
                          self.mainspaln['pgs_gaps'] + \
                          ':' + str(-(start_intron - 1)) + gap_phase + \
                          ':' + str(len_intron)
                    in_acc_gap = in_gap = in_intron = was_in_intron = False
                    dist_from_Psite = safety_space_ct = 0
                    queryposition += 1

                else:
                    raise xGDBgthxmlPGSError, \
                      'Error: Encountered invalid symbol in query sequence: %s' % \
                      self.prot_strand[i]

                if dist_from_Psite > 2 or safety_space_ct > 2:
                    raise xGDBgthxmlPGSError, \
                      'Error: Bad query alignment string (too many contiguous spaces!): \n\
 See position %i of \n%s' % \
                      (i+1, self.prot_strand)
                start_qprot = False

        else:
            raise xGDBgthxmlPGSError, \
              'Error: (Internal) Invalid query type after verification?!!'

        self.mainspaln['pgs_gaps'] = self.mainspaln['pgs_gaps'][1:] # chop leading ':'

    def __recordmatch(self):
        '''Record SQL statements for populating xGDB backend.'''

        if (self.querytype == 'cDNA' and \
            float(self.mainspaln['sim']) >= self.nuclloscore and \
            float(self.mainspaln['sim']) <= self.nuclhiscore and \
            float(self.mainspaln['cov']) >= self.nucllocov and \
            float(self.mainspaln['cov']) <= self.nuclhicov) \
           or \
           (self.querytype == 'Protein' and \
            float(self.mainspaln['sim']) >= self.protloscore and \
            float(self.mainspaln['sim']) <= self.prothiscore and \
            float(self.mainspaln['cov']) >= self.protlocov and \
            float(self.mainspaln['cov']) <= self.prothicov):

            # Populate the main spliced alignment table
            if self.querytype == 'Protein':
                if not (self.pfilexrefp and \
                        self.xref_prots.has_key(self.mainspaln['gi'])):
                    raise xGDBgthxmlPGSError, \
                      'Error: Protein %s not documented in xref file!' % \
                      self.mainspaln['gi']
                else:
                    sqlstring = 'INSERT INTO ' + self.ptable
            elif self.querytype == 'cDNA' and \
                 self.cfilexrefp and \
                 self.xref_cdnas.has_key(self.mainspaln['gi']):
                sqlstring = 'INSERT INTO ' + self.ctable
            else:
                sqlstring = 'INSERT INTO ' + self.table

            if self.chrp:
                sqlstring = \
                  sqlstring + \
                  ' (uid,gi,E_O,sim,mlength,cov,chr,G_O,l_pos,r_pos,' + \
                  'pgs,pgs_lpos,pgs_rpos,gseg_gaps,pgs_gaps) VALUES (' + \
                  self.mainspaln['uid'] + ',' + \
                  '\'' + self.mainspaln['gi'] + '\',' + \
                  '\'' + self.mainspaln['E_O'] + '\',' + \
                  self.mainspaln['sim'] + ',' + \
                  self.mainspaln['mlength'] + ',' + \
                  self.mainspaln['cov'] + ',' + \
                  str(self.mainspaln['chr']) + ','
            else:
                sqlstring = \
                  sqlstring + \
                  ' (uid,gi,E_O,sim,mlength,cov,gseg_gi,G_O,l_pos,r_pos,' + \
                  'pgs,pgs_lpos,pgs_rpos,gseg_gaps,pgs_gaps) VALUES (' + \
                  self.mainspaln['uid'] + ',' + \
                  '\'' + self.mainspaln['gi'] + '\',' + \
                  '\'' + self.mainspaln['E_O'] + '\',' + \
                  self.mainspaln['sim'] + ',' + \
                  self.mainspaln['mlength'] + ',' + \
                  self.mainspaln['cov'] + ',' + \
                  str(self.mainspaln['gseg_gi']) + ','

            sqlstring = \
              sqlstring + \
              '\'' + self.mainspaln['G_O'] + '\',' + \
              str(self.mainspaln['l_pos']) + ',' + \
              str(self.mainspaln['r_pos']) + ',' + \
              '\'' + self.mainspaln['pgs'] + '\',' + \
              str(self.mainspaln['pgs_lpos']) + ',' + \
              str(self.mainspaln['pgs_rpos']) + ',' + \
              '\'' + self.mainspaln['gseg_gaps'] + '\',' + \
              '\'' + self.mainspaln['pgs_gaps'] + '\');'

            self.outstream.write(sqlstring + '\n')

            # Now, for the exons table
            if self.querytype == 'Protein':
                if not (self.pfilexrefp and \
                        self.xref_prots.has_key(self.mainspaln['gi'])):
                    raise xGDBgthxmlPGSError, \
                      'Error: Protein %s not documented in xref file!' % \
                      self.mainspaln['gi']
                else:
                    sqlstring = 'INSERT INTO ' + self.ptable
            elif self.querytype == 'cDNA' and \
                 self.cfilexrefp and \
                 self.xref_cdnas.has_key(self.mainspaln['gi']):
                sqlstring = 'INSERT INTO ' + self.ctable
            else:
                sqlstring = 'INSERT INTO ' + self.table

            sqlstring = sqlstring + \
              '_exons (pgs_uid,num,gseg_start,gseg_stop,pgs_start,pgs_stop,score) VALUES '
            for i in range(len(self.exon_list)):
                self.exonspaln['num'] = str(i + 1)
                self.exonspaln['gseg_start'] = self.exon_list[i][0]
                self.exonspaln['gseg_stop'] = self.exon_list[i][1]
                self.exonspaln['pgs_start'] = self.exon_list[i][2]
                self.exonspaln['pgs_stop'] = self.exon_list[i][3]
                self.exonspaln['pgs_score'] = self.exon_list[i][4]

                sqlstring = \
                  sqlstring + '(' + \
                  self.exonspaln['pgs_uid'] + ',' + \
                  self.exonspaln['num'] + ',' + \
                  self.exonspaln['gseg_start'] + ',' + \
                  self.exonspaln['gseg_stop'] + ',' + \
                  self.exonspaln['pgs_start'] + ',' + \
                  self.exonspaln['pgs_stop'] + ',' + \
                  self.exonspaln['pgs_score'] + '),'

            sqlstring = sqlstring[:-1] + ';' # replace trailing ',' with ';'
            self.outstream.write(sqlstring + '\n')

            # Now, for the introns table, if applicable
            if self.intronp:
                if self.querytype == 'Protein':
                    if not (self.pfilexrefp and \
                            self.xref_prots.has_key(self.mainspaln['gi'])):
                        raise xGDBgthxmlPGSError, \
                          'Error: Protein %s not documented in xref file!' % \
                          self.mainspaln['gi']
                    else:
                        sqlstring = 'INSERT INTO ' + self.ptable
                elif self.querytype == 'cDNA' and \
                     self.cfilexrefp and \
                     self.xref_cdnas.has_key(self.mainspaln['gi']):
                    sqlstring = 'INSERT INTO ' + self.ctable
                else:
                    sqlstring = 'INSERT INTO ' + self.table

                sqlstring = sqlstring + \
                  '_introns (pgs_uid,num,gseg_start,gseg_stop,Dscore,Dsim,Ascore,Asim) VALUES '
                for i in range(len(self.intron_list)):
                    self.intronspaln['num'] = str(i + 1)
                    self.intronspaln['gseg_start'] = self.intron_list[i][0]
                    self.intronspaln['gseg_stop'] = self.intron_list[i][1]
                    self.intronspaln['Dscore'] = self.intron_list[i][2]
                    self.intronspaln['Dsim'] = self.intron_list[i][3]
                    self.intronspaln['Ascore'] = self.intron_list[i][4]
                    self.intronspaln['Asim'] = self.intron_list[i][5]

                    sqlstring = \
                      sqlstring + '(' + \
                      self.intronspaln['pgs_uid'] + ',' + \
                      self.intronspaln['num'] + ',' + \
                      self.intronspaln['gseg_start'] + ',' + \
                      self.intronspaln['gseg_stop'] + ',' + \
                      self.intronspaln['Dscore'] + ',' + \
                      self.intronspaln['Dsim'] + ',' + \
                      self.intronspaln['Ascore'] + ',' + \
                      self.intronspaln['Asim'] + '),'

                sqlstring = sqlstring[:-1] + ';' # replace trailing ',' with ';'
                self.outstream.write(sqlstring + '\n')

    def __del__(self):
        '''Class "destructor" method'''

        # close input stream
        if self.instreamp:
            self.instream.close()
            if self.webp:
                self.compressed_gthxml_stream.close()

        # close output stream
        if self.outstreamp:
            self.outstream.close()

if __name__ == '__main__':
    gthxml_handler = xGDBgthxmlPGS(sys.argv[1:])

    from xml.sax import make_parser
    from xml.sax.handler import ContentHandler

    genericsax = make_parser()
    genericsax.setContentHandler(gthxml_handler)
    try:
        genericsax.parse(gthxml_handler.instream)
    except xGDBgthxmlPGSError, inst:
        print inst.args[0]
