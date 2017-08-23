#!/usr/bin/env python
'''Convert spliced alignments in gthxml format to GFF

This module contains the gthxmlToGFF class definition.
It parses results stored in the gthxml format, as
defined at http://www.public.iastate.edu/~mespar1/gthxml/GenomeThreader.rng.txt .
(Note that GenomeThreader can generate this output directly using the
-xmlout option, whereas GeneSeqer plain text output can be
retrofitted to this standard via the GSQ2XML.pl tool, available
at http://www.public.iastate.edu/~mespar1/gthxml/GSQ2XML.pl .)

Usage:
    python gthxmlToGFF.py (aka XML2GFF.py)
        -h --help           # prints this message and exits
        -t --type qtype     # one of "EST", "cDNA", "prot", or "mixed"
        -i --infile file    # Input file identifier. Input stream can be gthxml data
                            # piped by stdin (default) or a file on the resident file
                            # system
        -o --outfile file   # Output file identifier. Can be stdout (default) or a
                            # local file, which will be overwritten if it exists
        --loscore float     # minimum overall alignment score
        --hiscore float     # maximum overall alignment score
        --locov float       # minimum coverage of G,C,or P
        --hicov float       # maximum coverage of G,C,or P
'''

__author__ = 'Michael E Sparks (mespar1@iastate.edu)'
__version__ = '1.2'
__copyright__ = 'Copyright (c) 2004 Michael E Sparks, all rights reserved'
__license__ = 'BSD-style license'

import getopt, sys
from xml.sax.handler import ContentHandler

ABS_MIN_SCORE = 0.00
ABS_MAX_SCORE = 2.00 # I don't like this, but it works...mes
ABS_MIN_COV = 0.00
ABS_MAX_COV = 2.00 # I don't like this, but it works...mes
QTYPES = ('mixed', 'EST', 'cDNA', 'prot') # If you do change this, please
                                          # leave 'mixed' at the front!

class gthxmlToGFFError(Exception):
    '''Exception subclass for addressing fatal errors in xGDBgthxmlPGS objects'''
    pass
                                                                                                                             
class gthxmlToGFF(ContentHandler):
    '''A handler to convert spliced alignments in gthxml to GFF

    Allows one to parse a gthxml v1.0 file from a local filesystem or over the
    Web and print spliced alignments, per minimum/maximum score and coverage 
    criteria, in GFF formats for use in, e.g., gbrowse
    '''
    def __init__(self,argv):
        '''Class "constructor" method'''
                                                                                                                             
        # Set default values prior to parsing argument list
        # Note that these may trigger exceptions later in the class
        # methods, but this is the intended effect!  In short, don't
        # monkey with this part of the code.
        self.infile = sys.stdin      # Input file identifier. Input stream can be gthxml data
                                     # piped by stdin or a file on the resident file system.
        self.instreamp = False       # Have we opened the input stream yet?
        self.outfile = sys.stdout    # Output file identifier. Can be stdout or a local file,
                                     # which will be overwritten if it exists
        self.outstreamp = False      # Have we opened the output stream yet?
        self.loscore = 0.80          # minimum alignment score
        self.hiscore = ABS_MAX_SCORE # maximum alignment score
        self.locov = 0.80            # min coverage of G,C,or P
        self.hicov = ABS_MAX_COV     # max coverage of G,C,or P
        self.qtype = 'mixed'         # type of query used for spliced alignment
                                                                                                                             
        # Parse option list
        try:
            opts, args = getopt.getopt(argv, 'ht:i:o:',
                                            ['help',
                                             'type=',
                                             'infile=',
                                             'outfile=',
                                             'loscore=', 'hiscore=',
                                             'locov=', 'hicov='])

        except getopt.GetoptError:
            print __doc__
            sys.exit(2)
                                                                                                                             
        for o, a in opts:
            if o in ('-h', '--help'):
                print __doc__
                sys.exit(None)
            elif o in ('-t', '--type'):
                self.qtype = a
            elif o in ('-i', '--infile'):
                self.infile = a
            elif o in ('-o', '--outfile'):
                self.outfile = a
            elif o == '--loscore':
                self.loscore = float(a)
            elif o == '--hiscore':
                self.hiscore = float(a)
            elif o == '--locov':
                self.locov = float(a)
            elif o == '--hicov':
                self.hicov = float(a)

        # Verify score values were legitimate
        if self.loscore > self.hiscore:
            sys.stderr.write('Warning: loscore > hiscore, swapping instead!??\n')
            garbage = self.loscore
            self.loscore = self.hiscore
            self.hiscore = garbage
        if self.loscore < ABS_MIN_SCORE:
            sys.stderr.write('Warning: Bad --loscore, using %s instead\n' % str(ABS_MIN_SCORE))
            self.loscore = ABS_MIN_SCORE
        if self.hiscore > ABS_MAX_SCORE:
            sys.stderr.write('Warning: Bad --hiscore, using %s instead\n' % str(ABS_MAX_SCORE))
            self.hiscore = ABS_MAX_SCORE
                                                                                                                             
        # Verify coverage values were legitimate
        if self.locov > self.hicov:
            sys.stderr.write('Warning: locov > hicov, swapping instead!??\n')
            self.garbage = self.locov
            self.locov = self.hicov
            self.hicov = self.garbage
        if self.locov < ABS_MIN_COV:
            sys.stderr.write('Warning: Bad --locov, using %s instead\n' % str(ABS_MIN_COV))
            self.locov = ABS_MIN_COV
        if self.hicov > ABS_MAX_COV:
            sys.stderr.write('Warning: Bad --hicov, using %s instead\n' % str(ABS_MAX_COV))
            self.hicov = ABS_MAX_COV

        # Verify query type is legitimate
        if self.qtype not in QTYPES:
            sys.stderr.write('Warning: Incorrect qtype, using %s instead\n' % QTYPES[0])
            self.qtype = QTYPES[0]

        # Attempt to open input stream
        if self.infile is not sys.stdin:
            self.instream = open(self.infile,'r')
        else:
            self.instream = self.infile
        self.instreamp = True

        # Attempt to open output stream
        if self.outfile is not sys.stdout:
            self.outstream = open(self.outfile,'w')
        else:
            self.outstream = self.outfile
        self.outstreamp = True

        # Bools for triggering SAX events
        self.in_header = False
        self.in_splaln = False
        self.in_totscor = False
        self.is_start = False
        self.qtypeset = False

        # 9 fields for GFF dictionary
        self.GFF = { 'template': '', \
                     'source'  : '', \
                     'type'    : '', \
                     'start'   : '', \
                     'stop'    : '', \
                     'score'   : '', \
                     'strand'  : '', \
                     'phase'   : '.', \
                     'group'   : '' \
                   }

        # Store current (though temporary!) query type designations
        # if processing gthxml with mixed query type
        self.tempqtype = ''

        # Direction components
        self.dir_exon = ''

        # Store critical values for whether we print the alignment or not
        self.score = ''
        self.coverage = ''

        # Count exons
        self.exon_ct = 1

        # Start and stop of aligned regions
        self.gene_start = 0
        self.gene_stop = 0

        # String to enqueue
        self.pre_queue = ''

        # List of exon lines, possibly a PPA line
        self.exon_queue = []

    def startElement(self, name, attrs):
        if name == 'header':
            self.in_header = True
        if name == 'spliced_alignment':
            self.exon_queue = []
            self.in_splaln = True
        if self.in_header and name == 'source':
            if self.qtype != 'mixed':
                self.GFF['source'] = attrs['program'] + self.qtype
                self.qtypeset = True
            else:
                self.GFF['source'] = attrs['program']
        if self.in_splaln:
            if name == 'reference':
                self.GFF['group'] = 'transcript \"' + attrs['ref_id'] + \
                                    '\" ; Note \"' + attrs['ref_id']
            if name == 'template':
                self.GFF['template'] = attrs['temp_id']
                self.GFF['strand'] = attrs['temp_strand']
            if name == 'predicted_gene_structure':
                self.is_start = True
                self.exon_ct = 1
            if name == 'gDNA_exon_boundary':
                # Record exon start, stop
                self.GFF['start'] = str(min(int(attrs['g_start']), \
                                            int(attrs['g_stop'])))
                self.GFF['stop'] = str(max(int(attrs['g_start']), \
                                           int(attrs['g_stop'])))

                # Record overall gene start
                if self.is_start:
                    if self.GFF['strand'] == '+':
                        self.gene_start = self.GFF['start']
                    else:
                        self.gene_start = self.GFF['stop']
                    self.is_start = False

                # Record overall gene stop for each exon...
                # eventually it will reach its correct value.
                if self.GFF['strand'] == '+':
                    self.gene_stop = self.GFF['stop']
                else:
                    self.gene_stop = self.GFF['start']
            if name == 'reference_exon_boundary':
                # Determine query type, if necessary
                if not self.qtypeset:
                    if attrs['r_type'] == 'Protein':
                        self.tempqtype = 'prot'
                    else:
                        self.tempqtype = 'EST' # Reasonable default!
                    self.qtypeset = True
                self.GFF['score'] = attrs['r_score']
                # Ready to enqueue an exon line!
                self.pre_queue = self.GFF['template'] + '\t' + \
                                 self.GFF['source'] + \
                                 self.tempqtype + '\t' + \
                                 'exon\t' + \
                                 self.GFF['start'] + '\t' + \
                                 self.GFF['stop'] + '\t' + \
                                 self.GFF['score'] + '\t' + \
                                 self.GFF['strand'] + '\t' + \
                                 self.GFF['phase'] + '\t' + \
                                 self.GFF['group'] + '.exon' + \
                                 str(self.exon_ct) + '\"\n'
                self.exon_ct += 1
                self.exon_queue.append(self.pre_queue)
            if name == 'total_alignment_score':
                self.in_totscor = True
            if name == 'coverage':
                self.coverage = attrs['percentage']

    def characters(self, characters):
        if self.in_totscor:
            self.score = self.score + str(characters)
            self.in_totscor = False

    def endElement(self, name):
        if name == 'header':
            self.in_header = False
        if name == 'spliced_alignment':
            self.__recordmatch()
            self.in_splaln = False
            self.score = ''
            self.coverage = ''
            if self.qtype == 'mixed':
                self.qtypeset = False
                self.tempqtype = ''

    def __recordmatch(self):
        if float(self.score) >= self.loscore and \
           float(self.score) <= self.hiscore and \
           float(self.coverage) >= self.locov and \
           float(self.coverage) <= self.hicov:

            if self.GFF['strand'] == '+':
                self.dir_exon = self.gene_start + '\t' + \
                                self.gene_stop + '\t'
            else:
                self.dir_exon = self.gene_stop + '\t' + \
                                self.gene_start + '\t'
            # Print the transcript line
            self.outstream.write(self.GFF['template'] + '\t' + \
                                 self.GFF['source'] + \
                                 self.tempqtype + '\t' + \
                                 'transcript\t' + \
                                 self.dir_exon + \
                                 self.score + '\t' + \
                                 self.GFF['strand'] + '\t' + \
                                 self.GFF['phase'] + '\t' + \
                                 self.GFF['group'] + '\"\n')
            # Print the exon queue
            for x in self.exon_queue:
                self.outstream.write(x)

    def __del__(self):
        '''Class *destructor* method'''
                                                                                                                             
        # close input stream
        if self.instreamp:
            self.instream.close()

        # close output stream
        if self.outstreamp:
            self.outstream.close()

if __name__ == '__main__':
    gthxml_handler = gthxmlToGFF(sys.argv[1:])
                                                                                                                             
    from xml.sax import make_parser

    genericsax = make_parser()
    genericsax.setContentHandler(gthxml_handler)
    try:
        genericsax.parse(gthxml_handler.instream)
    except gthxmlToGFFError, inst:
        print inst.args[0]
