#!/usr/bin/perl -w

# GSQ2XML.pl

# Maintained by Michael E Sparks (mespar1@iastate.edu)
# Functionality last modified: 27 July 2006
#   (Script was updated to support gthXML v1.1)

# Copyright (c) 2003,2004,2005,2006 Michael E Sparks
# All rights reserved.

# You may distribute this software under the same terms as perl itself,
# described at http://www.perl.com/pub/a/language/misc/Artistic.html .
# Release under these terms effective starting 6 July 2006.

###########################################################################
# DESCRIPTION:                                                            #
# This script is a simple tool that can be used to convert GeneSeqer or   #
# GeneSeqer2 plain text output into the XML markup defined in             #
# http://www.genomethreader.org/GenomeThreader.rng.txt .  It can convert  #
# output corresponding to reference types of either ESTcDNA or Protein    #
# aligned to a genomic template, but not both.  Furthermore, it is        #
# doubtful that this particular functionality will ever be implemented... #
# by me, anyway.                                                          #
###########################################################################

# Import packages
use strict;
use Getopt::Std;

# Global "Constant" Variables
my $GTHXMLVERS = "1.1";
my $EXIT_ON_ERROR = 1;
my $HOPELESS_SITUATION = "Information cannot be resolved.";
my $MIN_FILES_TO_PARSE = 1; 
my ($TRUE,$FALSE) = (1,0);
my %ref_types = (
     "nucleic_acid" => "ESTcDNA",
     "protein"      => "Protein",
   );
my $SPARKS_AT_FAULT = 
"Please notify mespar1\@iastate.edu of an error in this script!";
my $USAGE = "
Usage: ./GSQ2XML.pl [-eE|-pP] [-sS (optional)] input_file(s)

Notes: 1) Please specify either -eE (ESTcDNA reference type)
                         --or-- -pP (protein reference type)
       2) If -sS is specified, your output will be directed to the
          STDOUT stream rather than a file; this is not a good idea
          if you are processing multiple input files!
";

# Main Application
Main: {
  # Determine type of reference files and output preference
  use vars qw($opt_e $opt_E $opt_p $opt_P $opt_s $opt_S);
  getopts('eEpPsS');

  # Pool responses for $opt_E (ESTcDNA reference type)
  if (($opt_e) || ($opt_E)) { $opt_E = $TRUE;  }
  else                      { $opt_E = $FALSE; }

  # Pool responses for $opt_P (Protein reference type)
  if (($opt_p) || ($opt_P)) { $opt_P = $TRUE;  }
  else                      { $opt_P = $FALSE; }

  # Pool responses for $opt_S (Output to STDOUT or file?)
  if (($opt_s) || ($opt_S)) { $opt_S = $TRUE;  }
  else                      { $opt_S = $FALSE; }

  # Determine overall reference type
  my $reference_type = "";
  if ($opt_E == $opt_P) { # XOR--either, but not both
    &fatal("$USAGE");
  }
  else {
    if    ($opt_E) { $reference_type = $ref_types{nucleic_acid}; }
    elsif ($opt_P) { $reference_type = $ref_types{protein};      }
    else           { &fatal($SPARKS_AT_FAULT);                   }
  }
    
  # Verify that the user specified at least $MIN_FILES_TO_PARSE
  # GSQ/GSQ2 plain text output files to process.  If so, process them...
  my @files = @ARGV;
  if ((@files < $MIN_FILES_TO_PARSE) && ($opt_S != $TRUE)) {
    &fatal("You must supply an input file(s) to parse!");
  }
  else {
    foreach my $file (@files) {
      my @ref_files = &find_all_reference_files($file);
    
      open(FIN, "<$file") or &fatal("Can't open $file for reading! $!");

      if ($opt_S == $FALSE) { # write to file(s)
        open(FOUT, ">${file}.xml") or &fatal(
          "Can't create ${file}.xml for writing! $!");
  
        &parse_header(\*FIN, \*FOUT, $reference_type, @ref_files);
        &parse_body(\*FIN, \*FOUT, $reference_type);

        close(FIN);
        close(FOUT); 
      }
      else { # write to STDOUT
        &parse_header(\*FIN, \*STDOUT, $reference_type, @ref_files);
        &parse_body(\*FIN, \*STDOUT, $reference_type);

        close(FIN);
      }
    } # end foreach
  }
} # end Main

###########################################################################

# Description: This subroutine issues the probable error
#              and aborts the program.
sub fatal {
  my $message = shift(@_);
  print STDERR "\a$message\n";
  exit($EXIT_ON_ERROR);
} # end fatal

###########################################################################

# Description: This subroutine writes the names of all unique reference
#              files (regardless of whether any significant alignments
#              mapped to them) to @ref_files and returns it.
sub find_all_reference_files {
  my $file = shift(@_);
  open(FIN, "<$file") or &fatal("Can't open $file for reading! $!");

  my @ref_files = ();

  while (my $line=<FIN>) {
    if ($line =~ m/library file/) {
      # If we found a library file line, record the unique name
      $line =~ m/:\s+([^;\s]+)[;]?/;
      my $ref_file = $1;
      my $is_unique = $TRUE; # Default value
      foreach my $previous (@ref_files) {
        if ($ref_file eq $previous) {
          $is_unique = $FALSE;
          last;
        }
      }
      if ($is_unique) {
        push(@ref_files, $ref_file);    
      }
    }
  } # end while
  close(FIN);
  return(@ref_files);
} # end find_all_reference_files

###########################################################################

# Description: Parses the content corresponding to the "header" element
#              as desribed in "GenomeThreader.RNG"
sub parse_header {
  my $in_fh = shift(@_);
  my $out_fh = shift(@_);
  my $reference_type = shift(@_);
  my @reference_files = @_;
  my ($BSSM,$PARMS)="";

  # Frankly, the code in the block below is quite self-explanatory,
  # so I won't add more commentary to it.
  HEADER_PARSE: {
    print $out_fh 
"<?xml version=\"1.0\" encoding=\"ISO-8859-1\"?>
<GTH_output xmlns=\"http://www.genomethreader.org/GTH_output/\" GTH_XML_version=\"$GTHXMLVERS\">
  <header xmlns=\"http://www.genomethreader.org/GTH_output/header/\">
";
 
    my $line = <$in_fh>;
    $line =~ m/(\S+)\.\s+Version of (.*?).$/;
    print $out_fh "    <source program=\"$1\" version=\"$2\" build_date=\"$2\" ";

    $line = <$in_fh>;
    $line =~ m/Date run: (.*?)$/;
    print $out_fh "run_date=\"$1\"/>\n";

    until ($line =~ m/model/) { $line = <$in_fh>; }
    $line =~ m/\((.*?)\).*?:\s(.*?)$/;
    $BSSM = "    <splice_site_parameters parameter_type=\"$1\" species=\"$2\"/>";

    # Get information on vmatch or fast search parameters
    until ($line =~ m/parameters/g) { $line = <$in_fh>; }
    $PARMS = "    <parameters>\n";

    while ($line =~ m/(\S+)\s+(\d+)[,\.]/g) {
      $PARMS .= "      <parameter name=\"$1\" value=\"$2\"/>\n";
    } # end while

    # Try to find info on substitution matrix if possible
    # Note: As of GSQ2, 21 March 2004 version, the header appears truncated.
    #       This code is updated to handle such instances.
    until (($line =~ m/matrix/)                         ||
           ($line =~ m/file/)                           ||
           ($line =~ m/^Protein similarity threshold:/) ||
           ($line =~ m/^__/)
          ) {
      $line = <$in_fh>;
    }
    # Truncated header case
    if ($line =~ m/^Prot/) {
      $line =~ m/^Protein similarity threshold:\s+(\S+)\./;
      $PARMS .= "      <parameter name=\"prot_sim_thresh\" value=\"$1\"/>\n";
    }

    # Case where we have data on the substitution matrix being used
    elsif ($line =~ m/matrix/) {
      $line =~ m/^.*?matrix:\s+(.*?)$/;
      $PARMS .= "      <parameter name=\"substitution_matrix\" value=\"$1\"/>\n";     

      # Get info on the template files used
      until ($line =~ m/file/g) { $line = <$in_fh>; }
      $line =~ m/:\s+([^,\s]+)[,]?/g;
      my $tmp = $1;
      chomp $tmp;
      print $out_fh
"    <gDNA_template_files>
      <temp_name>$tmp</temp_name>
";
      while ($line =~ m/\s+([^,\s]+)[,]?/g) {
        my $tmp = $1;
        chomp $tmp;
        print $out_fh "      <temp_name>$tmp</temp_name>\n";
      }
      print $out_fh 
"    </gDNA_template_files>
";
    }

    # This conditions signals the end of the header
    elsif ($line =~ m/^__/) {
      print $out_fh
"    <gDNA_template_files>
      <temp_name>$HOPELESS_SITUATION</temp_name>
";
      print $out_fh 
"    </gDNA_template_files>
";
    }

    else {
      # Get info on the template files used
      until ($line =~ m/file/g) { $line = <$in_fh>; }
      $line =~ m/:\s+([^,\s]+)[,]?/g;
      my $tmp = $1;
      chomp $tmp;
      print $out_fh
"    <gDNA_template_files>
      <temp_name>$tmp</temp_name>
";
      while ($line =~ m/\s+([^,\s]+)[,]?/g) {
        my $tmp = $1;
        chomp $tmp;
        print $out_fh "      <temp_name>$tmp</temp_name>\n";
      }
      print $out_fh 
"    </gDNA_template_files>
";
    }

    # End of conditionals
    print $out_fh 
"    <reference_files>
";
    # Report the reference file data
    foreach my $ref_file (@reference_files) {
      print $out_fh "      <file ref_name=\"$ref_file\" type=\"$reference_type\"/>\n";
    }
    print $out_fh "    </reference_files>
";
    print $out_fh "$BSSM\n";
    $PARMS .= "    </parameters>\n";
    print $out_fh "$PARMS";
    print $out_fh
"    <overall_reference_type>$reference_type</overall_reference_type>
";
    print $out_fh "  </header>
";
  } # end HEADER_PARSE
} # end parse_header

###########################################################################
 
# Description: Finds a block of GSQ output that corresponds to a
#              significant alignment and invokes the proper functions
#              to parse it.
sub parse_body {
  my $in_fh = shift(@_);
  my $out_fh = shift(@_);
  my $reference_type = shift(@_);

  while($TRUE) {
    my $line = <$in_fh>;
    if (!defined($line)) { last; }
    until ((!defined($line)) || ($line =~ m/^Sequence/)) { $line = <$in_fh>; }
    if (!defined($line)) { last; }

    $line =~ m/:\s+(\S+)\,/;
    my $template_id = $1; # We'll save this info for later...
    until ($line =~ m/library file:/) { $line = <$in_fh>; }
    $line =~ m/:\s+(\S+[^;])[;\s]/;
    my $reference_source_file = $1; # We'll save this info for later...

    # We now determine if we are dealing with a significant alignment or
    # not.  Note that this decision hinges on format charactersitics of
    # the GSQ plain text output, which is seemingly risky business, but
    # I think it is consistent in all such documents
    until ((!defined($line))             ||
           ($line =~ m/^No significant/) || 
           ($line =~ m/^[*]{1,}/)
          ) {
      $line = <$in_fh>;
    }
    if (!defined($line)) { last; }

    # Branch to handle the case when there is no significant 
    # information at hand
    elsif ($line =~ m/^No significant/) { next; }

    # Branch to handle the case when there is significant 
    # information at hand
    elsif ($line =~ m/^[*]/) {
      # Determine what is the correct set of parsing subroutines
      # to invoke on the data.
      if ($reference_type eq $ref_types{nucleic_acid}) {
        &parse_alignment_module_ESTcDNA($in_fh, $out_fh,
         $reference_source_file, $template_id);
        &parse_PGL_module_ESTcDNA($in_fh, $out_fh,
         $template_id);
      }
      elsif ($reference_type eq $ref_types{protein}) {
        &parse_alignment_module_Protein($in_fh, $out_fh,
         $reference_source_file, $template_id);
        &parse_PGL_module_Protein($in_fh, $out_fh,
         $template_id);
      }
      else { &fatal($SPARKS_AT_FAULT); }
    }
    else { &fatal($SPARKS_AT_FAULT); }
  } # end while

  print $out_fh
"</GTH_output>
";
} # end parse_body

###########################################################################

# Description: Parses content corresponding to the 
#              alignment_module as described in
#              "GenomeThreader.RNG"
sub parse_alignment_module_ESTcDNA {
  my $in_fh = shift(@_);
  my $out_fh = shift(@_);
  my $reference_source_file = shift(@_);
  my $template_id = shift(@_);
  my $line = "";
  my $ESTseq = "";
  my ($genome_strand,$mrna_strand) = "";

  # Declare the alignment module
  print $out_fh
"  <alignment_module>\n";

  # Frankly, the code in the block below is quite self-explanatory,
  # so I'll only add commentary where it seems necessary.
  PARSE_ALIGNMENT_MODULE_ESTCDNA: while($TRUE) {
    print $out_fh
"    <spliced_alignment xmlns=\"http://www.genomethreader.org/GTH_output/alignment_module/spliced_alignment/\">
      <reference ";

    until ($line =~ m/^EST sequence/) { $line = <$in_fh>; }
    $line =~ m/^.*?(\S+)([+-])\)$/;
    print $out_fh 
"ref_file=\"$reference_source_file\" ref_id=\"$1\" ref_strand=\"$2\" ref_description=\"$1\">\n";

    # Parse the reference's raw sequence data
    until ($line =~ m/^\s+\d+/) { $line = <$in_fh>; }
    
    $ESTseq = "";
    until ($line !~ m/\S/) {
      $line =~ s/[\s\d]//g;
      $ESTseq .= $line;
      $line = <$in_fh>;
    }
    print $out_fh "        <seq>$ESTseq</seq>\n";
    # Here is the only component of the parsing process that
    # cannot be completed using GSQ plain text output, namely identifying
    # which file the template sequence originated from.  The best I
    # can do at present is give the id of the template.
    print $out_fh
"      </reference>
      <gDNA_segment>
        <template temp_file=\"$HOPELESS_SITUATION\" ";
    # Get the template threading coordinates 
    until ($line =~ m/(\d+) to (\d+)\):$/) { $line = <$in_fh>; }
    # Determine template orientation
    $line =~ m/(\d+) to (\d+)\):$/;
    my $orientation = "";
    if ($1 < $2) { $orientation = "+"; }
    else         { $orientation = "-"; }
    print $out_fh "temp_id=\"$template_id\" temp_strand=\"$orientation\" temp_description=\"$template_id\">
          <position start=\"$1\" stop=\"$2\"/>
        </template>
      </gDNA_segment>
      <predicted_gene_structure>
        <exon-intron_info>
";
    until ($line =~ m/(Exon|Intron)/) { $line = <$in_fh>; }
    while ($line =~ m/(Exon|Intron)/) {
      if ($line =~ m/Exon/) {
        $line =~ m/(\d+)\s+(\d+)\s+(\d+).*?(\d+).*?;\s+(\S+)\s+(\d+)\s+(\d+).*?(\d+).*?;.*?(\S+)$/;
        print $out_fh
"          <exon e_serial=\"$1\">
            <gDNA_exon_boundary g_start=\"$2\" g_stop=\"$3\" g_length=\"$4\"/>
            <reference_exon_boundary r_type=\"$5\" r_start=\"$6\" r_stop=\"$7\" r_length=\"$8\" r_score=\"$9\"/>
          </exon>
";
    }
      elsif ($line =~ m/Intron/) {
        $line =~ m/(\d+)\s+(\d+)\s+(\d+).*?(\d+).*?;\s+Pd:\s+(\S+).*?(\S+)\),\s+Pa:\s+(\S+).*?(\S+)\).*?$/;
        print $out_fh
"          <intron i_serial=\"$1\">
            <gDNA_intron_boundary i_start=\"$2\" i_stop=\"$3\" i_length=\"$4\">
              <donor d_prob=\"$5\" d_score=\"$6\"/>
              <acceptor a_prob=\"$7\" a_score=\"$8\"/>
            </gDNA_intron_boundary>
          </intron>
";
      }
      else { &fatal($SPARKS_AT_FAULT); }
      $line = <$in_fh>;
    } # end while

    print $out_fh
"        </exon-intron_info>
";
    # Was a PPA predicted?
    if ($line =~ m/PPA/) {
      $line =~ m/PPA.*?(\d+)\s+(\d+)$/;
      print $out_fh "        <PPA_line polyA_start=\"$1\" polyA_stop=\"$2\"/>\n";
    }

    until ($line =~ m/^MATCH/) { $line = <$in_fh>; } # end until
    $line =~ m/\s+(.*?)([+-])\s+(.*?)([+-])\s+(\S+)\s+(\d+)\s+(\S+)\s+([CG])/;
    print $out_fh
"        <MATCH_line gen_id=\"$1\" gen_strand=\"$2\" ref_id=\"$3\" ref_strand=\"$4\">
          <total_alignment_score>$5</total_alignment_score>
          <cumulative_length_of_scored_exons>$6</cumulative_length_of_scored_exons>
          <coverage percentage=\"$7\" high_type=\"$8\"/>
        </MATCH_line>
";
    until ($line =~ m/^PGS/) { $line = <$in_fh>; } # end until
    $line =~ m/PGS_(\S+)([+-])_(\S+)([+-])\s+\(/g;
    print $out_fh
"        <PGS_line>
          <gDNA gen_id=\"$1\" gen_strand=\"$2\"/>
          <rDNA rDNA_id=\"$3\" rDNA_strand=\"$4\"/>
          <gDNA_exon_coordinates>
";
    while ($line =~ m/(\d+)\s+(\d+)/g) {
      print $out_fh
"            <exon e_start=\"$1\" e_stop=\"$2\"/>
";
    } # end while
    print $out_fh 
"          </gDNA_exon_coordinates>
        </PGS_line>
        <alignment>
";  

    until ($line =~ m/^Alignment/) { $line = <$in_fh>; } # end until
    $line = <$in_fh>;
  
    # Parse out the alignment
    ($genome_strand,$mrna_strand)="";
    while($TRUE) {
      $line = <$in_fh>;
      last if (($line !~ m/\S/)     ||
               ($line =~ m/^hqPGS/) ||
               ($line =~ m/\(.*?\)/));
      $line =~ s/^(.{1,65})//; $genome_strand .= $1;

      $line = <$in_fh>; $line = <$in_fh>;
      $line =~ s/^(.{1,65})//; $mrna_strand .= $1;

      $line = <$in_fh>; $line = <$in_fh>;
    } # end while 
    $genome_strand =~ s/[\s\d]//g;
    $mrna_strand =~ s/[\s\d]//g;

    print $out_fh
"          <genome_strand>$genome_strand</genome_strand>
            <mrna_strand>$mrna_strand</mrna_strand>
        </alignment>
";
    # Parse hqPGS, if it exists
    if ($line =~ m/^hqPGS/) {
      $line =~ m/hqPGS_(\S+)([+-])_(\S+)([+-])\s+\(/g;
      print $out_fh
"        <hqPGS_line>
          <hqgDNA gen_id=\"$1\" gen_strand=\"$2\"/>
          <hqrDNA rDNA_id=\"$3\" rDNA_strand=\"$4\"/>
          <hqgDNA_exon_coordinates>
";
    while ($line =~ m/(\d+)\s+(\d+)/g) {
      print $out_fh
"            <hqexon e_start=\"$1\" e_stop=\"$2\"/>
";
    } # end while
    print $out_fh
"          </hqgDNA_exon_coordinates>
        </hqPGS_line>
";
    }

    print $out_fh
"      </predicted_gene_structure>
    </spliced_alignment>
";
    until (($line =~ m/Total\s+number/) ||
           ($line =~ m/^[*]/)) { $line = <$in_fh>; }
    if    ($line =~ m/^[*]/)           { next PARSE_ALIGNMENT_MODULE_ESTCDNA; }
    elsif ($line =~ m/Total\s+number/) { last PARSE_ALIGNMENT_MODULE_ESTCDNA; }
    else                               { &fatal($SPARKS_AT_FAULT);            }
  } # end PARSE_ALIGNMENT_MODULE_ESTCDNA
  
  $line =~ m/:\s+(\d+)/;
  print $out_fh
"    <total_number_ESTs_reported>$1</total_number_ESTs_reported>
  </alignment_module>
";   
} # end parse_alignment_module_ESTcDNA

###########################################################################

# Description: Parses content corresponding to the 
#              PGL_module element as described in
#              "GenomeThreader.RNG"
sub parse_PGL_module_ESTcDNA {
  my $in_fh = shift(@_);
  my $out_fh = shift(@_);
  my $template_id = shift(@_);
  my $line = "";
  my ($gDNA_template,$first_frame,$second_frame,$third_frame) = "";
  my $pps = "";

  print $out_fh
"  <PGL_module xmlns=\"http://www.genomethreader.org/GTH_output/PGL_module/\">
";

  # Frankly, the code in the block below is quite self-explanatory,
  # so I'll only add commentary where it seems necessary.
  PARSE_PGL_MODULE_ESTCDNA: while($TRUE) {
    until ($line =~ m/^PGL/) {
      $line = <$in_fh>;
    }
    $line =~ m/PGL\s+(\d+).*?([+-]).*?(\d+)\s+(\d+)$/;
    my $PGL_serial_number = $1; # Use this now and also save for later
    print $out_fh
"    <predicted_gene_location>
      <PGL_line PGL_serial=\"$PGL_serial_number\" PGL_strand=\"$2\" PGL_start=\"$3\" PGL_stop=\"$4\"/>
";
    PARSE_AGS_SECTION: while($TRUE) {
      print $out_fh
"      <AGS_information>
        <AGS_line ";

      until ($line =~ m/^AGS/) { $line = <$in_fh>; }
      $line =~ m/AGS-(\d+)/g;
      print $out_fh "AGS_serial=\"$1\">\n";

      print $out_fh "          <exon_coordinates>\n";
      while ($line =~ m/(\d+)\s+(\d+)/g) {
        print $out_fh "            <exon e_start=\"$1\" e_stop=\"$2\"/>\n";
      } # end while
      print $out_fh
"          </exon_coordinates>
        </AGS_line>
        <SCR_line>
";

      until ($line =~ m/^SCR/) { $line = <$in_fh>; }
      if ($line =~ m/e.*?d.*?a/) {
        while ($line =~ m/e\s+(\S+)\s+d\s+(\S+)\s+a\s+(\S+),/g) {
          print $out_fh
"          <exon-intron don_prob=\"$2\" acc_prob=\"$3\" e_score=\"$1\"/>
";
        } # end while
      }
      # Parse out the last "exon token"
      $line =~ m/e\s+(\S+)\)$/;
      print $out_fh
"          <exon-only e_score=\"$1\"/>
        </SCR_line>
        <exon-intron_info xmlns=\"http://www.genomethreader.org/GTH_output/PGL_module/predicted_gene_location/AGS_information/exon-intron_info/\">
";

      until ($line =~ m/(Exon|Intron)/) { $line = <$in_fh>; }
      while ($line =~ m/(Exon|Intron)/) {
        if ($line =~ m/Exon/) {
          $line =~ m/(\d+)\s+(\d+)\s+(\d+).*?(\d+).*?;.*?:\s+(\S+)/;
          print $out_fh
"          <exon e_serial=\"$1\" e_score=\"$5\">
            <gDNA_exon_boundary e_start=\"$2\" e_stop=\"$3\" e_length=\"$4\"/>
          </exon>
";
        }
        elsif ($line =~ m/Intron/) {
          $line =~ m/(\d+)\s+(\d+)\s+(\d+).*?(\d+).*?;\s+Pd:\s+(\S+)\s+Pa:\s+(\S+)/;
          print $out_fh
"          <intron i_serial=\"$1\" don_prob=\"$5\" acc_prob=\"$6\">
            <gDNA_intron_boundary i_start=\"$2\" i_stop=\"$3\" i_length=\"$4\"/>
          </intron>
";
        }
        else { &fatal($SPARKS_AT_FAULT); }
        $line = <$in_fh>;
      } # end while

      print $out_fh
"        </exon-intron_info>
        <supporting_evidence xmlns=\"http://www.genomethreader.org/GTH_output/PGL_module/predicted_gene_location/AGS_information/supporting_evidence/\">
";
      until ($line =~ m/PGS/) { $line = <$in_fh>; }
      # Expect one-to-many PGS lines
      while ($line =~ m/PGS/) {
        print $out_fh
"          <PGS_line>
            <gDNA_exon_coordinates>
";
        $line =~ m/PGS.*?(\S+)([+-])$/;
        my ($id, $orientation) = ($1, $2);
        while ($line =~ m/(\d+)\s+(\d+)/g) {
          print $out_fh
"              <exon start=\"$1\" stop=\"$2\"/>
";
        } # end while 
        print $out_fh
"            </gDNA_exon_coordinates>
            <referenceDNA id=\"$id\" strand=\"$orientation\"/>
          </PGS_line>
";
        $line = <$in_fh>;
      } # end while
      print $out_fh
"        </supporting_evidence>
";
      # The following while-loop parses the 1 or 2 3-phase translation
      # sections of the output.
      PARSE_THREE_PHASE_TRANSLATION: while($TRUE) {
        until ((!defined($line))      ||
               ($line =~ m/^3-phase/) ||
               ($line =~ m/^AGS/)     ||
               ($line =~ m/^PGL/)     ||
               ($line =~ m/^[_]{1,}/)
              ) {
          $line = <$in_fh>;
        }
        # Case where there isn't more data to parse
        if ((!defined($line)) || ($line =~ m/^[_]{1,}/)) {
          print $out_fh
"      </AGS_information>
    </predicted_gene_location>
  </PGL_module>
";
          last PARSE_PGL_MODULE_ESTCDNA;
        }
        elsif ($line =~ m/^AGS/) {
          print $out_fh
"      </AGS_information>
";
          last PARSE_THREE_PHASE_TRANSLATION;
        }
        elsif ($line =~ m/^PGL/) {
          print $out_fh
"      </AGS_information>
    </predicted_gene_location>
";
          last PARSE_AGS_SECTION;
        }
        elsif ($line =~ m/^3-phase/) {
          $line =~ m/AGS-(\d+)\s+\(([+-])strand\):/;
          print $out_fh
"        <three_phase_translation xmlns=\"http://www.genomethreader.org/GTH_output/PGL_module/predicted_gene_location/AGS_information/three_phase_translation/\">
          <description PGL_serial=\"$PGL_serial_number\" AGS_serial=\"$1\" gDNA_strand=\"$2\"/>
          <translation>
";
          # Parse the translation
          ($gDNA_template,$first_frame,$second_frame,$third_frame)="";
          $line = <$in_fh>;
          while($TRUE) {
            $line = <$in_fh>;
            # Allows an exit when finished parsing the translation
            last if ($line =~ m/Maximal/);

            $line = <$in_fh>; chomp($line);
            $line =~ s/^(.{12})//; $gDNA_template .= $line;

            $line = <$in_fh>; chomp($line);
            $line =~ s/^(.{12})//; $first_frame .= $line;

            $line = <$in_fh>; chomp($line);
            $line =~ s/^(.{12})//; $line =~ s/(\s)$//; $second_frame .= $line;

            $line = <$in_fh>; chomp($line);
            $line =~ s/^(.{12})//; $third_frame .= "$line ";

            $line = <$in_fh>; $line = <$in_fh>;
          } # end while

          my $translen = length($gDNA_template);
          my @gDNA_templateTMP = split(//, $gDNA_template);
          $first_frame .= "      ";
          my @first_frameTMP   = split(//, $first_frame);
          $second_frame .= "      ";
          my @second_frameTMP  = split(//, $second_frame);
          $third_frame .= "      ";
          my @third_frameTMP  = split(//, $third_frame);
          ($gDNA_template,$first_frame,$second_frame,$third_frame)="";
          for (my $i=0;$i<$translen;++$i) {
            $gDNA_template .= $gDNA_templateTMP[$i];  
            $first_frame   .= $first_frameTMP[$i];
            $second_frame  .= $second_frameTMP[$i];
            $third_frame   .= $third_frameTMP[$i];
          }

          print $out_fh
"            <gDNA_template>$gDNA_template</gDNA_template>
              <first_frame>$first_frame</first_frame>
             <second_frame>$second_frame</second_frame>
              <third_frame>$third_frame</third_frame>
          </translation>
          <probable_ORFs xmlns=\"http://www.genomethreader.org/GTH_output/PGL_module/predicted_gene_location/AGS_information/three_phase_translation/probable_ORFs/\">
";
          # Start parsing the probable ORFs data
          until (($line =~ m/^>/) || ($line =~ m/none/)) { $line = <$in_fh>; }
          if ($line =~ m/^>/) {
            # Get info on maximal non-overlapping open reading frames
            while ($line =~ m/^>/) { 
              $line =~ m/>(.*?)([+-])_PGL-(\d+)_AGS-(\d+)_PPS_(\d+)/;
              print $out_fh
"            <orf_entry>
              <id_line>
                <gDNA id=\"$1\" strand=\"$2\"/>
                <serials PGL_serial=\"$3\" AGS_serial=\"$4\" PPS_serial=\"$5\"/>
                <orf_info>
                  <exon_boundaries>
";
              while ($line =~ m/(\d+)\s+(\d+)/g) {
                print $out_fh
"                    <exon start=\"$1\" stop=\"$2\"/>
";
              } # end while
              $line =~ m/'(\d)';\s+(\d+)\s+bp,\s+(\d+).*?$/;
              print $out_fh
"                  </exon_boundaries>
                  <frame>$1</frame>
                  <number_coding_nucleotides>$2</number_coding_nucleotides>
                  <number_encoded_amino_acids>$3</number_encoded_amino_acids>
                </orf_info>
              </id_line>
              <predicted_protein_sequence>";
              $pps="";
              while($TRUE) {
                $line = <$in_fh>;
                # Handles case where we are finished parsing rows of the ORF
                if ($line !~ m/\S/) { last; }
                # Parses a row of the ORF
                else {
                  chomp($line);
                  $line =~ s/[\d\s]//g;
                  $pps .= $line;
                }
              } # end while

              print $out_fh "$pps</predicted_protein_sequence>
            </orf_entry>
";
              $line = <$in_fh>;
            } # end while
          } # end if

          # Handles case where there aren't any maximal ORFs to consider
          elsif ($line =~ m/none/) {
            print $out_fh
"            <none gDNA_id=\"$template_id\"/>
";
          } # end elsif
          else { &fatal($SPARKS_AT_FAULT); }

          print $out_fh
"          </probable_ORFs>
        </three_phase_translation>
";
        } # end elsif
        else { &fatal($SPARKS_AT_FAULT); }
      } # end PARSE_THREE_PHASE_TRANSLATION
    } # end PARSE_AGS_SECTION
  } # end PARSE_PGL_MODULE_ESTCDNA
} # end parse_PGL_module_ESTcDNA

###########################################################################

# Description: Parses content corresponding to the 
#              alignment_module element as described in
#              "GenomeThreader.RNG"
sub parse_alignment_module_Protein {
  my $in_fh = shift(@_);
  my $out_fh = shift(@_);
  my $reference_source_file = shift(@_);
  my $template_id = shift(@_);
  my $line = "";
  my $PROTseq = "";
  my ($genomeDNA,$genomeProt,$queryProt) = "";

  # Declare the alignment module
  print $out_fh
"  <alignment_module>
";

  # Frankly, the code in the block below is quite self-explanatory,
  # so I'll only add commentary where it seems necessary.
  PARSE_ALIGNMENT_MODULE_PROTEIN: while($TRUE) {

    print $out_fh
"    <spliced_alignment xmlns=\"http://www.genomethreader.org/GTH_output/alignment_module/spliced_alignment/\">
";
    until ($line =~ m/protein sequence/g) { $line = <$in_fh>; }
    $line =~ m/\(File: (.*?)\)$/;

    print $out_fh
"      <reference ref_file=\"$reference_source_file\" ref_id=\"$1\" ref_description=\"$1\">
        <seq>";
    # Parse the reference's raw sequence data
    until ($line =~ m/^\s+\d+/) { $line = <$in_fh>; }
    $PROTseq = "";
    until ($line !~ m/\S/) {
      chomp($line);
      $line =~ s/[\s\d]//g;
      $PROTseq .= $line;
      $line = <$in_fh>;
    }
    print $out_fh "$PROTseq</seq>
      </reference>
      <gDNA_segment>
        <template temp_file=\"$HOPELESS_SITUATION\" temp_id=\"$template_id\" temp_description=\"$template_id\" ";

    # Get the template threading coordinates 
    until ($line =~ m/(\d+) to (\d+)\):$/) { $line = <$in_fh>; }
    # Determine template orientation
    $line =~ m/(\d+) to (\d+)\):$/;
    my $orientation = "";
    if ($1 < $2) { $orientation = "+"; }
    else         { $orientation = "-"; }

    print $out_fh "temp_strand=\"$orientation\">
          <position start=\"$1\" stop=\"$2\"/>
        </template>
      </gDNA_segment>
      <predicted_gene_structure>
        <exon-intron_info>
";
    until ($line =~ m/(Exon|Intron)/) { $line = <$in_fh>; }
    while ($line =~ m/(Exon|Intron)/) {
      if ($line =~ m/Exon/) {
        $line =~ m/(\d+)\s+(\d+)\s+(\d+).*?(\d+).*?;\s+(\S+)\s+(\d+)\s+(\d+).*?(\d+).*?;.*?(\S+)$/;
        print $out_fh
"          <exon e_serial=\"$1\">
            <gDNA_exon_boundary g_start=\"$2\" g_stop=\"$3\" g_length=\"$4\"/>
            <reference_exon_boundary r_type=\"$5\" r_start=\"$6\" r_stop=\"$7\" r_length=\"$8\" r_score=\"$9\"/>
          </exon>
";
      }
      elsif ($line =~ m/Intron/) {
        $line =~ m/(\d+)\s+(\d+)\s+(\d+).*?(\d+).*?;\s+Pd:\s+(\S+)\s+Pa:\s+(\S+)$/;
        print $out_fh
"          <intron i_serial=\"$1\">
            <gDNA_intron_boundary i_start=\"$2\" i_stop=\"$3\" i_length=\"$4\">
              <donor d_prob=\"$5\"/>
              <acceptor a_prob=\"$6\"/>
            </gDNA_intron_boundary>
          </intron>
";
      }
      else { &fatal($SPARKS_AT_FAULT); }
      $line = <$in_fh>;
    } # end while
    print $out_fh
"        </exon-intron_info>
";
    until ($line =~ m/^MATCH/) { $line = <$in_fh>; }
    $line =~ m/\s+(.*?)([+-])\s+(.*?)\s+(\S+)\s+(\d+)\s+(\S+)\s+([PG])/;
    print $out_fh
"        <MATCH_line gen_id=\"$1\" gen_strand=\"$2\" ref_id=\"$3\">
          <total_alignment_score>$4</total_alignment_score>
          <cumulative_length_of_scored_exons>$5</cumulative_length_of_scored_exons>
          <coverage percentage=\"$6\" high_type=\"$7\"/>
        </MATCH_line>
";
    until ($line =~ m/^PGS/) { $line = <$in_fh>; }
    $line =~ m/PGS_(\S+)([+-])_(\S+)\s+\(/g;
    print $out_fh
"        <PGS_line>
            <gDNA gen_id=\"$1\" gen_strand=\"$2\"/>
            <rProt rProt_id=\"$3\"/>
          <gDNA_exon_coordinates>
";
    while ($line =~ m/(\d+)\s+(\d+)/g) {
      print $out_fh
"            <exon e_start=\"$1\" e_stop=\"$2\"/>
";
    } # end while
    print $out_fh 
"          </gDNA_exon_coordinates>
        </PGS_line>
        <alignment>
";  

    until ($line =~ m/^Alignment/) {
      $line = <$in_fh>;
    } # end until
    $line = <$in_fh>;

    # Parse out the alignment
    ($genomeDNA,$genomeProt,$queryProt)="";
    while($TRUE) {
      $line = <$in_fh>;
      last if ($line !~ m/\S/);
      chomp($line);
      $line =~ s/^(.{1,65})//;
      $genomeDNA .= $1;

      $line = <$in_fh>; chomp($line);
      $line =~ s/^(.{1,65})//;
      $genomeProt .= $1;

      $line = <$in_fh>; $line = <$in_fh>; chomp($line);
      $line =~ s/^(.{1,65})//;
      $queryProt .= $1;

      $line = <$in_fh>; $line = <$in_fh>;
    } # end while 
    $genomeDNA =~ s/(.{13})$//;
    $queryProt =~ s/(.{13})$//;

    my $translen=length($genomeDNA);
    my @genomeDNATMP  = split(//, $genomeDNA);
    my @genomeProtTMP = split(//, $genomeProt);
    my @queryProtTMP  = split(//, $queryProt);
    ($genomeDNA,$genomeProt,$queryProt)="";
    for (my $i=0;$i<$translen;++$i) {
      if ($genomeDNATMP[$i] eq ' ') { next; }
      else {
        $genomeDNA  .= $genomeDNATMP[$i];  
        $genomeProt .= $genomeProtTMP[$i];
        $queryProt  .= $queryProtTMP[$i];
      }
    }

    print $out_fh
"          <genome_strand>$genomeDNA</genome_strand>
         <genomeProt>$genomeProt</genomeProt>
          <queryProt>$queryProt</queryProt>
        </alignment>
      </predicted_gene_structure>
    </spliced_alignment>
";
    until (($line =~ m/Total\s+number/) || ($line =~ m/^[*]/)) {
      $line = <$in_fh>;
    }
    if ($line =~ m/^[*]/) {
      next PARSE_ALIGNMENT_MODULE_PROTEIN;
    }
    elsif ($line =~ m/Total\s+number/) {
      last PARSE_ALIGNMENT_MODULE_PROTEIN;
    }
    else {
      &fatal($SPARKS_AT_FAULT);
    }
  } # end PARSE_ALIGNMENT_MODULE_PROTEIN

  $line =~ m/:\s+(\d+)/;
  print $out_fh
"    <total_number_proteins_reported>$1</total_number_proteins_reported>
  </alignment_module>
";   
} # end parse_alignment_module_Protein

###########################################################################

# Description: Parses content corresponding to the 
#              PGL_module_Protein element as described in
#              "GenomeThreader.RNG"
sub parse_PGL_module_Protein {
  my $in_fh = shift(@_);
  my $out_fh = shift(@_);
  my $template_id = shift(@_);
  my $line = "";
  my ($gDNA_template,$first_frame,$second_frame,$third_frame) = "";
  my $pps = "";

  print $out_fh
"  <PGL_module xmlns=\"http://www.genomethreader.org/GTH_output/PGL_module/\">
";
  # Frankly, the code in the block below is quite self-explanatory,
  # so I'll only add commentary where it seems necessary.
  PARSE_PGL_MODULE_PROTEIN: while($TRUE) {
    until ($line =~ m/^PGL/) { $line = <$in_fh>; }
    $line =~ m/PGL\s+(\d+).*?([+-]).*?(\d+)\s+(\d+)$/;
    my $PGL_serial_number = $1; # Use this now and also save for later
    print $out_fh
"    <predicted_gene_location>
      <PGL_line PGL_serial=\"$PGL_serial_number\" PGL_strand=\"$2\" PGL_start=\"$3\" PGL_stop=\"$4\"/>
";
    PARSE_AGS_SECTION: while($TRUE) {
      print $out_fh
"      <AGS_information>
        <AGS_line ";

      until ($line =~ m/^AGS/) { $line = <$in_fh>; }
      $line =~ m/AGS-(\d+)/g;
      print $out_fh "AGS_serial=\"$1\">
          <exon_coordinates>
";
      while ($line =~ m/(\d+)\s+(\d+)/g) {
        print $out_fh
"            <exon e_start=\"$1\" e_stop=\"$2\"/>
";    
      } # end while
      print $out_fh
"          </exon_coordinates>
        </AGS_line>
        <SCR_line>
";

      until ($line =~ m/^SCR/) { $line = <$in_fh>; }
      if ($line =~ m/e.*?d.*?a/) {
        while ($line =~ m/e\s+(\S+)\s+d\s+(\S+)\s+a\s+(\S+),/g) {
          print $out_fh
"          <exon-intron don_prob=\"$2\" acc_prob=\"$3\" e_score=\"$1\"/>
";
        } # end while
      }
      # Parse out the last "exon token"
      $line =~ m/e\s+(\S+)\)$/;
      print $out_fh
"          <exon-only e_score=\"$1\"/>
        </SCR_line>
        <exon-intron_info xmlns=\"http://www.genomethreader.org/GTH_output/PGL_module/predicted_gene_location/AGS_information/exon-intron_info/\">
";

      until ($line =~ m/(Exon|Intron)/) { $line = <$in_fh>; }
      while ($line =~ m/(Exon|Intron)/) {
        if ($line =~ m/Exon/) {
          $line =~ m/(\d+)\s+(\d+)\s+(\d+).*?(\d+).*?;.*?:\s+(\S+)/;
          print $out_fh
"          <exon e_serial=\"$1\" e_score=\"$5\">
            <gDNA_exon_boundary e_start=\"$2\" e_stop=\"$3\" e_length=\"$4\"/>
          </exon>
";
        }
        elsif ($line =~ m/Intron/) {
          $line =~ m/(\d+)\s+(\d+)\s+(\d+).*?(\d+).*?;\s+Pd:\s+(\S+)\s+Pa:\s+(\S+)/;
          print $out_fh
"          <intron i_serial=\"$1\" don_prob=\"$5\" acc_prob=\"$6\">
            <gDNA_intron_boundary i_start=\"$2\" i_stop=\"$3\" i_length=\"$4\"/>
          </intron>
";
        }
        else { &fatal($SPARKS_AT_FAULT); }
        $line = <$in_fh>;
      } # end while
      print $out_fh
"        </exon-intron_info>
        <supporting_evidence xmlns=\"http://www.genomethreader.org/GTH_output/PGL_module/predicted_gene_location/AGS_information/supporting_evidence/\">
";
      until ($line =~ m/PGS/) { $line = <$in_fh>; }
      # Expect one-to-many PGS lines
      while ($line =~ m/PGS/) {
        print $out_fh
"          <PGS_line>
            <gDNA_exon_coordinates>
";
        $line =~ m/PGS.*?\)\s+(\S+).*?$/;
        my $id = $1;
        while ($line =~ m/(\d+)\s+(\d+)/g) {
          print $out_fh
"              <exon start=\"$1\" stop=\"$2\"/>
";
        } # end while 
        print $out_fh
"            </gDNA_exon_coordinates>
            <referenceProtein id=\"$id\"/>
          </PGS_line>
";
        $line = <$in_fh>;
      } # end while
      print $out_fh
"        </supporting_evidence>
";
      # The following while-loop parses the 1 or 2 3-phase translation
      # sections of the output.
      PARSE_THREE_PHASE_TRANSLATION: while($TRUE) {
        until ((!defined($line))      ||
               ($line =~ m/^3-phase/) ||
               ($line =~ m/^AGS/)     ||
               ($line =~ m/^PGL/)     ||
               ($line =~ m/^[_]{1,}/)
              ) {
          $line = <$in_fh>;
        }
        # Case where there isn't more data to parse
        if ((!defined($line)) || ($line =~ m/^[_]{1,}/)) {
          print $out_fh
"      </AGS_information>
    </predicted_gene_location>
  </PGL_module>
";
          last PARSE_PGL_MODULE_PROTEIN;
        }
        elsif ($line =~ m/^AGS/) {
          print $out_fh
"      </AGS_information>
";
          last PARSE_THREE_PHASE_TRANSLATION;
        }
        elsif ($line =~ m/^PGL/) {
          print $out_fh
"      </AGS_information>
    </predicted_gene_location>
";
          last PARSE_AGS_SECTION;
        }
        elsif ($line =~ m/^3-phase/) {
          $line =~ m/AGS-(\d+)\s+\(([+-])strand\):/;
          print $out_fh
"        <three_phase_translation xmlns=\"http://www.genomethreader.org/GTH_output/PGL_module/predicted_gene_location/AGS_information/three_phase_translation/\">
          <description PGL_serial=\"$PGL_serial_number\" AGS_serial=\"$1\" gDNA_strand=\"$2\"/>
          <translation>
";
          # Parse the translation
          ($gDNA_template,$first_frame,$second_frame,$third_frame)="";
          $line = <$in_fh>;
          while($TRUE) {
            $line = <$in_fh>;
            # Allows an exit when finished parsing the translation
            last if ($line =~ m/Maximal/);

            $line = <$in_fh>; chomp($line);
            $line =~ s/^.{10}(.*?)$//g;
            $gDNA_template .= $1;

            $line = <$in_fh>; chomp($line);
            $line =~ s/^.{10}(.*?)$//g;
            $first_frame .= $1;

            $line = <$in_fh>; chomp($line);
            $line =~ s/^.{10}(.*?)$//g;
            $second_frame .= $1;
            chop($second_frame);

            $line = <$in_fh>; chomp($line);
            $line =~ s/^.{10}(.*?)$//g;
            $third_frame .= "$1 ";

            $line = <$in_fh>; $line = <$in_fh>;
          } # end while

          my $translen=length($gDNA_template);

          my @gDNA_templateTMP = split(//, $gDNA_template);
          # Buffer 3' end of translations with a little whitespace (harmless)
          $first_frame        .= "      ";
          $second_frame       .= "      ";
          $third_frame        .= "      ";
          my @first_frameTMP   = split(//, $first_frame);
          my @second_frameTMP  = split(//, $second_frame);
          my @third_frameTMP   = split(//, $third_frame);
          ($gDNA_template,$first_frame,$second_frame,$third_frame)="";
          for (my $i=0;$i<$translen;++$i) {
            if ($gDNA_templateTMP[$i] eq ' ') { next; }
            else {
              $gDNA_template .= $gDNA_templateTMP[$i];  
              $first_frame   .= $first_frameTMP[$i];
              $second_frame  .= $second_frameTMP[$i];
              $third_frame   .= $third_frameTMP[$i];
            }
          }

          print $out_fh
"            <gDNA_template>$gDNA_template</gDNA_template>
              <first_frame>$first_frame</first_frame>
             <second_frame>$second_frame</second_frame>
              <third_frame>$third_frame</third_frame>
          </translation>
          <probable_ORFs xmlns=\"http://www.genomethreader.org/GTH_output/PGL_module/predicted_gene_location/AGS_information/three_phase_translation/probable_ORFs/\">
";
          # Start parsing the probable ORFs data
          until (($line =~ m/^>/) || ($line =~ m/none/)) { $line = <$in_fh>; }
          if ($line =~ m/^>/) {
            # Get info on maximal non-overlapping open reading frames
            while ($line =~ m/^>/) { 
              $line =~ m/>(.*?)([+-])_PGL-(\d+)_AGS-(\d+)_PPS_(\d+)/;
              print $out_fh
"            <orf_entry>
              <id_line>
                <gDNA id=\"$1\" strand=\"$2\"/>
                <serials PGL_serial=\"$3\" AGS_serial=\"$4\" PPS_serial=\"$5\"/>
                <orf_info>
                  <exon_boundaries>
";
              while ($line =~ m/(\d+)\s+(\d+)/g) {
                print $out_fh
"                    <exon start=\"$1\" stop=\"$2\"/>
";
              } # end while
              $line =~ m/'(\d)';\s+(\d+)\s+bp,\s+(\d+).*?$/;
              print $out_fh
"                  </exon_boundaries>
                  <frame>$1</frame>
                  <number_coding_nucleotides>$2</number_coding_nucleotides>
                  <number_encoded_amino_acids>$3</number_encoded_amino_acids>
                </orf_info>
              </id_line>
              <predicted_protein_sequence>";

              $pps="";
              while($TRUE) {
                $line = <$in_fh>;
                # Handles case where we are finished parsing rows of the ORF
                if ($line !~ m/\S/) { last; }
                # Parses a row of the ORF
                else {
                  chomp($line);
                  $line =~ s/[\d\s]//g;
                  $pps .= $line;
                }
              } # end while

              print $out_fh "$pps</predicted_protein_sequence>
            </orf_entry>
";
              $line = <$in_fh>;
            } # end while
          } # end if

          # Handles case where there aren't any maximal ORFs to consider
          elsif ($line =~ m/none/) {
            print $out_fh
"            <none gDNA_id=\"$template_id\"/>
";
          } # end elsif
          else { &fatal($SPARKS_AT_FAULT); }

          print $out_fh
"          </probable_ORFs>
        </three_phase_translation>
";
        } # end elsif
        else { &fatal($SPARKS_AT_FAULT); }
      } # end PARSE_THREE_PHASE_TRANSLATION
    } # end PARSE_AGS_SECTION
  } # end PARSE_PGL_MODULE_PROTEIN
} # end parse_PGL_module_PROTEIN

###########################################################################
