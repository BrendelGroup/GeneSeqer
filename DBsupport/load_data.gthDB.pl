#!/usr/bin/perl -w
# Michael E Sparks (mespar1@iastate.edu)
#
# Last Modified: 27 July 2006
#   (Script was updated to support gthXML v1.1)
# 495 lines, 20094 characters

###########################################################################
# DESCRIPTION:                                                            #
# This script is used to take, from STDIN, [GeneSeqer|GenomeThreader]     #
# output in XML format and store it in a MySQL database, described in     #
# ``gthDB_make.sql".  Please note that this script will only handle the   #
# XML input corresponding to one run of the [GSQ|GTH] program--handling   #
# input corresponding to multiple runs will require multiple invocations  #
# of this script (consider writing BASH shell scripts to automate).       #
#                                                                         #
# Details for the construction of the gthDB database that this script     #
# uses are described in the file "gthDB_make.sql".                        #
#                                                                         #
# IMPORTANT!  The user will need to edit the variables $USER, $PASS, and  #
# $HOST below to allow this script to engage with his/her gthDB           #
# installation.                                                           #
###########################################################################

# Copyright (c) 2003,2004,2005,2006 Michael E Sparks
# All rights reserved.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions
# are met:
#
#     * Redistributions of source code must retain the above copyright
#       notice, this list of conditions and the following disclaimer.
#     * Redistributions in binary form must reproduce the above copyright
#       notice, this list of conditions and the following disclaimer in
#       the documentation and/or other materials provided with the
#       distribution.
#     * Neither the name of the Iowa State University nor the names of
#       its contributors may be used to endorse or promote products
#       derived from this software without specific prior written permission.
#
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
# ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
# LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
# CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
# SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
# INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
# CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
# ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
# POSSIBILITY OF SUCH DAMAGE.

# Import packages
use DBI;
use strict;

# Does user wish to echo SQL load statements?
my $USAGE = "
  Usage: cat data.gthxml | script.pl SQLverbose
    where SQLverbose = T -> echo SQL load statements
                       F -> load tables silently

";

my $SQLverbose = shift or die $USAGE;
if    (($SQLverbose eq 'T') || ($SQLverbose eq 't')) { $SQLverbose=1; }
elsif (($SQLverbose eq 'F') || ($SQLverbose eq 'f')) { $SQLverbose=0; }
else  { die $USAGE; }

# User-modifiable parameters...YOU SHOULD MODIFY THESE!!!!!
my $USER = "gthuser";        # MySQL user name
my $PASS = "gthpass";        # User's password
my $HOST = "localhost";      # Host serving the gthDB database

# Non-modifiable parameters
my ($EXIT_ON_FAILURE,$EXIT_ON_SUCCESS) = (1,0);
my $DNE_flag = -100; # don_score & acc_score do not exist for protein ref's
my $HOPELESS_SITUATION = "Data not available";

# Program variables
my $dsn = "DBI:mysql:gthDB:" . "${HOST}"; # Data Source Name
my ($dbh,$sth);                           # Database and statement handles
my $query;                                # Store database query strings
my $line = "";                            # Store input from STDIN stream
my @ary = ();                             # Array to store returned query rows

# Variables for storing data items --------
# Table linkers
my ($big_run,$big_alignID,$big_pglID,$big_agsID) = "";
# header table values
my ($prog,$version,$daterun,$spliceparm,$parms,
    $temp_files,$ref_files,$ref_type) = "";
# alignments table values
my ($ref_id,$ref_strand,$ref_file,$temp_ID,$temp_strand,$tempfile,$temp_start,
    $temp_stop,$tot_score,$length_exons,$tot_exons,$coverage,$type_high) = "";
# align_exons table values
my ($exonserial,$genstart,$genstop,$refstart,$refstop,$refscore) = "";
# align_introns table values
my ($intronserial,$intstart,$intstop,$donprob,
    $donscore,$accprob,$accscore) = "";
# predgenloc table values
my ($pgl_start,$pgl_stop) = "";
# altgenstruct table values
my ($support,$tot_ags_exons) = "";
# ags_exons table values
my ($ags_exon_serial,$estart,$estop,$escore) = "";
# ags_introns table values
my ($ags_intron_serial,$istart,$istop,$donor_prob,$acceptor_prob) = "";
# prob_orfs table values
my ($pps_temp_id,$pps_temp_strand,$e_coors,$ppsframe,$peplength,$ppseq) = "";
# store truth value of high quality PGS presence
my $hqPGS = 0;

# Connect to database
$dbh = DBI->connect($dsn,$USER,$PASS,{RaiseError => $EXIT_ON_FAILURE});

# Parse in the header information
HEADER_PARSE: {
  until ($line =~ m/<source /) { $line = <>; }
  $line =~ m/<source program=\"(.*?)\" version=\"(.*?)\".*?run_date=\"(.*?)\"\/>/;
  ($prog,$version,$daterun) = ($1,$2,$3);

  until (($line =~ m/<gDNA_template_files>/) ||
         ($line =~m/<reference_files/)) { $line = <>; }
  if ($line =~ m/<gDNA_template_files>/) {
    for($line=<>;($line !~ m/<\/gDNA_template_files>/);$line=<>) {
      $line =~ m/<temp_name>(.*?)<\/temp_name>/; $temp_files .= "$1 ";
    }
  }
  else { $temp_files = $HOPELESS_SITUATION; }

  until ($line =~ m/<reference_files>/) { $line = <>; }
  for($line=<>;($line !~ m/<\/reference_files>/);$line=<>) {
    $line =~ m/ ref_name=\"(.*?)\"/; $ref_files .= "$1 ";
  }

  until ($line =~ m/<splice_site_parameters /) { $line = <>; }
  $line =~ m/ parameter_type=\"(.*?)\"/; $spliceparm = $1;
  $line =~ m/ species=\"(.*?)\"/;        $spliceparm .= " $1";

  until ($line =~ m/<parameters>/) { $line = <>; }
  for($line=<>;($line !~ m/<\/parameters/);$line=<>) {
    $line =~ m/ name=\"(.*?)\"/;  $parms .= "($1";
    $line =~ m/ value=\"(.*?)\"/; $parms .= "=$1)";
  }

  until ($line =~ m/<overall_reference_type>/) { $line = <>; }
  $line =~ m/<overall_reference_type>(.*?)<\/overall_reference_type>/; $ref_type = $1;

  $query = "INSERT INTO header(run,program,version,daterun,splicemodel,parms,template_files,ref_files,ref_type) VALUES(NULL,'$prog','$version','$daterun','$spliceparm','$parms','$temp_files','$ref_files','$ref_type')";
  if($SQLverbose) { print $query,"\n"; }
  $dbh->do($query);
} # end HEADER_PARSE

# Get $big_run, used to tie together the tables
$sth = $dbh->prepare("SELECT LAST_INSERT_ID() FROM header");
$sth->execute();
@ary = $sth->fetchrow_array();
$big_run = shift(@ary);
$sth->finish();

MAIN: while(1) {
  until (($line =~ m/<alignment_module/) ||
         ($line =~ m/<\/GTH_output>/)
        ) { $line = <>; }

  if ($line =~ m/<\/GTH_output>/) { # end of input
    # Disconnect from database, report success, and return to operating system
    $dbh->disconnect();
    my $date = `date`;
    chomp($date);
    print STDERR "Finished loading file at $date.\n";
    exit($EXIT_ON_SUCCESS);
  }

  # else -> Parse [ESTcDNA|Protein]_information
  PARSE_ALIGNMENT_MODULE: while(1) {
    until (($line =~ m/<spliced_alignment/) ||
           ($line =~ m/<\/alignment_module/)
          ) { $line = <>; }

    if ($line =~ m/<\/alignment_module/) { last PARSE_ALIGNMENT_MODULE; }

    until ($line =~ m/<reference /) { $line = <>; }

    if ($line =~ m/ref_strand/) {
      $line =~ m/ ref_file=\"(.*?)\" ref_id=\"(.*?)\" ref_strand=\"(.*?)\"/;
      ($ref_file,$ref_id,$ref_strand) = ($1,$2,$3);
    }
    else {
      $line =~ m/ ref_file=\"(.*?)\" ref_id=\"(.*?)\"/;
      ($ref_file,$ref_id,$ref_strand) = ($1,$2,"+");
    }

    until ($line =~ m/<template /) { $line = <>; }
    $line =~ m/ temp_file=\"(.*?)\" temp_id=\"(.*?)\".*?temp_strand=\"(.*?)\"/;
    ($tempfile,$temp_ID,$temp_strand) = ($1,$2,$3);

    until ($line =~ m/<position /) { $line = <>; }
    $line =~ m/ start=\"(.*?)\" stop=\"(.*?)\"/;
    ($temp_start,$temp_stop) = ($1,$2);

    # insert what values we can into alignments table
    $query = "INSERT INTO alignments(run,alignID,refid,refstrand,ref_file,templateID,templatestrand,template_file,temp_start,temp_stop) VALUES('$big_run',NULL,'$ref_id','$ref_strand','$ref_file','$temp_ID','$temp_strand','$tempfile','$temp_start','$temp_stop')";
    if($SQLverbose) { print $query,"\n"; }
    $dbh->do($query);

    # Get $big_alignID, used to tie together the tables
    $sth = $dbh->prepare("SELECT LAST_INSERT_ID() FROM alignments");
    $sth->execute();
    @ary = $sth->fetchrow_array();
    $big_alignID = shift(@ary);
    $sth->finish();

    # Fill align_exons and align_introns table
    until ($line =~ m/<exon-intron_info/) { $line = <>; }
    while(1) {
      # Fill align_exons table
      until (($line =~ m/<exon /) ||
             ($line =~ m/<\/exon-intron_info/)
            ) { $line = <>; }

      if ($line =~ m/<exon /) {
        $line =~ m/ e_serial=\"(.*?)\"/; $tot_exons = $exonserial = $1;

        until ($line =~ m/<gDNA_exon_boundary /) { $line = <>; }
        $line =~ m/ g_start=\"(.*?)\" g_stop=\"(.*?)\"/;
        ($genstart,$genstop) = ($1,$2);

        until ($line =~ m/<reference_exon_boundary /) { $line = <>; }
        $line =~ m/ r_start=\"(.*?)\" r_stop=\"(.*?)\".*?r_score=\"(.*?)\"/;
        ($refstart,$refstop,$refscore) = ($1,$2,$3);

        $query = "INSERT INTO align_exons(alignID,serial,gen_start,gen_stop,ref_start,ref_stop,ref_score) VALUES('$big_alignID','$exonserial','$genstart','$genstop','$refstart','$refstop','$refscore')";
        if($SQLverbose) { print $query,"\n"; }
        $dbh->do($query);
      }
      else { last; }

      # Fill align_introns table
      until (($line =~ m/<intron /) ||
             ($line =~ m/<\/exon-intron_info/)
            ) { $line = <>; }

      if ($line =~ m/<intron /) {
        $line =~ m/ i_serial=\"(.*?)\"/; $intronserial = $1;

        until ($line =~ m/<gDNA_intron_boundary /) { $line = <>; }
        $line =~ m/ i_start=\"(.*?)\" i_stop=\"(.*?)\"/;
        ($intstart,$intstop) = ($1,$2);

        until ($line =~ m/<donor /) { $line = <>; }
        if ($line =~ m/d_score/) {
          $line =~ m/ d_prob=\"(.*?)\" d_score=\"(.*?)\"/;
          ($donprob,$donscore) = ($1,$2);
        }
        else {
          $line =~ m/d_prob=\"(.*?)\"/;
          $donprob=$1;
          $donscore = $DNE_flag;
        } # protein case

        until ($line =~ m/<acceptor /) { $line = <>; }
        if ($line =~ m/a_score/) {
          $line =~ m/a_prob=\"(.*?)\" a_score=\"(.*?)\"/;
          ($accprob,$accscore) = ($1,$2);
        }
        else { 
          $line =~ m/a_prob=\"(.*?)\"/;
          $accprob=$1;
          $accscore = $DNE_flag;
        }

        $query = "INSERT INTO align_introns(alignID,serial,int_start,int_stop,don_prob,don_score,acc_prob,acc_score) VALUES('$big_alignID','$intronserial','$intstart','$intstop','$donprob','$donscore','$accprob','$accscore')";   
        if($SQLverbose) { print $query,"\n"; }
        $dbh->do($query);
      }
      else { last; }
    } # end while

    # Finish parsing information for the alignments table record
    until ($line =~ m/<total_alignment_score/) { $line = <>; }
    $line =~ m/<total_alignment_score>(.*?)<\/total_alignment_score>/;
    $tot_score = $1;

    until ($line =~ m/<cumulative_length_of_scored_exons/) { $line = <>; }
    $line =~ m/_scored_exons>(.*?)<\/cumulative_length/;
    $length_exons = $1;

    until ($line =~ m/<coverage /) { $line = <>; }
    $line =~ m/ percentage=\"(.*?)\" high_type=\"(.*?)\"/;
    ($coverage,$type_high) = ($1,$2);

    # If there is high quality PGS info, load that, too.
    $hqPGS="FALSE";
    until (($line =~ m/<hqPGS_line/)        ||
           ($line =~ m/<spliced_alignment/) ||
           ($line =~ m/<\/alignment_module/)) { $line = <>; }

    if ($line =~ m/<hqPGS_line/) {
      $hqPGS="TRUE";
      until ($line =~ m/<hqexon /) { $line = <>; }
      my $exonct=0;
      while ($line =~ m/<hqexon /) {
        $line =~ m/e_start=\"(.*?)\" e_stop=\"(.*?)\"/;
        ($estart,$estop) = ($1,$2);
        ++$exonct;
        $query = "INSERT INTO hiqual_align_exons(alignID,serial,gen_start,gen_stop) VALUES('$big_alignID','$exonct','$estart','$estop')";
        if($SQLverbose) { print $query,"\n"; }
        $dbh->do($query);
        $line=<>;
      }
    }

    $query = "UPDATE alignments SET total_score='$tot_score',length_exons='$length_exons',tot_exons='$tot_exons',coverage='$coverage',type_high_cov='$type_high',hiqual_bool='$hqPGS' WHERE alignID = '$big_alignID'";
    if($SQLverbose) { print $query,"\n"; }
    $dbh->do($query);
  } # end PARSE_ALIGNMENT_MODULE

  # Parse in the predgenloc information
  PARSE_PGL_MODULE: while(1) {
    until ($line =~ m/<PGL_line /) { $line = <>; }
    $line =~ m/PGL_start=\"(.*?)\" PGL_stop=\"(.*?)\"/;
    ($pgl_start,$pgl_stop) = ($1,$2);

    $query = "INSERT INTO predgenloc(run,pglID,pgl_start,pgl_stop) VALUES('$big_run',NULL,'$pgl_start','$pgl_stop')";
    if($SQLverbose) { print $query,"\n"; }
    $dbh->do($query);

    # Get $big_pglID, used to tie together the tables
    $sth = $dbh->prepare("SELECT LAST_INSERT_ID() FROM predgenloc");
    $sth->execute();
    @ary = $sth->fetchrow_array();
    $big_pglID = shift(@ary);
    $sth->finish();

    # Add what information we can into the altgenstruct table
    $query = "INSERT INTO altgenstruct(pglID,agsID) VALUES('$big_pglID',NULL)";
    if($SQLverbose) { print $query,"\n"; }
    $dbh->do($query);

    # Get $big_agsID, used to tie together the tables
    $sth = $dbh->prepare("SELECT LAST_INSERT_ID() FROM altgenstruct");
    $sth->execute();
    @ary = $sth->fetchrow_array();
    $big_agsID = shift(@ary);
    $sth->finish();

    # ags_exon code module
    until ($line =~ m/<exon-intron_info /) { $line = <>; }
    until ($line =~ m/<exon /) { $line = <>; }
    $line =~ m/e_serial=\"(.*?)\" e_score=\"(.*?)\"/;
    ($ags_exon_serial,$escore) = ($1,$2);

    until ($line =~ m/<gDNA_exon_boundary /) { $line = <>; }
    $line =~ m/e_start=\"(.*?)\" e_stop=\"(.*?)\"/;
    ($estart,$estop) = ($1,$2);
    $query = "INSERT INTO ags_exons(agsID,serial,e_start,e_stop,e_score) VALUES('$big_agsID','$ags_exon_serial','$estart','$estop','$escore')";
    if($SQLverbose) { print $query,"\n"; }
    $dbh->do($query);
 
    # Intron-exon code module
    while(1) {
      until (($line =~ m/<intron /) ||
             ($line =~ m/<\/exon-intron_info>/)
            ) { $line = <>; }

      if ($line =~ m/<\/exon-intron_info>/) { last; }
      else {
        # ags_intron code module
        if ($line =~ m/don_prob/) { # ESTcDNA case
          $line =~ m/i_serial=\"(.*?)\" don_prob=\"(.*?)\" acc_prob=\"(.*?)\"/;
          ($ags_intron_serial,$donor_prob,$acceptor_prob) = ($1,$2,$3); 
        }
        else { # Protein case
          $line =~ m/i_serial=\"(.*?)\"/;
          $ags_intron_serial=$1;
          ($donor_prob,$acceptor_prob) = $DNE_flag;
        }

        until ($line =~ m/<gDNA_intron_boundary /) { $line=<>; }
        $line =~ m/i_start=\"(.*?)\" i_stop=\"(.*?)\"/;
        ($istart,$istop) = ($1,$2);

        $query = "INSERT INTO ags_introns(agsID,serial,i_start,i_stop,don_prob,acc_prob) VALUES('$big_agsID','$ags_intron_serial','$istart','$istop','$donor_prob','$acceptor_prob')";
        if($SQLverbose) { print $query,"\n"; }
        $dbh->do($query);

        # ags_exon code module
        until ($line =~ m/<exon /) { $line = <>; }
        $line =~ m/e_serial=\"(.*?)\" e_score=\"(.*?)\"/;
        ($ags_exon_serial,$escore) = ($1,$2);

        until ($line =~ m/<gDNA_exon_boundary /) { $line = <>; }
        $line =~ m/e_start=\"(.*?)\" e_stop=\"(.*?)\"/;
        ($estart,$estop) = ($1,$2);
        $query = "INSERT INTO ags_exons(agsID,serial,e_start,e_stop,e_score) VALUES('$big_agsID','$ags_exon_serial','$estart','$estop','$escore')";
        if($SQLverbose) { print $query,"\n"; }
        $dbh->do($query);
      }
    } # end while

    # Find the supporting expression data
    $support = "";
    GET_SUPPORT: while(1) {
      until (($line =~ m/<PGS_line>/) ||
             ($line =~ m/<\/supporting_evidence>/)
            ) { $line = <>; }

      if ($line =~ m/<\/supporting_evidence>/) { last GET_SUPPORT; }
      else {
        until ($line =~ m/<reference/) { $line = <>; }
        if ($line =~ m/<referenceDNA/) {
          $line =~ m/ id=\"(.*?)\" strand=\"(.*?)\"/;
          $support .= "$1$2,";
        }
        else { # ($line =~ m/<referenceProtein/)
          $line =~ m/ id=\"(.*?)\"/;
          $support .= "$1,";
        }
      }
    } # end GET_SUPPORT 
    chop($support); # Remove trailing comma
   
    # Finish adding information to the altgenstruct record
    $query = "UPDATE altgenstruct SET support='$support',tot_exons='$ags_exon_serial' WHERE pglID = '$big_pglID'";
    if($SQLverbose) { print $query,"\n"; }
    $dbh->do($query);

    # Parse information for the prob_orfs table
    PARSE_THREE_PHASE:  while(1) { 
      until (($line =~ m/<three_phase_translation /) ||
             ($line =~ m/<\/AGS_information>/)
            ) { $line = <>; }

      if ($line =~ m/<\/AGS_information>/) { last PARSE_THREE_PHASE; }
      PARSE_PROB_ORF: while(1) {
        until (($line =~ m/<probable_ORFs /) ||
               ($line =~ m/<\/three_phase_translation>/)
              ) { $line = <>; } 

        if ($line =~ m/<\/three_phase_translation>/) { next PARSE_THREE_PHASE; }
        $line = <>;
        if ($line =~ m/<none/) { last PARSE_PROB_ORF; }
        else { # A probable_ORF(s) presents!
          PARSE_ORF_ENTRY: while(1) {
            until (($line =~ m/<orf_entry>/) ||
                   ($line =~ m/<\/probable_ORFs>/)
                  ) { $line = <>; }

            if ($line =~ m/<\/probable_ORFs>/) { last PARSE_ORF_ENTRY; }
            until ($line =~ m/<gDNA /) { $line = <>; }
            $line =~ m/id=\"(.*?)\" strand=\"(.*?)\"/;
            ($pps_temp_id,$pps_temp_strand) = ($1,$2);
            $e_coors = "";
            until (($line =~ m/<exon /) ||
                   ($line =~ m/<\/exon_boundaries>/)
                  ) { $line = <>; }

            while(1) {
              if ($line =~ m/<\/exon_boundaries>/) { last; }
              else {  
                $line =~ m/start=\"(.*?)\" stop=\"(.*?)\"/;
                $e_coors .= "($1 $2)";
                $line = <>;
              }
            }
            until ($line =~ m/<frame/) { $line = <>; }    
            $line =~ m/<frame>(.*?)<\/frame>/; $ppsframe = $1;
            until ($line =~ m/amino_acids/) { $line = <>; }
            $line =~ m/_amino_acids>(.*?)<\/number_/; $peplength = $1;

            until ($line =~ m/<predicted_protein_sequence>/) { $line = <>; }
            $line =~ m/<predicted_protein_sequence>(.*?)<\/predicted_protein_sequence>/; $ppseq = $1;

            $query = "INSERT INTO prob_orfs(agsID,templateID,template_strand,exon_coors,frame,pep_length,pep_seq) VALUES('$big_agsID','$pps_temp_id','$pps_temp_strand','$e_coors','$ppsframe','$peplength','$ppseq')";
            if($SQLverbose) { print $query,"\n"; }
            $dbh->do($query);
          } # end PARSE_ORF_ENTRY
        } # end else
      } # end PARSE_PROB_ORF
    } # end PARSE_THREE_PHASE

    until (($line =~ m/<predicted_gene_location>/) ||
           ($line =~ m/<\/PGL_module/)
          ) { $line = <>; }

    if ($line =~ m/<predicted_gene_location>/) { next PARSE_PGL_MODULE; }
    else { # ($line =~ m/<\/PGL_module/)
      next MAIN;
    }
  } # end PARSE_PGL_MODULE
} # end MAIN

# exit occurs in MAIN loop
