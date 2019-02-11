#! /usr/bin/perl -w
#
# fastasplit.pl
# Version of January 8, 2014.
#
# Please direct all communications related to this software to
#
#  Volker Brendel
#  Indiana University
#  Department of Biology
#  212 South Hawthorne Drive
#  Simon Hall 205C
#  Bloomington, IN 47405, U.S.A.
#
#  email:         vbrendel@indiana.edu
#  Tel.:          (812) 855-7074


use strict;
use Getopt::Std;


#------------------------------------------------------------------------------
my $USAGE="\nUsage: $0 -i fastafile -s psize [-n npieces] [-o outputdir]\n

** This script splits the input FASTA file fastafile into multiple pieces
   fastafile1, fastafile2, ... .  How many piecess are generated depends on the
   arguments of either the -s option (each piece is at most psize MB) or the
   -n option (npieces files will be generated).

   The original fastafile remains untouched.

   If an output directory is specified with the -o option, then the pieces
   will be written to the specified directory.
   \n";


my %args;
getopts('i:s:n:o:', \%args);

if (!defined($args{i})) {
  print "\n!! Please specify an input file.\n\n";
  die $USAGE;
}
my $fastafile = $args{i};
if (! -e $fastafile) {
  print "\n!! Input file $fastafile does not exist.\n\n";
  die $USAGE;
}


# Initiate variables:
# 
my $outputdir = ".";
my $seqHeader = "";		#The current sequence header
my $sequence = "";		#The current sequence
my $seqLength = 0.0;		#The current sequence length in Mb
my $curFileLength = 0.0;	#The current file length in Mb
my $fileCount = 1;		#The count of output files made
my $maxLength;			#The specified maximal file length in Mb
my $wasSeq = 0;			#Flag for sequence construction
my $npieces = 0;		#Number of files to be generated

if (defined($args{s})) {
  if ($args{s} !~ /\d+/) {
    print "\nArgument to -s should be integer (piece size in MB).\n\n";
    die $USAGE;
  }
  $maxLength = $args{s};
}
else {
  if (defined($args{n})) {
    if ($args{n} !~ /\d+/) {
      print "\nArgument to -n should be integer (number of pieces).\n\n";
      die $USAGE;
    }
    $npieces = $args{n};
  }
  else {
    print "\nSpecify either -s or -n option argument.\n\n";
    die $USAGE;
  }
}
if (defined($args{o})) {
  $outputdir = $args{o};
  if (! -e $outputdir) {
    print "\n!! Output directory $outputdir does not exist.\n\n";
    die $USAGE;
  }
}


# Ok, now either the -s or the -n argument is well-formed.  In the latter case,
# set the piece size to 1/npieces of the input file size:
#
if (!defined($args{s})) {
	open (IN, "<$fastafile");
	while (<IN>) {
		if ($_ =~ /^[>;\s]/) {
			next;
		}
		elsif ($_ =~ /(\w+)/) {
			$sequence .=  $&;
			next;
		}
		else {
			print "I don't know how the code got here.  Please check the input file for FASTA-format compatibility.\n";
			exit;
		}
	}
	close (IN);
	$maxLength =  1.005 * length($sequence) / ($args{n} * 1000000);
	#Note: We make $maxLength a bit longer, as otherwise we might
	#      inadvertently create an extra tiny piece due to rounding effects.
}
printf "\n\n$fastafile pieces of approximate cumulative length %6.2f Mb will be deposited in directory $outputdir/ ... \n\n", $maxLength;


#Open input file for real:
#
open (IN, "<$fastafile");
$sequence = "";

#Open first output file:
#
open (OUT, ">$outputdir/$fastafile$fileCount");
	
while (<IN>) {
	if ($_ =~ /^[>;\s]/) {		#Gets header and comments -> $seqHeader
		if ($wasSeq) {
			$wasSeq = 0;
			printToFile();
		}
		$seqHeader .= $_;	#Captures entire line, including white space
		next;
	}
	elsif ($_ =~ /(\w+)/) {		#Gets just the sequence, no whitespace
		$sequence .=  $&;
		$wasSeq = 1;
		next;
	}
}

#Print final sequence and close final file:
#
printToFile();
close OUT;
printf "File: $fastafile$fileCount"." created. Approximate cumulative length: %6.2f Mb\n", $curFileLength;



#SUBROUTINES:
#
sub printToFile {
	$seqLength = length($sequence)/1000000.0;
	if ( $curFileLength > 0.0  &&
	     ($seqLength + $curFileLength) > $maxLength  &&	#If maxLength will be exceeded, close current outFile and open a new one,
	     !($npieces > 0  &&  $fileCount == $npieces)   ) {	# unless the specified number of output files has already been reached.
		close OUT;
		printf "File: $fastafile$fileCount"." created. Approximate cumulative length: %6.2f Mb\n", $curFileLength;
		$curFileLength = 0.0;
		$fileCount++;
		open (OUT, ">$outputdir/$fastafile$fileCount");
	}
	
	#Write sequence to file and update $curFileLength:
	#
	print OUT "$seqHeader$sequence\n";
	$curFileLength += $seqLength;	
	$sequence = "";
	$seqHeader = "";
}
