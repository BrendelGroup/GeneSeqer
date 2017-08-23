__GeneSeqer/SplicePredictor.                       Version 5.1 (February 15, 2014)__
________________________________________________________________________________

__UTILITY:__ Gene prediction by spliced alignment

__LIMITATIONS:__
Currently implemented are splice site models for 'human',
'mouse', 'rat', 'chicken', 'Drosophila', 'nematode', 'yeast',
'Aspergillus', 'Arabidopsis', 'maize', 'rice', and "Medicago'.
For highly diverged  species, these models may not work equally
well, however we recommend to always use the most closely
related species (the "generic" species choice should work in
case of strong sequence similarity, but the performance gain
because of differential splice site scoring will be lost).

Potential gene structures are limited to an overall length of
15 kb in the genomic sequence [default setting; the limit may be
changed at compilation time by changing parameters in
_include/sahmt.h_, provided the CPU has sufficient RAM].  The
spliced alignment display of genes containing long introns may
be fragmented.

Default performance is optimized for typical EST lengths (less
than 2 kb). cDNAs up to 8 kp can be aligned [default limit; the
limit may be changed at compilation time by changing parameters
in include/sahmt.h], but the alignment algorithm will be slow.
Also, chance matches leading to insignificant alignments will be
more prevalent for long cDNAs.  For such applications, use of
the options "-x", "-y", and "-z" with higher than default values
is recommended.  

__DATABASE SUPPORT:__
Please see 0README in the DBsupport directory for
scripts and examples of database storage of GeneSeqer results.


Please direct all communications related to this software to

 Volker Brendel
 Indiana University
 Department of Biology
 212 South Hawthorne Drive
 Simon Hall 205C
 Bloomington, IN 47405, U.S.A.

 email:		vbrendel@indiana.edu
 Tel.:		(812) 855-7074


__INSTALLATION:__
```bash
cd src
make [system-type]	#(options: aix, hp_ux, linux, sgi, sun-os, win)
make install		#(default location for binaries: /usr/local/bin)
cd ../demo		#(see 0README_EXAMPLES)

#(Note: To install the binaries in a different directory, change the value
#       of the INSTALLDIR variable in Makefile. Alternatively, provide the
#       value on the commandline, e.g.
#
$make INSTALLDIR=~/bin install
$)
```

__BINARIES:__
The file _include/sahmt.h_ sets a number of parameters for compilation. The
provided makefiles invoke three different choices: _sahmt.hLRG_, _sahmt.hMED_
(= _sahmt.hDEF_), and _sahmt.SML_ and generate corresponding binaries
_GeneSeqerL_, _GeneSeqerM_ (= _GeneSeqer_), and _GeneSeqerS_ (or
_GeneSeqerMPI[l|m|s]_; see below).
Explore what works best for your system and application.
The difference is the size of spliced alignment attempted: \*L will allow a gene
size up to 100000 nucleotides and cDNA length up to 30000 nucleotides, whereas
\*S limits the sizes to 15000 and 8000, respectively.
_GeneSeqer_ will print out a WARNING message if it reaches the specified limits.
If you see WARNING messages in your output, increased sizes could be tried -
which will take more memory and time.


__OPTIONAL PARALLEL IMPLEMENTATION USING Open MPI:__
 The GeneSeqer MPI version should run on any MPI implementation.  Here we
 describe installation based on Open MPI (http://www.open-mpi.org/).

 On a Linux system, you may be able to install the required libraries as
 follows:

```bash
sudo dnf install openmpi openmpi-devel
```
 or
```bash1
sudo apt-get install openmpi openmpi-devel
```
 Then
```bash
locate mpicc
mpicc -showme:compile
mpicc -showme:link
```
should tell you where the installed binary is located and how to compile and
link _GeneSeqer_ on your system.  On our system the correct commands as
recorded in _src/makefile.lnxMPI_ are

```bash
MPICC    = /usr/lib64/openmpi/bin/mpicc -I/usr/include/openmpi-x86_64 -pthread -m64 -L/usr/lib64/openmpi/lib -lmpi -ldl
```

Thus, to compile the MPI version of _GeneSeqer_, do the following:

```bash
cd src
#(edit makefile.lnxMPI to set MPIHOME to the correct directory for your system
#and make other changes as appropriate for your installation)
make linux
make clean
make -f makefile.lnxMPI
cd ../mpidata		(see README_EXAMPLES)
```

To make things work, you need to have the PATH variables set right.  For
example, if you use the bash shell, add to _.bash_profile_:

```bash
PATH=$PATH:/usr/lib64/openmpi/bin
LD_LIBRARY_PATH=/usr/local/lib
LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/usr/lib64/openmpi/lib
```

Almost there!  Also execute (or add to your .bashrc file):

```bash
export LD_PRELOAD=/usr/lib64/openmpi/lib/libmpi.so
```

After making these changes active (source .bash_profile; source .bashrc), you
should be able to follow the example as per 0README_EXAMPLES instructions.
