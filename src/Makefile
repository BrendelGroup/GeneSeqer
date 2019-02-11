#Makefile for GeneSeqer/SplicePredictor
#Last update: February 11, 2019.

BINDIR     = ../bin
INSTALLDIR = /usr/local/bin

what:
	clear;
	@echo ""
	@echo "GeneSeqer/SplicePredictor make. Version 5.2 (February 11, 2019)."
	@echo ""
	@echo "Usage: make system-type"
	@echo "       make install"
	@echo "       make clean"
	@echo ""
	@echo "Known system types: aix, hp_ux, linux, sgi, sun-os, win"
	@echo ""
	@echo " Shortcuts: OPENMPI"
	@echo ""
	@echo "Default install directory:"
	@echo "  $(INSTALLDIR)"
	@echo ""
	@exit 1;

aix:
	mkdir -p $(BINDIR)
	@make -f makefile.aix

hp_ux:
	mkdir -p $(BINDIR)
	@make -f makefile.hpx

linux:
	mkdir -p $(BINDIR)
	@make -f makefile.lnx

sgi:
	mkdir -p $(BINDIR)
	@make -f makefile.sgi

sun-os:
	mkdir -p $(BINDIR)
	@make -f makefile.sun

win:
	mkdir -p $(BINDIR)
	@make -f makefile.win

OPENMPI:
	mkdir -p $(BINDIR)
	@make -f makefile.lnx
	@make clean
	@make -f makefile.lnxMPI


install:
	cp $(BINDIR)/* $(INSTALLDIR)/

clean:
	\rm *.o

allclean:
	\rm *.o ../bin/GeneSeqer* ../bin/MakeArray ../bin/SplicePredictor*
