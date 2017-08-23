#Makefile for GeneSeqer/SplicePredictor
#Last update: February 15, 2014.

BINDIR     = ../bin
INSTALLDIR = /usr/local/bin

what:
	clear;
	@echo ""
	@echo "GeneSeqer/SplicePredictor make. Version 5.1 (February 15, 2014)."
	@echo ""
	@echo "Usage: make system-type"
	@echo "       make install"
	@echo "       make clean"
	@echo ""
	@echo "Known system types: aix, hp_ux, linux, sgi, sun-os, win"
	@echo ""
	@echo " Shortcuts: OPENMPI BIGRED2 MASON"
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

BIGRED2:
	mkdir -p $(BINDIR)
	@make -f makefile.lnxBIGRED2
	@make clean
	@make -f makefile.lnxBIGRED2MPI

MASON:
	mkdir -p $(BINDIR)
	@make -f makefile.lnx
	@make clean
	@make -f makefile.lnxMASONMPI


install:
	cp $(BINDIR)/* $(INSTALLDIR)/

clean:
	\rm *.o

allclean:
	\rm *.o ../bin/*
