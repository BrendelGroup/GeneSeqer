#Makefile for GeneSeqer/SplicePredictor
#Last update: February 11, 2019.

ifeq ($(OS),Windows_NT)
	BINDIR     = ..\bin
	INSTALLDIR = "%programfiles(x86)%\GeneSeqer\bin"
else
	BINDIR     = ../bin
	INSTALLDIR = /usr/local/bin
endif

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
# Check if file exists, create it if it doesn't
ifeq ("$(wildcard $(BINDIR))", "")
	mkdir $(BINDIR)
endif
	@make -f makefile.win

OPENMPI:
	mkdir -p $(BINDIR)
	@make -f makefile.lnx
	@make clean
	@make -f makefile.lnxMPI


install:
ifeq ($(OS),Windows_NT)
	xcopy /f $(BINDIR)\* $(INSTALLDIR)
else
	cp $(BINDIR)/* $(INSTALLDIR)/
endif

clean:
ifeq ($(OS),Windows_NT)
	del /f *.obj
else
	\rm *.o
endif

allclean:
ifeq ($(OS),Windows_NT)
	del /f *.o ..\bin\GeneSeqer* ..\bin\MakeArray ..\bin\SplicePredictor*
else
	\rm *.o ../bin/GeneSeqer* ../bin/MakeArray ../bin/SplicePredictor*
endif
