#Makefile VBtools (subset used in mRNAmarkup)
#Last update: April 8, 2015.

BINDIR     = ../bin
SRCDIR     = $(shell pwd)
INCLDIR    = $(SRCDIR:src=include)
INSTALLDIR = /usr/local/bin
# or the following, to install in /bin under your home directory:
#INSTALLDIR = $(HOME)/bin

CC     = gcc
CFLAGS = -DINCLDIR=\"$(INCLDIR)/\"
IFLAGS = -I$(INCLDIR)
LFLAGS = -lm

.KEEP_STATE:

PROGRAMS=	dnatopro \
		genestat

.c.o:
	$(CC) $(CFLAGS) $(IFLAGS) -c $<

all: $(PROGRAMS)


dnatopro:	dnatopro.o getlibn.o getgbs.o getindn.o find_orfs.o \
	paste_exons.o pr_seq.o pr_pro.o resuse.o coduse.o propuse.o
	$(CC) $(CFLAGS) $(IFLAGS) dnatopro.o getlibn.o getgbs.o getindn.o \
	find_orfs.o paste_exons.o pr_seq.o pr_pro.o resuse.o coduse.o propuse.o \
	$(LFLAGS) -o $(BINDIR)/dnatopro

genestat:	genestat.o getlibn.o getindn.o getems.o find_orfs.o \
	paste_exons.o pr_seq.o pr_pro.o resuse.o coduse.o propuse.o
	$(CC) $(CFLAGS) $(IFLAGS) genestat.o getlibn.o getindn.o getems.o \
	find_orfs.o paste_exons.o pr_seq.o pr_pro.o resuse.o coduse.o propuse.o \
	$(LFLAGS) -o $(BINDIR)/genestat


install:
	cp $(BINDIR)/dnatopro $(BINDIR)/genestat $(BINDIR)/*.pl $(BINDIR)/ESTScan.conf  $(INSTALLDIR)/

clean:
	\rm *.o

allclean:
	\rm *.o ../bin/dnatopro ../bin/genestat
