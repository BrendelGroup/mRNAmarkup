/* GENESTAT.C;                                Last update: December 19, 2010. */
/*   - a program to analyze codon usage in genes and open reading frames.     */
/* Dependencies:   getlibn.c getindn.c find_orfs.c paste_exons.c              */
/* Bugs:                                                                      */


#define USAGE "Usage:\
 %s [-tv] [-f x] [-M] [-b libfname] [-e] [-l lstfname] seqfname(s)\n\n\
  -t: terse output;   -v: verbose output (print gene sequence)\n\
  -f x:   find and process phase-1 open reading frames exceeding x codons in\n\
          length [by default open reading frames are interactively assembled\n\
	  from user specified exon coordinates]\n\
  -M:  begin open reading frames with AUG codons only\n\
  -b libfname:  read nucleic acid sequence data from library file libfname.\n\
  -e: subsequent files specify EMBL formatted data files.\n\
  -l lstfname:  read nucleic acid sequence data from files specified in\n\
                the list file LST_lstfname.\n\
  seqfname(s):  read nucleic acid sequence data from individual file(s).\n"
 

/*   Volker Brendel, Department of Genetics, Development and Cell Biology     */
/*   Iowa State University, Ames, IA 50010-3260;                              */
/*   (515) 294-9884, vbrendel@iastate.edu                                     */

/*******************************************************************************
Copyright (c) 2000 Volker Brendel
All Rights Reserved. E-mail: vbrendel@iastate.edu

Permission to use, copy, modify, and distribute this software and its
documentation for educational, research and non-profit purposes, without fee,
and without a written agreement is hereby granted, provided that the above
copyright notice, this paragraph and the following three paragraphs appear in
all copies. If you modify this file or included files you must cause the
modified files to carry prominent notices stating that you changed the files.

Inqueries for permission to incorporate this software into commercial products
should be directed to the Office of Intellectual Property and Technology
Transfer, 310 Lab of Mechanics, Iowa State University, Ames, IA 50011, phone:
(515) 294-4740, E-mail: Licensing@iastate.edu.

IN NO EVENT SHALL THE AUTHOR OR IOWA STATE UNIVERSITY BE LIABLE TO ANY PARTY FOR
DIRECT, INDIRECT, SPECIAL, INCIDENTAL, OR CONSEQUENTIAL DAMAGES, INCLUDING LOST
PROFITS, ARISING OUT OF THE USE OF THIS SOFTWARE AND ITS DOCUMENTATION, EVEN IF
IOWA STATE UNIVERSITY HAS BEEN ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

IOWA STATE UNIVERSITY SPECIFICALLY DISCLAIMS ANY WARRANTIES, INCLUDING, BUT NOT
LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
PARTICULAR PURPOSE.  THE SOFTWARE PROVIDED HEREUNDER IS ON AN "AS IS" BASIS,
AND IOWA STATE UNIVERSITY HAS NO OBLIGATIONS TO PROVIDE MAINTENANCE, SUPPORT,
UPDATES, ENHANCEMENTS, OR MODIFICATIONS.
*******************************************************************************/



#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

#include "abc.h"
#include "def.h"

#define MAXORFLGTH 300000

int pstyle= 2, tabflag= 0;
int dna[DNALGTH];
int protein[PROTLGTH];
FILE *infp, *outfp, *tabfp;
int af[23], pf[9];
float aq[20][4], pq[9][4];
char sfname[256], outfname[60];
int fcount= 0, gfcount= 0, emfcount= 0;
int orfnum= 0, orftlgth= 0;
int Mflag= 0, fdorfs= 0, emflag= 0;
int MINORFL;

main(argc,argv)
int argc;
char *argv[];
{
int i, numbp;
int anum= 0;
FILE *listfp;
static char lstname[60], libname[60];
static char emfname[60], inpline[LINELGTH];
int libop= 0, listop= 0;
time_t tlc;
outfp= stdout;

time(&tlc);

if (argc<2)   {fprintf(stderr,USAGE,argv[0]); exit(0);}
 
for( i=1; i<argc; ++i )
  {if ( *argv[i]=='-' )   switch ( *(argv[i]+1) )   {
        case 't':
                pstyle= 1;
                ++anum;
                break;
        case 'v':
                pstyle= 4;
                ++anum;
                break;
        case 'f':
		fdorfs= 1;
		MINORFL= atoi(argv[i+1]);
                anum+= 2;
                break;
        case 'M':
                Mflag= 1;
                ++anum;
                break;
        case 'b':
                libop= 1;
                strcat(libname,argv[i+1]);
                if ( (infp= fopen( libname, "r" )) == NULL )
                   {fprintf(stderr,"File %s cannot be opened.\n", libname);
                    perror(libname); exit(-1);}
                fprintf(stderr,
		"\nLibrary file %s has been opened for reading.\n\n", libname);
                anum+= 2;
                break;
        case 'e':
                emflag= 1;
                ++anum;
                break;
        case 'l':
                listop= 1;
                strcat(lstname,"LST_");
                strcat(lstname,argv[i+1]);
                if ( (listfp = fopen( lstname, "r" )) == NULL )
                   {fprintf(stderr,"File %s cannot be opened.\n", lstname);
                    perror(lstname); exit(-1);}
                fprintf(stderr,
		  "\nList %s has been opened for reading.\n\n", lstname);
                fgets( inpline, LINELGTH, listfp );
                fgets( inpline, LINELGTH, listfp );
                anum+= 2;
                break;
        default:
                break;
        }
  }


fprintf(stderr,"\nNOW EXECUTING:  ");
for( i=0; i<argc; ++i )   fprintf(stderr," %s", argv[i]);
fprintf(stderr,"\n");

fprintf(outfp,"GENESTAT.   Version of December 19, 2010.\n");
fflush(outfp);
fprintf(outfp,"Date run: %s\n", ctime(&tlc) );

if (fdorfs) fprintf(outfp,"Processing ORFs in frame 1 of at least %d nucleotides",MINORFL);
if (Mflag)  fprintf(outfp," and starting with ATG");
            fprintf(outfp,".\n\n");
fflush(outfp);

if (libop)
  {fprintf(outfp, "Library file: %s \n\n", libname );
   while (numbp=getlibn(dna,1,infp))
     {if (++fcount%100 == 0)   fprintf(stderr,
	" ... now processing nucleotide sequence %4d (%s)\n", fcount, sfname);
      doit(numbp);
     }
  }
 
if (emflag) while (1)
  {if (anum+1==argc)   break;
   ++anum;
   if ( (infp = fopen( argv[anum], "r" )) == NULL )
    {fprintf(stderr,"EMBL file %s cannot be opened.\n", argv[anum]);
     perror(argv[anum]); continue;}
   fprintf(stderr,"EMBL file %s has been opened for reading.\n", argv[anum]);
   ++emfcount;
   strcpy(emfname,argv[anum]);
   if (anum+1==argc && emfcount==1)
    fprintf(outfp,"EMBL file: %s\n", emfname);
   else
    {fprintf(outfp,"\n\n");
     for( i=0;i<80;++i )   fprintf(outfp,"*");
     fprintf(outfp,"\nEMBL file: %s (# %4d)\n", emfname, emfcount);
    }
   while (numbp=getems(dna,infp))
    {if (++fcount%100 == 0)   fprintf(stderr,
        " ... now processing nucleotide sequence %4d (%s)\n", fcount, sfname);
     doit(numbp);
    }
   fflush(outfp); fclose(infp);
  }
else while (1)
  {if (listop)
     {if ( fscanf(listfp, "%s", sfname) == EOF )
        {fclose(listfp); listop= 0; continue;}
      if ( (infp = fopen( sfname, "r" )) == NULL )
        {fprintf(stderr,"File %s cannot be opened.\n", sfname);
	 perror(sfname); continue;}
      if (++fcount%25 == 0)
	fprintf(stderr,"File %4d (%s) has been opened for reading.\n",
		fcount, sfname);
     }
   else
     {if (anum+1==argc)   break;
      ++anum;
      if ( (infp = fopen( argv[anum], "r" )) == NULL )
        {fprintf(stderr,"File %s cannot be opened.\n", argv[anum]);
	 perror(argv[anum]); continue;}
      fprintf(stderr,"File %s has been opened for reading.\n", argv[anum]);
      ++fcount;
      strcpy(sfname,argv[anum]);
     }
   numbp= getindn();
   doit(numbp);
   fclose(infp);
  }
 
if (fcount>1)   fprintf(stderr,"\n%5d files analyzed.\n", fcount);

fprintf(outfp,
  "\n________________________________________________________________________________\n");
fprintf(outfp,"Number of ORFs: %d.  Total length: %d nt.  Number of codons: %d.\n",
	orfnum,orftlgth,orftlgth/3);
coduse(outfp,dna,0,0,2);

fprintf(outfp,"\n");
exit(0);
}

 
 
doit(numbp)
int numbp;
{
 
if (fdorfs)   find_orfs(dna,numbp,1,numbp,MINORFL,0,1,2,pstyle);
else   paste_exons(numbp,2);

}
