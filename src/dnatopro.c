/* DNATOPRO.C;                                Last update: January 1, 2011. */
/*   - a program to find and translate open reading frames in nucleic acid    */
/*   sequences.                                                               */
/* Dependencies:   getlibn.c getgbs.c getindn.c find_orfs.c paste_exons.c     */
/*                 pr_pro.c                                                   */
/* Bugs:                                                                      */


#define USAGE "Usage:\
 %s [-tv] [-a from] [-b to] [-rR] [-fF x] [-p x] [-M] [-o out_fname]\n\
    [-l libfname] [-s lstfname] [-g] seqfname(s)\n\n\
  -t: print nucleotide sequence with ORF coordinates indicated in the\n\
       header [default: print translated sequences]; in combination with\n\
       -F, the coding strand of the sequence is displayed, whereas in in\n\
       combination with -f the original sequence is displayed.\n\
  -v: in combination with -f, print coding sequence in addition to translated\n\
       sequence; in combination with -F, print coding sequence only\n\
  -a from:   from position\n\
  -b to:     to position\n\
  -r: analyze reverse strand only\n\
  -R: analyze both strands\n\
  -f x:   find and process open reading frames exceeding x codons in length\n\
          [by default open reading frames are interactively assembled from\n\
	   user specified exon coordinates]\n\
  -F x:   as -f, but output will be a single library file of input sequences\n\
	   annotated with the maximal qualifying open reading frame coordinates\n\
  -p x:   find open reading frames in sequence phase x only;\n\
	  [options: x= 1, 2, 0; default: all phases]\n\
  -M:  begin open reading frames with AUG codons only\n\
  -o out_fname:  print proteins into out_fname_xyz [default: PRO_xyz]\n\
  -l libfname:  read nucleic acid sequence data from library file libfname.\n\
  -s lstfname:  read nucleic acid sequence data from files specified in\n\
                the list file LST_lstfname.\n\
  -g         :  the following sequences are in GenBank format.\n\
  seqfname(s):  read nucleic acid sequence data from individual file(s).\n"


/*   Volker Brendel, Department of Genetics, Development and Cell Biology,    */
/*   Iowa State University,                                                   */
/*   Ames, IA 50010-3260; (515) 294-9884, vbrendel@iastate.edu                */

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

int pstyle= 2, tabflag= 0;
int dna[DNALGTH];
int protein[PROTLGTH];
FILE *infp, *outfp, *tabfp;
int af[23], pf[9];
float aq[20][4], pq[9][4];
char sfname[256], outfname[60]= "PRO";
int fcount= 0, orfnum= 0, orftlgth= 0;
int Mflag= 0, fdorfs= 0;
int MINORFL;

main(argc,argv)
int argc;
char *argv[];
{
int i, numbp;
int anum= 0;
int fposset= 0, FROMPOS, tposset= 0, TOPOS, rflag= 0, frflag= 9;
FILE *listfp;
static char lstname[60], libname[60];
static char inpline[LINELGTH];
int libop= 0, listop= 0;
static char gbftstrg[60]= "all";
int gbflag= 0, calln, nbaf, ftc[MAXNBAF][2];
time_t tlc;
outfp= stdout;

time(&tlc);

if (argc<2)   {fprintf(stderr,USAGE,argv[0]); exit(0);}
 
for( i=1; i<argc; ++i )
  {if ( *argv[i]=='-' )   switch ( *(argv[i]+1) )   {
        case 't':
                pstyle= 1; ++anum; break;
        case 'v':
                pstyle= 4; ++anum; break;
        case 'a':
                FROMPOS= atoi(argv[i+1])-1; fposset= 1; anum+= 2; break;
        case 'b':
                TOPOS= atoi(argv[i+1])-1; tposset= 1; anum+= 2; break;
        case 'r':
                rflag= 1; ++anum; break;
        case 'R':
                rflag= 2; ++anum; break;
        case 'f':
		fdorfs= 1; MINORFL= atoi(argv[i+1]); anum+= 2; break;
        case 'F':
		fdorfs= 2; MINORFL= atoi(argv[i+1]); anum+= 2; break;
        case 'p':
		frflag= atoi(argv[i+1]); anum+= 2; break;
        case 'M':
                Mflag= 1; ++anum; break;
	case 'o':
		strcpy(outfname,argv[i+1]); anum+= 2; break;
        case 'l':
                libop= 1;
                strcat(libname,argv[i+1]);
                if ( (infp= fopen( libname, "r" )) == NULL )
                   {fprintf(stderr,"File %s cannot be opened.\n", libname);
                    perror(libname); exit(-1);}
                fprintf(stderr,
		  "Library file %s has been opened for reading.\n\n", libname);
                anum+= 2; break;
        case 's':
                listop= 1;
                strcat(lstname,"LST_"); strcat(lstname,argv[i+1]);
                if ( (listfp = fopen( lstname, "r" )) == NULL )
                   {fprintf(stderr,"File %s cannot be opened.\n", lstname);
                    perror(lstname); exit(-1);}
                fprintf(stderr,
		  "File %s has been opened for reading.\n\n", lstname);
                fgets( inpline, LINELGTH, listfp );
                fgets( inpline, LINELGTH, listfp );
                anum+= 2; break;
        case 'g':
                gbflag= 1; ++anum; break;
        default:
                break;
        }
  }

fprintf(stderr,"\nNOW EXECUTING:  ");
for( i=0; i<argc; ++i )   fprintf(stderr," %s", argv[i]);
fprintf(stderr,"\n");

fprintf(outfp,"\nDNATOPRO.   Version of January 1, 2011.\n");
fprintf(outfp,"Date run: %s\n", ctime(&tlc) );

if (fdorfs==2)
 {if ( (outfp= fopen( outfname, "w" )) == NULL )
   {fprintf(stderr,"File %s cannot be opened.\n", outfname);
    perror(outfname); exit(-1);}
  fprintf(stderr,"File %s has been opened for writing.\n", outfname);
 }

if (libop)
 {if (fdorfs<2) fprintf(outfp, "Library file: %s \n\n", libname );
  while (numbp=getlibn(dna,1,infp))
   {if (++fcount%100 == 0)   fprintf(stderr,
	" ... now processing nucleotide sequence %4d (%s)\n", fcount, sfname);
    if (!fposset || FROMPOS < 1) FROMPOS= 1;
    if (!tposset || TOPOS > numbp) TOPOS= numbp;
    doit(numbp,FROMPOS,TOPOS,rflag,frflag);
   }
 }
 
while (1)
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
  if (gbflag)
   {calln= 0;
    while (numbp= getgbs(dna,0,&calln,infp,sfname,gbftstrg,&nbaf,ftc,0,0,0))
     {if (numbp<0) continue;
      if (!fposset || FROMPOS < 1) FROMPOS= 1;
      if (!tposset || TOPOS > numbp) TOPOS= numbp;
      doit(numbp,FROMPOS,TOPOS,rflag,frflag);
     }
   }
  else
   {numbp= getindn();
    if (!fposset || FROMPOS < 1) FROMPOS= 1;
    if (!tposset || TOPOS > numbp) TOPOS= numbp;
    doit(numbp,FROMPOS,TOPOS,rflag,frflag);
   }
  fclose(infp);
 }
 
if (fcount>1)   fprintf(stderr,"\n%5d files analyzed.\n", fcount);

fprintf(stderr,"\n");
exit(0);

} /* end main() */

 
 
doit(numbp,ia,ib,rflag,frflag)
int numbp, ia, ib, rflag, frflag;
{
int maxorflgth;
 
if (fdorfs) {
 if (fdorfs == 1) find_orfs(dna,numbp,ia,ib,MINORFL,rflag,frflag,1,pstyle);
 else {
  maxorflgth= find_orfs(dna,numbp,ia,ib,MINORFL,rflag,frflag,9,pstyle);
  if (maxorflgth > 0)
   find_orfs(dna,numbp,ia,ib,maxorflgth,rflag,frflag,0,pstyle);
 }
}
else   paste_exons(numbp,1);

} /* end doit() */
