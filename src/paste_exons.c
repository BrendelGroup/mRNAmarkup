/* PASTE_EXONS.C;                              Last update: December 19, 2010 */
/*   - a subroutine to assemble exons into cDNA for translation.              */
/* Dependencies:   called by dnatopro.c, genestat.c.                          */
/* Bugs:                                                                      */


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



#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "def.h"
extern FILE *outfp, *tabfp;
extern int tabflag;
extern char sfname[81], outfname[60];
extern int dna[DNALGTH];
extern int protein[PROTLGTH];
extern int af[23], pf[9];
extern float aq[20][4], pq[9][4];
extern int Mflag, pstyle;
extern char NAUC[12], AAUC[25], CHPN[25];
extern int MINORFL;
extern int orfnum;
extern int codtoaa[4][4][4];
#define MAXORFLGTH 300000

paste_exons(numbp,cflag)   /* assembles exons into the mature message */
int numbp, cflag;
{
int i, id, ip;
int cod[MAXORFLGTH];
int iexon, exonst, exoned;
int gnum= 0;
int glgth, numaa;
int lmod;
int answer;
static char corfn[4]={'0','0','0','\0'};
char proname[60];
FILE *oldfp;

oldfp= outfp;

newsegm:

fprintf(stderr,"\nDo you wish to assemble an open reading frame ?");
fprintf(stderr,"\n ('q' to exit, 'return' to continue) \n");
answer= getchar();
if (answer=='q')   return;

glgth= 0; ip= 0; iexon= 0;
++gnum; ++orfnum;

if (cflag==1)
  {strcpy(proname,outfname);
   corfn[0]= (char)(orfnum/100 +48);
   corfn[1]= (char)((orfnum-(orfnum/100)*100)/10 +48);
   corfn[2]= (char)(orfnum%10 +48);
   strcat(proname,corfn);
   if ( (outfp= fopen( proname, "w" )) == NULL )
     {fprintf(stderr,"File %s cannot be opened.\n", proname);
      perror(proname); exit(-1);}
   fprintf(stderr,"File %s has been opened for writing.\n", proname);
  }
  
if (cflag==2)   fprintf(outfp,"\n\n________________________________________\n");
fprintf(outfp,"%s gene %3d:", sfname, gnum);

while (1)
  {fprintf(stderr,
	"\n Enter position of first base of exon ('0' to end) !\n");
  scanf("%d",&exonst);
  if (exonst==0)   {answer= getchar(); break;}
  fprintf(stderr,"\n Enter position of last base of exon !\n");
  scanf("%d",&exoned);
  answer= getchar();   /* ... this absorbs the return character ... */

  if (exonst<exoned)
    {fprintf(outfp,"\nExon %2d:  from %6d to %6d; length: %4d",
                  ++iexon, exonst, exoned, exoned-exonst+1);
     for( i=exonst; i<=exoned; ++i )   {cod[ip]= dna[i-1]; ip++;}
     glgth+= exoned-exonst+1;
    }
  else
    {fprintf(outfp,"\nExon %2d:  from %6d to %6d; length: %4d",
                  ++iexon, exonst, exoned, exonst-exoned+1);
     for( i=exonst; i>=exoned; --i )   {cod[ip]= (dna[i-1]+2)%4; ip++;}
     glgth+= exonst-exoned+1;
    }
 }

lmod= glgth%3;
if (lmod)
  {fprintf(stderr,"\n... exons do not add up to 0 modulo 3 !\n");
   fprintf(outfp,"\n... exons do not add up to 0 modulo 3 !\n");
   --gnum; --orfnum;
   goto newsegm;
  }

numaa= glgth/3;
if (numaa==0)   goto newsegm;

fprintf(outfp,"\nTotal length:  %5d bp, %4d residues\n", glgth, numaa);
if (cflag==1)   fprintf(outfp,"SQ");

id= 0;
for( ip=0;ip<numaa;++ip )
  {protein[ip]= codtoaa[cod[id]][cod[id+1]][cod[id+2]];
   id+=3;
  }
 
if (pstyle==4)
  {fprintf(outfp,"\n");
   for( i=0;i<glgth;++i )
     {if(i%10 == 0)  fprintf(outfp, " ");
      if(i%60 == 0)  fprintf(outfp, "\n%8d  ",i+1);
      fprintf(outfp, "%c", NAUC[cod[i]] );
     }
   fprintf(outfp,"\n");
  }

pr_pro(numaa,AAUC);

if (numaa>PROTLGTH)
  {fprintf(outfp,
   "\nORF exceeds PROTLGTH; redefine PROTLGTH in genestat.c.\n\n"); exit(0);}
if (cflag==2)
 {resuse(protein,numaa,af,pf,outfp,pstyle,1,aq,pq,0,tabfp,tabflag,AAUC,CHPN);
  coduse(outfp,cod,glgth,glgth,1);
  propuse(outfp,numaa,af,1);
 }
fprintf(outfp,"//\n");
if (cflag==1)   {fclose(outfp);   outfp= oldfp;}

goto newsegm;

} /* end paste_exons() */
