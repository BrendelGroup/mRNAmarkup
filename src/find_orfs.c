/* FIND_ORFS.C;                                Last update: January 7, 2011 */
/*   - a subroutine to find open reading frames in nucleic acid sequences.    */
/* Dependencies:   called by seqtopro.c, genestat.c; calls coduse.c           */
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
extern char sfname[256], outfname[60];
extern int protein[PROTLGTH];
extern int af[23], pf[9];
extern float aq[20][4], pq[9][4];
extern int Mflag;
extern char NAUC[12], AAUC[25], CHPN[25];
extern int orfnum, orftlgth;
extern int codtoaa[4][4][4];

#define MAXORFLGTH 300000
int cod[MAXORFLGTH];
int cflag;
char proname[60];
int gnum;
FILE *oldfp;


find_orfs(seq,numbp,ia,ib,minorflgth,rflag,frflag,flag,pstyle)
	/* finds the maximal open reading frames >= minorflgth in frame frflag */
	/* 	frflag==9:	all frames */
	/* 	flag==9:	determine length of maximal ORF only */
int seq[], numbp, ia, ib, minorflgth, rflag, frflag, flag, pstyle;
{
int i, j, aa;
int last_stp[3];
int frame, l_orf, maxl_orf= -1;

cflag= flag;
gnum= 0;
oldfp= outfp;

if (rflag%2==0)
 {for (i=0;i<=2;++i)  last_stp[i]= -1;
  for (i=ia-1;i<=ib-minorflgth*3;++i)
   {if (seq[i]>3 || seq[i+1]>3 || seq[i+2]>3) aa= 22;
    else aa= codtoaa[seq[i]][seq[i+1]][seq[i+2]];
    if (Mflag && aa!=17)   continue;
    frame= (i+1)%3;
    if ( i<=last_stp[frame] )   continue;
    j= i;
    while (aa!=23)
     {j+=3;
      if (j>ib-3)   break;   
      if (seq[j]>3 || seq[j+1]>3 || seq[j+2]>3) aa= 22;
      else aa= codtoaa[seq[j]][seq[j+1]][seq[j+2]];
     }
    last_stp[frame]= j;
    l_orf= (j-i)/3;
    if (l_orf > maxl_orf) maxl_orf= l_orf;
    if (flag < 9) {
     if ( l_orf >= minorflgth && (frflag==9 || frame==frflag))
      orf(seq,numbp,i,j,frame,pstyle); 
    }
   }
 }

if (rflag>0)
 {for (i=0;i<=2;++i)  last_stp[i]= numbp;
  for (i=ib-1;i>=ia-1+minorflgth*3;--i)
   {if (seq[i]>3 || seq[i+1]>3 || seq[i+2]>3) aa= 22;
    else aa= codtoaa[(seq[i]+2)%4][(seq[i-1]+2)%4][(seq[i-2]+2)%4];
    if (Mflag && aa!=17)   continue;
    frame= (numbp-i)%3;
    if ( i>= last_stp[frame] )   continue;
    j= i;
    while (aa!=23)
     {j-=3;
      if (j<ia+1)   break;
      if (seq[j]>3 || seq[j-1]>3 || seq[j-2]>3) aa= 22;
      else aa= codtoaa[(seq[j]+2)%4][(seq[j-1]+2)%4][(seq[j-2]+2)%4];
     }
    last_stp[frame]= j;
    l_orf= (i-j)/3;
    if (l_orf > maxl_orf) maxl_orf= l_orf;
    if (flag < 9) {
     if ( l_orf >= minorflgth  && (frflag==9 || -frame==frflag))
      orf(seq,numbp,i,j,frame,pstyle);
    }
   }
 }

if (maxl_orf < minorflgth) maxl_orf= -1;
return(maxl_orf);
} /* end find_orfs() */



orf(seq,numbp,firstbp,lastbp,frame,pstyle)
	/* ... processes open reading frames ... */
int seq[], numbp, firstbp, lastbp, frame, pstyle;
{
int i;
int glgth, numaa= 0;
char zeichen;
int prlgth;
static char corfn[7]={'0','0','0','0','0','0','\0'};

++gnum; ++orfnum;

if (firstbp<lastbp)
 {if (lastbp<numbp-2) {lastbp+= 3; numaa= -1;}/* include stop codon */
  glgth= lastbp-firstbp; zeichen= '+';
  numaa+= glgth/3;
 }
else
 {if (lastbp>1) {lastbp-= 3; numaa= -1;}/* include stop codon */
  glgth= firstbp-lastbp; lastbp+= 2; zeichen= '-';
  numaa+= glgth/3;
 }
orftlgth+= 3*numaa;

strcpy(proname,outfname);
corfn[0]= (char)(orfnum/100000 +48);
corfn[1]= (char)((orfnum-(orfnum/100000)*100000)/10000 +48);
corfn[2]= (char)((orfnum-(orfnum/10000)*10000)/1000 +48);
corfn[3]= (char)((orfnum-(orfnum/1000)*1000)/100 +48);
corfn[4]= (char)((orfnum-(orfnum/100)*100)/10 +48);
corfn[5]= (char)(orfnum%10 +48);
strcat(proname,corfn);

if (cflag==1)
 {if ( (outfp= fopen( proname, "w" )) == NULL )
   {fprintf(stderr,"File %s cannot be opened.\n", proname);
    perror(proname); exit(-1);}
  fprintf(stderr,"File %s has been opened for writing.\n", proname);
 }
     
if (cflag==2 && pstyle!=1)   fprintf(outfp,"________________________________________\n");


if (cflag!=2 && pstyle==1)
 {if (cflag == 1  ||  zeichen == '+')
   {fprintf(outfp,">%s   CDS: %d %d (frame '%c%1d'; %5d bp, %4d residues)\n",
	sfname, firstbp+1, lastbp, zeichen, frame, glgth, numaa );
    pr_seq(outfp,seq,numbp,NAUC,0,2);
   }
  else
   {fprintf(outfp,">%s   CDS: %d %d (frame '%c%1d'; %5d bp, %4d residues) Reverse complementary sequence displayed\n",
	sfname, numbp-firstbp, numbp-lastbp+1, zeichen, frame, glgth, numaa );
    pr_seq(outfp,seq,numbp,NAUC,1,2);
   }
  return;
 }

if (cflag==0 || (cflag==2 && pstyle!=1))
 {fprintf(outfp,">gnl|%s|%s\tORF #%3d:", proname, sfname, gnum );
  fprintf(outfp," from %6d to %6d (frame '%c%1d'; %5d bp, %4d residues)\n",
	firstbp+1, lastbp, zeichen, frame, glgth, numaa );
 }
else if (cflag==1)
 {fprintf(outfp,"ID   %s|%s\tORF #%3d:", proname, sfname, gnum );
  fprintf(outfp," from %6d to %6d (frame '%c%1d'; %5d bp, %4d residues)\n",
	firstbp+1, lastbp, zeichen, frame, glgth, numaa );
 }

prlgth= glgth;
if (firstbp<lastbp)
 {for( i=firstbp; i<lastbp; ++i )   cod[i-firstbp]= seq[i];
 }
else
 {lastbp-= 2;
  for( i=firstbp; i>lastbp; --i )   cod[firstbp-i]= (seq[i]+2)%4;
 }

translate(seq,firstbp,lastbp,numaa);

if (cflag==0)
 {if (pstyle==4) pr_seq(outfp,cod,prlgth,NAUC,0,2);
  else           pr_seq(outfp,protein,numaa,AAUC,0,2);
 }
if (cflag==1)
 {if (pstyle==4) {fprintf(outfp,"CQ"); pr_seq(outfp,cod,prlgth,NAUC,0,1);}
  fprintf(outfp,"SQ");
  pr_seq(outfp,protein,numaa,AAUC,0,1);
  fprintf(outfp,"//\n");
 }
else if (cflag==2 && pstyle!=1)
  pr_seq(outfp,protein,numaa,AAUC,0,2);

if (numaa>PROTLGTH)
 {fprintf(outfp,
  "\nORF exceeds PROTLGTH; redefine PROTLGTH in def.h.\n\n"); exit(0);}
if (cflag==2)
 {if (pstyle==2)
   resuse(protein,numaa,af,pf,outfp,1,1,aq,pq,0,tabfp,tabflag,AAUC,CHPN);
  else if (pstyle==4)
   resuse(protein,numaa,af,pf,outfp,4,1,aq,pq,0,tabfp,tabflag,AAUC,CHPN);
  coduse(outfp,cod,glgth,prlgth,1);
  propuse(outfp,numaa,af,1);
 }
if (cflag==1)   {fclose(outfp);   outfp= oldfp;}

} /* end orf() */



translate(seq,firstbp,lastbp,numaa)   /* ... creates protein[]= protein sequence
			 		 in numerical translation ... */
int seq[], firstbp, lastbp, numaa;
{
int id, ip;

id= firstbp;
if (firstbp <= lastbp)
  {for( ip=0; ip<numaa;++ip )
     {if (seq[id]>3 || seq[id+1]>3 || seq[id+2]>3) protein[ip]= 22;
      else protein[ip]= codtoaa[seq[id]][seq[id+1]][seq[id+2]];
      id+=3;
     }
  }
else
  {for( ip=0; ip<numaa;++ip )
     {if (seq[id]>3 || seq[id-1]>3 || seq[id-2]>3) protein[ip]= 22;
      else
	protein[ip]= codtoaa[(seq[id]+2)%4][(seq[id-1]+2)%4][(seq[id-2]+2)%4];
      id-=3;
     }
  }

} /* end translate() */
