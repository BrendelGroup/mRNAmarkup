/* GETGBS.C;                                  Last update: December 19, 2010. */
/*   - a subroutine to read nucleic acid sequence(s) in GenBank file format.  */
/* Dependencies:   called by exdomino.c, gbfeata.c, sgap.c                    */
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
extern FILE *outfp;
#define SHRTNFLAG 1

getgbs(seq,rflag,calln,fp,sfname,fstr,nbaf,ftc,lflag,deltal,deltar)
int *seq, rflag, *calln, *nbaf, ftc[MAXNBAF][2], lflag, deltal, deltar;
FILE *fp;
char sfname[], *fstr;
{
int i, j, k, ch, sline= 0, fdflag= 0, tmpn;
char buf[LINELGTH];
long offset;
static long c0offset;
int inb;
static char corfn[5]= {'0','0','0','0','\0'};

if (*calln == 0)   c0offset= ftell(fp);
else   fseek(fp,c0offset,0);
	/* IF *calln IS GREATER THAN 0, READ THE FILE AGAIN (FROM THE
	   LAST "LOCUS" LINE */

for (i=0;i<257;++i)   sfname[i]= '\0';
*nbaf= 0;

if (fgets(buf,LINELGTH,fp)==NULL) return(0);
while (buf[0]=='\n')
  if (fgets(buf,LINELGTH,fp)==NULL) return(0);

if (buf[0]=='L' && buf[1]=='O' && buf[2]=='C' && buf[3]=='U' && buf[4]=='S')
  {for (i=12;buf[i]!='\n' && buf[i]!=' ';++i)   sfname[i-12]= buf[i];}

if (*calln)
 {inb= *calln +1;
  corfn[0]= (char)(inb/1000 +48);
  corfn[1]= (char)((inb-(inb/1000)*1000)/100 +48);
  corfn[2]= (char)((inb-(inb/100)*100)/10 +48);
  corfn[3]= (char)(inb%10 +48);
  strcat(sfname,"_");
  strcat(sfname,corfn);
 }	/* APPEND THE *calln INTEGER TO THE FILE NAME */

while (sline==0)
 {if (fgets(buf,LINELGTH,fp)==NULL)
   {fprintf(stderr,"\nSequence file not in GenBank format!\n"); exit(-1);}
  if (buf[0]=='F' && buf[1]=='E' && buf[2]=='A' && buf[3]=='T' && buf[4]=='U')
   {while (1)
     {fscanf(fp,"%s",buf);
      if (strcmp(buf,"FEATDONE")==0) {++fdflag; continue;}
	/* "FEATDONE" IS THE SEPARATER BETWEEN INDEPENDENT DETERMINATIONS
	   OF FEATURES; E.G. BETWEEN SEPARATE CDSs UNDER ONE LOCUS NAME */
      if (strcmp(buf,fstr)==0)
       {if (*calln != fdflag) continue;
		/* ONLY RECORD THE FEATURE INFORMATION PERTINENT TO THE
		   CURRENT *calln READING OF THE SEQUENCE */
        if (*nbaf == MAXNBAF)
	 {fprintf(stderr,"\n\nMAXNBAF (%d) exceeded; recompile!\n", MAXNBAF);
	  exit(-1);
	 }
        fscanf(fp,"%d %d",&ftc[*nbaf][0],&ftc[*nbaf][1]);
	++(*nbaf);
       }
      if (strcmp(buf,"ORIGIN")==0) {fgets(buf,LINELGTH,fp); sline= 1; break;}
     }
    if (sline==1) break;
   }
  if (buf[0]=='O' && buf[1]=='R' && buf[2]=='I' && buf[3]=='G' && buf[4]=='I')
    sline= 1;
 }
if (fdflag)
 {if (*calln == 0 && strcmp(fstr,"all")!=0) strcat(sfname,"_0001");
  ++(*calln);
 }
offset= ftell(fp);
fseek(fp,offset,0);

i= 0;
while ( (ch=fgetc(fp)) != '/' )
 {if (i==DNALGTH)
   {fprintf(stderr,"\n\nDNALGTH (%d) exceeded in file %s. Skipped.\n",
			DNALGTH,sfname);
    fprintf(outfp,"\n\nDNALGTH (%d) exceeded in file %s. Skipped.\n",
			DNALGTH,sfname);
    *calln= 0;   return(-1);
   }
  if (ch=='T' || ch=='U' || ch=='t' || ch=='u' )  {seq[i]= 0; ++i; continue;}
  if (ch=='C' || ch=='c' )   {seq[i]= 1; ++i; continue;}
  if (ch=='A' || ch=='a' )   {seq[i]= 2; ++i; continue;}
  if (ch=='G' || ch=='g' )   {seq[i]= 3; ++i; continue;}
  if (ch=='Y' || ch=='y' )   {seq[i]=10; ++i; continue;}
  if (ch=='R' || ch=='r' )   {seq[i]=10; ++i; continue;}
  if (ch=='S' || ch=='s' )   {seq[i]=10; ++i; continue;}
  if (ch=='W' || ch=='w' )   {seq[i]=10; ++i; continue;}
  if (ch=='K' || ch=='k' )   {seq[i]=10; ++i; continue;}
  if (ch=='M' || ch=='m' )   {seq[i]=10; ++i; continue;}
  if (ch=='N' || ch=='n' )   {seq[i]=10; ++i; continue;}
  if (ch=='V' || ch=='v' )   {seq[i]=10; ++i; continue;}
 }
ch=fgetc(fp); ch=fgetc(fp);
if (rflag)
 {for (j=i-1;j>=(i+1)/2;--j)
   {if (seq[j]<=3) tmpn= (seq[j]+2)%4;   else tmpn= 10;
    if (seq[i-1-j]<=3) seq[j]= (seq[i-1-j]+2)%4;   else seq[j]= 10;
    seq[i-1-j]= tmpn;
   }
  if (i%2==1) seq[i/2]= (seq[i/2]+2)%4;
  if (rflag==1)
   {for (j=0;j<*nbaf;++j)
     {tmpn= ftc[j][0]; ftc[j][0]= i-ftc[j][1]+1; ftc[j][1]= i-tmpn+1;
     }
    for (j=*nbaf-1;j>=(*nbaf+1)/2;--j)
     {tmpn=  ftc[j][0]; ftc[j][0]= ftc[*nbaf-1-j][0]; ftc[*nbaf-1-j][0]= tmpn;
     }
    for (j=*nbaf-1;j>=(*nbaf+1)/2;--j)
     {tmpn=  ftc[j][1]; ftc[j][1]= ftc[*nbaf-1-j][1]; ftc[*nbaf-1-j][1]= tmpn;
     }
   }
 }

if (strcmp(fstr,"all")==0)
 {*nbaf= 1; ftc[0][0]= 1; ftc[0][1]= i; *calln= 0; return(i);}

if (*calln > fdflag+1) {*calln= 0; return(-1);}
	/* THIS CASE SIGNIFIES THAT THE LAST FEATURE BLOCK HAS ALREADY BEEN
	   READ; THUS -1 IS RETURNED TO THE CALLING PROGRAM, WHICH
	   WILL SKIP ON TO THE NEXT "LOCUS" BLOCK WITH RESET *calln */

for (j=k=0;j<*nbaf;++j,++k)
 {if (lflag<0)
   {if (ftc[j][0]-deltal < 1 || ftc[j][0]+deltar > i) {--k; continue;}
    if (ftc[j][0]<ftc[j][1])
     {ftc[k][1]= ftc[j][0]+deltar;
      ftc[k][0]= ftc[j][0]-deltal;
     }
    else
     {ftc[k][1]= ftc[j][0]-deltal;
      ftc[k][0]= ftc[j][0]+deltar;
     }
   }
  if (lflag==0)
   {if (ftc[j][0]-deltal < 1 || ftc[j][1]+deltar > i)
     {if (SHRTNFLAG)
       {if (ftc[j][0]-deltal < 1) ftc[k][0]= 1;
	else ftc[k][0]= ftc[j][0]-deltal;
        if (ftc[j][1]+deltar > i) ftc[k][1]= i;
	else ftc[k][1]= ftc[j][1]+deltar;
       }
      else  --k;
      continue;
     }
    ftc[k][0]= ftc[j][0]-deltal;
    ftc[k][1]= ftc[j][1]+deltar;
   }
  if (lflag>0)
   {if (ftc[j][1]-deltal < 1 || ftc[j][1]+deltar > i) {--k; continue;}
    if (ftc[j][0]<ftc[j][1])
     {ftc[k][0]= ftc[j][1]-deltal;
      ftc[k][1]= ftc[j][1]+deltar;
     }
    else
     {ftc[k][0]= ftc[j][1]+deltar;
      ftc[k][1]= ftc[j][1]-deltal;
     }
   }
 }
*nbaf= k;

return(i);

} /* end getgbs() */
