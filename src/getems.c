/* GETEMS.C;                                  Last update: December 19, 2010. */
/*   - a subroutine to read nucleic acid sequence(s) in EMBL file format.     */
/* Dependencies:   called by intron.c                                         */
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
#include "def.h"
extern FILE *outfp;
extern char sfname[257];

getems(seq,fp)
int *seq;
FILE *fp;
{
int i, j, k, ch, sline= 0;
char buf[LINELGTH];
long offset;

for (i=0;i<257;++i)   sfname[i]= '\0';

if (fgets(buf,LINELGTH,fp)==NULL) return(0);
if (buf[0]=='I' && buf[1]=='D' && buf[2]==' ' && buf[3]==' ' && buf[4]==' ')
  {for (i=5; buf[i]!='\n' && buf[i]!=' ';++i)   sfname[i-5]= buf[i];}

while (sline==0)
  {if (fgets(buf,LINELGTH,fp)==NULL)
    {fprintf(stderr,"\nSequence file not in EMBL format!\n"); exit(-1);}
   if (buf[0]=='S' && buf[1]=='Q' && buf[2]==' ' && buf[3]==' ' && buf[4]==' ')
    sline= 1;
  }
offset= ftell(fp);
fseek(fp,offset,0);

i= 0;
while ( (ch=fgetc(fp)) != '/' )
  {if (ch=='T' || ch=='U' || ch=='t' || ch=='u' )  {seq[i]= 0; ++i; continue;}
   if (ch=='C' || ch=='c' )   {seq[i]= 1; ++i; continue;}
   if (ch=='A' || ch=='a' )   {seq[i]= 2; ++i; continue;}
   if (ch=='G' || ch=='g' )   {seq[i]= 3; ++i; continue;}
   if (ch=='Y' || ch=='y' )   {seq[i]=10; ++i; continue;}
   if (ch=='R' || ch=='r' )   {seq[i]=10; ++i; continue;}
   if (ch=='S' || ch=='s' )   {seq[i]=10; ++i; continue;}
   if (ch=='W' || ch=='w' )   {seq[i]=10; ++i; continue;}
   if (ch=='K' || ch=='k' )   {seq[i]=10; ++i; continue;}
   if (ch=='M' || ch=='m' )   {seq[i]=10; ++i; continue;}
   if (ch=='B' || ch=='b' )   {seq[i]=10; ++i; continue;}
   if (ch=='D' || ch=='d' )   {seq[i]=10; ++i; continue;}
   if (ch=='H' || ch=='h' )   {seq[i]=10; ++i; continue;}
   if (ch=='V' || ch=='v' )   {seq[i]=10; ++i; continue;}
   if (ch=='N' || ch=='n' )   {seq[i]=10; ++i; continue;}
  }
ch=fgetc(fp); ch=fgetc(fp);
 
return(i);

} /* end getems() */
