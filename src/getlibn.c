/* GETLIBN.C;                                 Last update: December 19, 2010. */
/*   - a subroutine to read a nucleic acid sequence in library file format.   */
/* Dependencies:   called by dnacomp.c, dnaplot.c, dnatopro.c, dscore.c,      */
/*   genestat.c, getlibn.c                                                    */
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



#include <stdio.h>
#include <string.h>

#include "def.h"
extern FILE *outfp;
extern char sfname[256];

getlibn(seq,jbeg, fp)
int *seq, jbeg;
FILE *fp;
{
char buf[LINELGTH];
int n,j,l;

n=0;
for (j=0;j<257;++j)   sfname[j]= '\0';

fgets(buf,LINELGTH,fp);
if (buf[0]=='>')
  {for (j=jbeg;buf[j]!='\n' && buf[j]!=':';++j)   sfname[j-jbeg]= buf[j];}

while (fgets(buf,LINELGTH,fp) && buf[0]!='>')
  {l= strlen(buf);
   for (j=0;j<l;j++)
     {if ( buf[j]=='T' || buf[j]=='U' || buf[j]=='t' || buf[j]=='u' )
	{ seq[n++] =  0; continue; }
      if ( buf[j]=='C' || buf[j]=='c' )
	{ seq[n++] =  1; continue; }
      if ( buf[j]=='A' || buf[j]=='a' )
	{ seq[n++] =  2; continue; }
      if ( buf[j]=='G' || buf[j]=='g' )
	{ seq[n++] =  3; continue; }
      if ( buf[j]=='N' || buf[j]=='n' )
	{ seq[n++] = 10; continue; }
     }
  }

if (buf[0]=='>') fseek(fp,-strlen(buf),1);

return(n);

} /* end getlibn() */
