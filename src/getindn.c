/* GETINDN.C;                                 Last update: December 19, 2010. */
/*   - a subroutine to read a nucleic acid sequence in individual file format.*/
/* Dependencies:   called by aaruns.c, mproc.c, pgrep.c, pmagic.c, pmotif.c,  */
/*   propeat.c, prorep.c, prostat.c, pscore.c, psearch.c, pword.c, saps.c,    */
/*   zipscan.c                                                                */
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
extern FILE *infp, *outfp;
extern int pstyle;
extern int dna[DNALGTH];

getindn()
{
int l, ch, numbp;
char inpline[LINELGTH];
int sline= 0;
long offset;
int i= 0;
for( l=0; l<LINELGTH; ++l )   inpline[l]= 'x';

/*FORMAT:
  The first line of the file is a header to be printed, the first line
  of the sequence is preceded by a line beginning with 'SQ', and subsequent
  symbols are A-U (one-letter-symbols) as part of the sequence or irrelevant
  characters (like numbers and blanks);  non-standard symbols for ambiguous
  or missing monomers are ignored.
*/

fgets(inpline,LINELGTH,infp);	
if (pstyle%2==0)   fputs(inpline,outfp);
while (sline==0)		
  {if (fgets(inpline,LINELGTH,infp)==NULL)
    {fprintf(stderr,"\nSequence file not in minimal EMBL format!\n"); exit(-1);}
   if (inpline[0]=='S' && inpline[1]=='Q')   sline= 1;
  }						
offset= ftell(infp);
fseek(infp,offset,0);

while ( (ch=fgetc(infp)) != EOF )
  {if (ch=='T' || ch=='U' || ch=='t' || ch=='u' )  {dna[i]= 0; ++i; continue;}
   if (ch=='C' || ch=='c' )   {dna[i]= 1; ++i; continue;}
   if (ch=='A' || ch=='a' )   {dna[i]= 2; ++i; continue;}
   if (ch=='G' || ch=='g' )   {dna[i]= 3; ++i; continue;}
  }

numbp= i--;

if (pstyle%2==0)   fprintf(outfp,"\nnumber of nucleotides: %4d\n", numbp );

return (numbp);

} /* end getindn() */
