/* PR_PRO.C;                                  Last update: December 19, 2010. */
/*   - a subroutine to print a protein sequence in a specified alphabet.      */
/* Dependencies:   called by aaruns.c, dnatopro.c, find_orfs.c, genestat.c,   */
/*   mproc.c, paste_exons.c, pdbtoind.c, pmagic.c, pmotif.c, propeat.c,       */
/*   prorep.c, prostat.c, pscore.c, ranseq.c, saps.c                          */
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

#include "def.h"
extern FILE *outfp;
extern protein[PROTLGTH];


pr_pro(numaa,ABC)   /* prints protein[] in ABC translation */
int numaa;
char ABC[];
{
int i;

for( i=0; i<numaa ; ++i )
  {if (i%10 == 0)   fprintf(outfp," ");
   if(i%60== 0)  fprintf(outfp, "\n%8d  ",i+1);
   fprintf(outfp,"%c", ABC[protein[i]] );
  }
fprintf(outfp,"\n");

} /* end pr_pro() */
