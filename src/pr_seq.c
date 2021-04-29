/* PR_SEQ.C;                                    Last update: January 1, 2011. */
/*   - a subroutine to print a sequence array.                                */
/* Dependencies:   called by codcmp.c find_orfs.c prostat.c saps.c ppat.c     */
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


pr_seq(fp,seq,length,ABC,rflag,ftflag)   /* prints seq[] in ABC translation */
FILE *fp; int *seq, length, rflag, ftflag; char ABC[];
{
int i;

for (i=0;i<length;++i)
 {if (ftflag==1 && i%10==0 && i>0)   fprintf(fp," ");
  if (ftflag==1 && i%60==0)   fprintf(fp,"\n  %8d  ",i+1);
  if (ftflag==2 && i%80==0 && i>0)   fprintf(fp,"\n");
  if (rflag) {	/* print reverse complementary sequence; valid only for
		   nucleotide sequences */
    if (seq[length-1-i]<=3) fprintf(fp,"%c", ABC[(seq[length-1-i]+2)%4]);
    else                    fprintf(fp,"%c", ABC[10]);
  } else
    fprintf(fp,"%c", ABC[seq[i]] );
 }
fprintf(fp,"\n");

} /* end pr_seq() */
