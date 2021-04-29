/* RESUSE.C;                                  Last update: December 19, 2010. */
/*   - a subroutine to determine the residue usage of a protein sequence.     */
/* Dependencies:   called by aaruns.c, dnatopro.c, find_orfs.c, genestat.c,   */
/*   paste_exons.c, pgrep.c, prostat.c, pscore.c, ranseq.c, resuse.c, saps.c  */
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

resuse(seq,numaa,af,pf,outfp,pstyle,pflag,aq,pq,qflag,tabfp,tflag,AAUC,CHPN)
int seq[], numaa, af[23], pf[9], pstyle, pflag, qflag, tflag;
float aq[20][4], pq[9][4];
FILE *outfp, *tabfp;
char AAUC[], CHPN[];
  /* qflag= 1: compare composition with quantile points aq[][] and pq[][] */
{
int i;
static char ps[9][7]= {"KR    ","ED    ","KRED  ","'0'   ","KR-ED ","LVIFM ",
	"ST    ","AGP   ","FIKMNY"};
char c[29][2];
int lct= 0, hct= 0;
int eu[8];

if (CHPN[16]=='+')
  {strcpy(ps[0],"KRH   "); strcpy(ps[2],"KRHED "); strcpy(ps[4],"KRH-ED");}

for (i=0; i< 8; ++i)   eu[i]=0;
for (i=0; i<=22; ++i)   af[i]=0;
for (i=0; i<numaa; ++i)
  if (seq[i]<=22) ++af[seq[i]]; /* amino acid frequencies */

numaa= numaa-af[20]-af[21]-af[22]; /* ignore B,Z,X in calculating proportions */

pf[0]= af[5]+af[10]; /* frequency of K+R[+H] */
  if (CHPN[16]=='+')   pf[0]+= af[16];  /* add H to positive charge group */
pf[1]= af[6]+af[8]; /* frequency of E+D */
pf[2]= pf[0]+pf[1]; /* frequency of K+R[+H]+E+D */
pf[3]= numaa-pf[2]; /* frequency of non-(K+R[+H]+E+D) */
pf[4]= pf[0]-pf[1]; /* frequency of K+R[+H]-E-D */
pf[5]= af[0]+af[4]+af[9]+af[13]+af[17]; /* frequency of L+V+I+F+M */
pf[6]= af[3]+af[7]; /* frequency of S+T */
pf[7]= af[1]+af[2]+af[11]; /* frequency of A+G+P */
pf[8]= af[13]+af[9]+af[5]+af[17]+af[12]+af[15]; /* frequency of F+I+K+M+N+Y */

/* If qflag= 1 and numaa>=200, observed proportions are compared with the
   quantile points aq[][], pq[][]; upper extremes are labeled
   ++ (.99 point) or + (.95 point), lower extremes - (.05 point) or
   -- (.01 point):
*/
for (i=0;i<=19;++i)
  {c[i][0]=' ';   c[i][1]=' ';
   if (numaa<200 || qflag==0)   continue;
   if (100.*(float)af[i]/(float)numaa < aq[i][1])   
     {c[i][0]='-';   ++eu[1];
      if (100.*(float)af[i]/(float)numaa < aq[i][0])
	{c[i][1]='-';   --eu[1]; ++eu[0];}
      if (tflag)
	{++lct;
	 if (lct+hct==1)   fprintf(tabfp,"RE ");
	 else if ((lct+hct)%7==0)   fprintf(tabfp,"\nRE ");
	 fprintf(tabfp,"%c%c%c: %4.1f; ",
		AAUC[i], c[i][0], c[i][1], 100.*(float)af[i]/(float)numaa );
	}
      continue;
     }
   if (100.*(float)af[i]/(float)numaa > aq[i][2])
     {c[i][0]='+';   ++eu[2];
      if (100.*(float)af[i]/(float)numaa > aq[i][3])
	{c[i][1]='+';   --eu[2]; ++eu[3];}
      if (tflag)
	{++hct;
	 if (lct+hct==1)   fprintf(tabfp,"RE ");
	 else if ((lct+hct)%7==0)   fprintf(tabfp,"\nRE ");
	 fprintf(tabfp,"%c%c%c: %4.1f; ",
		AAUC[i], c[i][0], c[i][1], 100.*(float)af[i]/(float)numaa );
	}
      continue;
     }
  }
if (tflag && lct+hct>0)   fprintf(tabfp,"\n");
lct= hct= 0;
for (i=0;i<9;++i)  
  {c[20+i][0]=' ';   c[20+i][1]=' ';
   if (numaa<200 || qflag==0)   continue;
   if (CHPN[16]=='+')   continue; /* quantile points valid for += K,R only */
   if (100.*(float)pf[i]/(float)numaa < pq[i][1])   
     {c[20+i][0]='-';   ++eu[5];
      if (100.*(float)pf[i]/(float)numaa < pq[i][0])
	{c[20+i][1]='-';   --eu[5]; ++eu[4];}
      if (tflag)
	{++lct;
	 if (lct+hct==1)   fprintf(tabfp,"PE ");
	 else if ((lct+hct)%5==0)   fprintf(tabfp,"\nPE ");
	 fprintf(tabfp,"%s%c%c: %4.1f; ",
		ps[i], c[20+i][0], c[20+i][1], 100.*(float)pf[i]/(float)numaa );
	}
      continue;
     }
   if (100.*(float)pf[i]/(float)numaa > pq[i][2])
     {c[20+i][0]='+';   ++eu[6];
      if (100.*(float)pf[i]/(float)numaa > pq[i][3])
	{c[20+i][1]='+';   --eu[6]; ++eu[7];}
      if (tflag)
	{++hct;
	 if (lct+hct==1)   fprintf(tabfp,"PE ");
	 else if ((lct+hct)%5==0)   fprintf(tabfp,"\nPE ");
	 fprintf(tabfp,"%s%c%c: %4.1f; ",
		ps[i], c[20+i][0], c[20+i][1], 100.*(float)pf[i]/(float)numaa );
	}
      continue;
     }
  }
if (tflag && lct+hct>0)   fprintf(tabfp,"\n");
if (tflag)   fprintf(tabfp,"EU %2d %2d %2d %2d     %2d %2d %2d %2d\n",
		       eu[0], eu[1], eu[2], eu[3],  eu[4], eu[5], eu[6], eu[7]);

if (pflag && pstyle%2==0)
 {fprintf(outfp,"\n");
  fprintf(outfp,"A%c%c:%3d(%4.1f%%); C%c%c:%3d(%4.1f%%); D%c%c:%3d(%4.1f%%);",
	c[1][0], c[1][1], af[1], 100.*(float)af[1]/(float)numaa,
	c[18][0], c[18][1], af[18], 100.*(float)af[18]/(float)numaa,
	c[8][0], c[8][1], af[8], 100.*(float)af[8]/(float)numaa );
  fprintf(outfp," E%c%c:%3d(%4.1f%%); F%c%c:%3d(%4.1f%%)\n",
	c[6][0], c[6][1], af[6], 100.*(float)af[6]/(float)numaa,
	c[13][0], c[13][1], af[13], 100.*(float)af[13]/(float)numaa );
  fprintf(outfp,"G%c%c:%3d(%4.1f%%); H%c%c:%3d(%4.1f%%); I%c%c:%3d(%4.1f%%);",
	c[2][0], c[2][1], af[2], 100.*(float)af[2]/(float)numaa,
	c[16][0], c[16][1], af[16], 100.*(float)af[16]/(float)numaa,
	c[9][0], c[9][1], af[9], 100.*(float)af[9]/(float)numaa );
  fprintf(outfp," K%c%c:%3d(%4.1f%%); L%c%c:%3d(%4.1f%%)\n",
	c[5][0], c[5][1], af[5], 100.*(float)af[5]/(float)numaa,
	c[0][0], c[0][1], af[0], 100.*(float)af[0]/(float)numaa );
  fprintf(outfp,"M%c%c:%3d(%4.1f%%); N%c%c:%3d(%4.1f%%); P%c%c:%3d(%4.1f%%);",
	c[17][0], c[17][1], af[17], 100.*(float)af[17]/(float)numaa,
	c[12][0], c[12][1], af[12], 100.*(float)af[12]/(float)numaa,
	c[11][0], c[11][1], af[11], 100.*(float)af[11]/(float)numaa );
  fprintf(outfp," Q%c%c:%3d(%4.1f%%); R%c%c:%3d(%4.1f%%)\n",
	c[14][0], c[14][1], af[14], 100.*(float)af[14]/(float)numaa,
	c[10][0], c[10][1], af[10], 100.*(float)af[10]/(float)numaa );
  fprintf(outfp,"S%c%c:%3d(%4.1f%%); T%c%c:%3d(%4.1f%%); V%c%c:%3d(%4.1f%%);",
	c[3][0], c[3][1], af[3], 100.*(float)af[3]/(float)numaa,
	c[7][0], c[7][1], af[7], 100.*(float)af[7]/(float)numaa,
	c[4][0], c[4][1], af[4], 100.*(float)af[4]/(float)numaa );
  fprintf(outfp," W%c%c:%3d(%4.1f%%); Y%c%c:%3d(%4.1f%%)\n",
	c[19][0], c[19][1], af[19], 100.*(float)af[19]/(float)numaa,
	c[15][0], c[15][1], af[15], 100.*(float)af[15]/(float)numaa );

  if (af[20]+af[21]+af[22])   fprintf(outfp,
	"B: %3d; Z: %3d; X: %3d   (ignored in calculating proportions)\n",
	af[20], af[21], af[22]);
  fprintf(outfp,"\n");

  fprintf(outfp,
  "%s%c%c: %4d (%5.1f%%);   %s%c%c: %4d (%5.1f%%);   %s%c%c: %4d (%5.1f%%);\n",
	ps[0], c[20][0], c[20][1], pf[0], 100.*(float)pf[0]/(float)numaa,
	ps[1], c[21][0], c[21][1], pf[1], 100.*(float)pf[1]/(float)numaa,
	ps[7], c[27][0], c[27][1], pf[7], 100.*(float)pf[7]/(float)numaa);
  fprintf(outfp,
  "%s%c%c: %4d (%5.1f%%);   %s%c%c: %4d (%5.1f%%);   %s%c%c: %4d (%5.1f%%);\n",
	ps[2], c[22][0], c[22][1], pf[2], 100.*(float)pf[2]/(float)numaa,
	ps[4], c[24][0], c[24][1], pf[4], 100.*(float)pf[4]/(float)numaa,
	ps[8], c[28][0], c[28][1], pf[8], 100.*(float)pf[8]/(float)numaa );
  fprintf(outfp,
  "%s%c%c: %4d (%5.1f%%);   %s%c%c: %4d (%5.1f%%).\n",
	ps[5], c[25][0], c[25][1], pf[5], 100.*(float)pf[5]/(float)numaa,
	ps[6], c[26][0], c[26][1], pf[6], 100.*(float)pf[6]/(float)numaa );
  fprintf(outfp,"\n");
 }

} /* end resuse() */
