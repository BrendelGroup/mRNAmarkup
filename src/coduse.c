/* CODUSE.C;                                  Last update: December 19, 2010. */
/*   - a subroutine to display codon usage in a coding region.                */
/* Dependencies:   called by dnatopro.c, genestat.c, sgap.c                   */
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
extern char AAUC[25], NAUC[12];
extern int codtoaa[4][4][4];

#define NOT_SYNCODUSAGE


coduse(outfp,cseq,glgth,fglgth,Tflag)
FILE *outfp;
int *cseq, glgth, fglgth, Tflag;
			/* Tflag==0: print codon usage numbers for cseq
			   Tflag==1: only tally, do not print codon usage
				      numbers for cseq
			   Tflag==2: print tally of codon usage (Tcf etc.)
			*/
{
int i, j, k;
int cf[4][4][4], df[28][3], mf[11][3], rf[20], aaf[4][4][4];
static int Tglgth;
static int Tcf[4][4][4], Tdf[28][3], Tmf[11][3], Trf[20], Taaf[4][4][4];

if (Tflag!=2)
 {for (i=0;i< 4;++i) for (j=0;j<4;++j) for (k=0;k<4;++k)   cf[i][j][k]= 0;
  for (i=0;i<28;++i) for (j=0;j<3;++j)   df[i][j]= 0;
  for (i=0;i<11;++i) for (j=0;j<3;++j)   mf[i][j]= 0;

  for (i=0;i<=glgth-3;i+=3)
   {++mf[cseq[i]][0]; ++mf[cseq[i+1]][1]; ++mf[cseq[i+2]][2];
    if (cseq[i]!=10 && cseq[i+1]!=10 && cseq[i+2]!=10)
      ++cf[cseq[i]][cseq[i+1]][cseq[i+2]];
    if (cseq[i]!=10 && cseq[i+1]!=10) ++df[cseq[i]*4 + cseq[i+1]][0];
    if (cseq[i+1]!=10 && cseq[i+2]!=10) ++df[cseq[i+1]*4 + cseq[i+2]][1];
    if (i==glgth-3)   break;
    if (cseq[i+2]!=10 && cseq[i+3]!=10) ++df[cseq[i+2]*4 + cseq[i+3]][2];
   }
  if (fglgth>glgth &&
      cseq[fglgth-3]!=10 && cseq[fglgth-2]!=10 && cseq[fglgth-1]!=10)
    ++cf[cseq[fglgth-3]][cseq[fglgth-2]][cseq[fglgth-1]];

  rf[0]= cf[0][0][2]+cf[0][0][3]+cf[1][0][0]+cf[1][0][1]+cf[1][0][2]+
	cf[1][0][3]; Trf[0]+= rf[0];
  rf[1]= cf[3][1][0]+cf[3][1][1]+cf[3][1][2]+cf[3][1][3]; Trf[1]+= rf[1];
  rf[2]= cf[3][3][0]+cf[3][3][1]+cf[3][3][2]+cf[3][3][3]; Trf[2]+= rf[2];
  rf[3]= cf[0][1][0]+cf[0][1][1]+cf[0][1][2]+cf[0][1][3]+cf[2][3][0]+
	cf[2][3][1]; Trf[3]+= rf[3];
  rf[4]= cf[3][0][0]+cf[3][0][1]+cf[3][0][2]+cf[3][0][3]; Trf[4]+= rf[4];
  rf[5]= cf[2][2][2]+cf[2][2][3]; Trf[5]+= rf[5];
  rf[6]= cf[3][2][2]+cf[3][2][3]; Trf[6]+= rf[6];
  rf[7]= cf[2][1][0]+cf[2][1][1]+cf[2][1][2]+cf[2][1][3]; Trf[7]+= rf[7];
  rf[8]= cf[3][2][0]+cf[3][2][1]; Trf[8]+= rf[8];
  rf[9]= cf[2][0][0]+cf[2][0][1]+cf[2][0][2]; Trf[9]+= rf[9];
  rf[10]= cf[1][3][0]+cf[1][3][1]+cf[1][3][2]+cf[1][3][3]+cf[2][3][2]+
	cf[2][3][3]; Trf[10]+= rf[10];
  rf[11]= cf[1][1][0]+cf[1][1][1]+cf[1][1][2]+cf[1][1][3]; Trf[11]+= rf[11];
  rf[12]= cf[2][2][0]+cf[2][2][1]; Trf[12]+= rf[12];
  rf[13]= cf[0][0][0]+cf[0][0][1]; Trf[13]+= rf[13];
  rf[14]= cf[1][2][2]+cf[1][2][3]; Trf[14]+= rf[14];
  rf[15]= cf[0][2][0]+cf[0][2][1]; Trf[15]+= rf[15];
  rf[16]= cf[1][2][0]+cf[1][2][1]; Trf[16]+= rf[16];
  rf[17]= cf[2][0][3]; Trf[17]+= rf[17];
  rf[18]= cf[0][3][0]+cf[0][3][1]; Trf[18]+= rf[18];
  rf[19]= cf[0][3][3]; Trf[19]+= rf[19];
  
  Tglgth+= glgth;
  for (i=0;i< 4;++i) for (j=0;j<4;++j) for (k=0;k<4;++k)
	Tcf[i][j][k]+= cf[i][j][k];
  for (i=0;i<28;++i) for (j=0;j<3;++j)   Tdf[i][j]+= df[i][j];
  for (i=0;i<11;++i) for (j=0;j<3;++j)   Tmf[i][j]+= mf[i][j];
  
  for (i=0;i< 4;++i) for (j=0;j<4;++j) for (k=0;k<4;++k)
    aaf[i][j][k]= rf[codtoaa[i][j][k]];
  aaf[0][2][2]= aaf[0][2][3]= aaf[0][3][2]= 0;

#ifdef NOT_SYNCODUSAGE
  aaf[0][0][2]= aaf[0][0][3]= cf[0][0][2] + cf[0][0][3];
  aaf[1][0][0]= aaf[1][0][1]= aaf[1][0][2]= aaf[1][0][3]= rf[0] - aaf[0][0][2];
  aaf[2][3][0]= aaf[2][3][1]= cf[2][3][0] + cf[2][3][1];
  aaf[0][1][0]= aaf[0][1][1]= aaf[0][1][2]= aaf[0][1][3]= rf[3] - aaf[2][3][0];
  aaf[2][3][2]= aaf[2][3][3]= cf[2][3][2] + cf[2][3][3];
  aaf[1][3][0]= aaf[1][3][1]= aaf[1][3][2]= aaf[1][3][3]= rf[10] - aaf[2][3][2];
#endif
  
  for (i=0;i< 4;++i) for (j=0;j<4;++j) for (k=0;k<4;++k)
	Taaf[i][j][k]+= aaf[i][j][k];
 }

#ifdef SYNCODUSAGE
if (Tflag==0) print_syncodonusage(outfp,rf,cf,aaf,glgth);
if (Tflag==2) print_syncodonusage(outfp,Trf,Tcf,Taaf,Tglgth);
return;
#endif

if (Tflag==0) print_codonusage(outfp,rf,cf,aaf,glgth);
if (Tflag==2) print_codonusage(outfp,Trf,Tcf,Taaf,Tglgth);

if (Tflag!=2)
 {for (i=0;i<3;++i)
   {mf[4][i]= mf[0][i] + mf[1][i];   Tmf[4][i]+= mf[4][i];
    mf[5][i]= mf[2][i] + mf[3][i];   Tmf[5][i]+= mf[5][i];
    mf[6][i]= mf[1][i] + mf[3][i];   Tmf[6][i]+= mf[6][i];
    mf[7][i]= mf[0][i] + mf[2][i];   Tmf[7][i]+= mf[7][i];
    mf[8][i]= mf[0][i] + mf[3][i];   Tmf[8][i]+= mf[8][i];
    mf[9][i]= mf[1][i] + mf[2][i];   Tmf[9][i]+= mf[9][i];
   }
 }

if (Tflag==0) print_mononusage(outfp,mf,glgth);
if (Tflag==2) print_mononusage(outfp,Tmf,Tglgth);

if (Tflag!=2)
 {for (i=0;i<3;++i)
   {df[16][i]=  df[0][i] + df[1][i] + df[4][i] + df[5][i];
	Tdf[16][i]+= df[16][i];
    df[17][i]=  df[2][i] + df[3][i] + df[6][i] + df[7][i];
	Tdf[17][i]+= df[17][i];
    df[18][i]=  df[8][i] + df[9][i] + df[12][i] + df[13][i];
	Tdf[18][i]+= df[18][i];
    df[19][i]=  df[10][i] + df[11][i] + df[14][i] + df[15][i];
	Tdf[19][i]+= df[19][i];
    df[20][i]=  df[5][i] + df[7][i] + df[13][i] + df[15][i];
	Tdf[20][i]+= df[20][i];
    df[21][i]=  df[4][i] + df[6][i] + df[12][i] + df[14][i];
	Tdf[21][i]+= df[21][i];
    df[22][i]=  df[1][i] + df[3][i] + df[9][i] + df[11][i];
	Tdf[22][i]+= df[22][i];
    df[23][i]=  df[0][i] + df[2][i] + df[8][i] + df[10][i];
	Tdf[23][i]+= df[23][i];
   df[24][i]=  df[0][i] + df[3][i] + df[12][i] + df[15][i];
	Tdf[24][i]+= df[24][i];
   df[25][i]=  df[1][i] + df[2][i] + df[13][i] + df[14][i];
	Tdf[25][i]+= df[25][i];
   df[26][i]=  df[4][i] + df[7][i] + df[8][i] + df[11][i];
	Tdf[26][i]+= df[26][i];
   df[27][i]=  df[5][i] + df[6][i] + df[9][i] + df[10][i];
	Tdf[27][i]+= df[27][i];
  }
 }

if (Tflag==0) print_dinusage(outfp,df,glgth);
if (Tflag==2) print_dinusage(outfp,Tdf,Tglgth);

if (Tflag==0) print_resusage(outfp,rf,cf,glgth);
if (Tflag==2) print_resusage(outfp,Trf,Tcf,Tglgth);

} /* end coduse() */



print_codonusage(outfp,rf,cf,aaf,glgth)
FILE *outfp;
int rf[20], cf[4][4][4], aaf[4][4][4], glgth;
{
int i,j,k;

fprintf(outfp,"\nCODON USAGE\n");
fprintf(outfp,"\n   amino acid   codon count      codon freq. / gene\n");
fprintf(outfp,"            (syn. codon c.) (c. freq./ syn. codons)\n\n");

fprintf(outfp,"\n  *        T       *        C       *        A       *");
fprintf(outfp,"        G       *  \n");
fprintf(outfp,"*************************************************");
fprintf(outfp,"*****************************\n");

for (i=0; i<=3; ++i)
 {for (k=0; k<=3; ++k)
   {fprintf(outfp,"%c *  ", NAUC[i]);
    for (j=0; j<=3; ++j)
     {fprintf(outfp,"%c %3d  ", AAUC[codtoaa[i][j][k]], cf[i][j][k]);
      if (aaf[i][j][k]==0 || AAUC[codtoaa[i][j][k]]=='-')
		fprintf(outfp,"       |  ");
      else
		fprintf(outfp,"%5.1f%% |  ",300.*(float)cf[i][j][k]/(float)glgth);
     }
    fprintf(outfp,"  %c\n  *", NAUC[k]);
    for (j=0; j<=3; ++j)
     {if (aaf[i][j][k]==0 || AAUC[codtoaa[i][j][k]]=='-')
		fprintf(outfp,"                |");
      else   fprintf(outfp,"  (%4d)(%5.1f%%)|", aaf[i][j][k],
				100.*(float)cf[i][j][k]/(float)aaf[i][j][k]);
     }
    fprintf(outfp,"\n");
   }
  fprintf(outfp,"\n");
 }

fprintf(outfp,"\nLEU-2 :  %4d of %4d, or %5.1f%%", aaf[0][0][2], rf[0],
		100.*(float)aaf[0][0][2]/(float)rf[0]);
fprintf(outfp,"\nLEU-4 :  %4d of %4d, or %5.1f%%", aaf[1][0][0], rf[0],
		100.*(float)aaf[1][0][0]/(float)rf[0]);
fprintf(outfp,"\n\nARG-2 :  %4d of %4d, or %5.1f%%", aaf[2][3][2], rf[10],
		100.*(float)aaf[2][3][2]/(float)rf[10]);
fprintf(outfp,"\nARG-4 :  %4d of %4d, or %5.1f%%", aaf[1][3][0], rf[10],
		100.*(float)aaf[1][3][0]/(float)rf[10]);
fprintf(outfp,"\n\nSER-2 :  %4d of %4d, or %5.1f%%", aaf[2][3][0], rf[3],
		100.*(float)aaf[2][3][0]/(float)rf[3]);
fprintf(outfp,"\nSER-4 :  %4d of %4d, or %5.1f%%\n", aaf[0][1][0], rf[3],
		100.*(float)aaf[0][1][0]/(float)rf[3]);

} /* end print_codonusage() */



print_mononusage(outfp,mf,glgth)
FILE *outfp;
int mf[11][3], glgth;
{
int i;

fprintf(outfp,"\n\n(DEGENERATE) MONONUCLEOTIDE DISTRIBUTION\n\n");
fprintf(outfp,"   count   frequency per column group\n");
fprintf(outfp,"             (frequency per row group)\n\n");
fprintf(outfp,"  *  position 1  *  position 2  *  position 3  *");
fprintf(outfp,"      all      * \n");
fprintf(outfp,"*****************************************************");
fprintf(outfp,"**************\n");

for (i=0; i<=9; ++i)
 {fprintf(outfp,"%c * %5d  %4.1f%% | %5d  %4.1f%% | %5d  %4.1f%% |",
            NAUC[i],
	    mf[i][0], 300.*(float)mf[i][0]/(float)glgth,
	    mf[i][1], 300.*(float)mf[i][1]/(float)glgth,
	    mf[i][2], 300.*(float)mf[i][2]/(float)glgth);
  fprintf(outfp,"  %5d  %4.1f%% * \n",
	    mf[i][0]+mf[i][1]+mf[i][2], 
	    100.*(float)(mf[i][0]+mf[i][1]+mf[i][2])/(float)glgth); 
  fprintf(outfp,"  *       (%4.1f%%)|       (%4.1f%%)|       (%4.1f%%)|",
	    100.*(float)mf[i][0]/(float)(mf[i][0]+mf[i][1]+mf[i][2]),
	    100.*(float)mf[i][1]/(float)(mf[i][0]+mf[i][1]+mf[i][2]),
	    100.*(float)mf[i][2]/(float)(mf[i][0]+mf[i][1]+mf[i][2]));
  fprintf(outfp,"               * \n");
  if (i==3 || i==5 || i==7 || i==9)
   {fprintf(outfp,"*****************************************************");
    fprintf(outfp,"**************\n");
   }
 }
if (mf[10][0]+mf[10][1]+mf[10][2])
 {fprintf(outfp,"%c * %5d        | %5d        | %5d        |  %5d        * \n",
		NAUC[10],mf[10][0],mf[10][1],mf[10][2],
		mf[10][0]+mf[10][1]+mf[10][2]          );
  fprintf(outfp,"*****************************************************");
  fprintf(outfp,"**************\n");
 }

} /* end print_mononusage() */



print_dinusage(outfp,df,glgth)
FILE *outfp;
int df[28][3], glgth;
{
int i;

fprintf(outfp,"\n\nDINUCLEOTIDE DISTRIBUTION\n\n");
fprintf(outfp,"   count   frequency per column group\n");
fprintf(outfp,"             (frequency per row group)\n\n");
fprintf(outfp,"  *  position 1,2*  position 2,3*  position 3,1*");
fprintf(outfp,"      all      * \n");
fprintf(outfp,"*****************************************************");
fprintf(outfp,"**************\n");

for (i=0;i<16;++i)
 {fprintf(outfp,"%c%c* %5d  %4.1f%% | %5d  %4.1f%% | %5d  %4.1f%% |",
            NAUC[i/4], NAUC[i%4],
	    df[i][0], 300.*(float)df[i][0]/(float)glgth,
	    df[i][1], 300.*(float)df[i][1]/(float)glgth,
	    df[i][2], 300.*(float)df[i][2]/(float)(glgth-3));
  fprintf(outfp,"  %5d  %4.1f%% * \n",
	    df[i][0]+df[i][1]+df[i][2], 
	    100.*(float)(df[i][0]+df[i][1]+df[i][2])/(float)(glgth-1)); 
  fprintf(outfp,"  *       (%4.1f%%)|       (%4.1f%%)|       (%4.1f%%)|",
	    100.*(float)df[i][0]/(float)(df[i][0]+df[i][1]+df[i][2]),
	    100.*(float)df[i][1]/(float)(df[i][0]+df[i][1]+df[i][2]),
	    100.*(float)df[i][2]/(float)(df[i][0]+df[i][1]+df[i][2]));
  fprintf(outfp,"               * \n");
 }
fprintf(outfp,"*****************************************************");
fprintf(outfp,"**************\n");


fprintf(outfp,"\n\nDEGENERATE DINUCLEOTIDE DISTRIBUTION\n\n");
fprintf(outfp,"   count   frequency per column group\n");
fprintf(outfp,"             (frequency per row group)\n\n");
fprintf(outfp,"  *  position 1,2*  position 2,3*  position 3,1*");
fprintf(outfp,"      all      * \n");
fprintf(outfp,"*****************************************************");
fprintf(outfp,"**************\n");

for (i=16;i<20;++i)
 {fprintf(outfp,"%c%c* %5d  %4.1f%% | %5d  %4.1f%% | %5d  %4.1f%% |",
            NAUC[(i-16)/2 +4], NAUC[i-12 -2*((i-16)/2)],
	    df[i][0], 300.*(float)df[i][0]/(float)glgth,
	    df[i][1], 300.*(float)df[i][1]/(float)glgth,
	    df[i][2], 300.*(float)df[i][2]/(float)(glgth-3));
  fprintf(outfp,"  %5d  %4.1f%% * \n",
	    df[i][0]+df[i][1]+df[i][2], 
	    100.*(float)(df[i][0]+df[i][1]+df[i][2])/(float)(glgth-1)); 
  fprintf(outfp,"  *       (%4.1f%%)|       (%4.1f%%)|       (%4.1f%%)|",
	    100.*(float)df[i][0]/(float)(df[i][0]+df[i][1]+df[i][2]),
	    100.*(float)df[i][1]/(float)(df[i][0]+df[i][1]+df[i][2]),
	    100.*(float)df[i][2]/(float)(df[i][0]+df[i][1]+df[i][2]));
  fprintf(outfp,"               * \n");
 }
fprintf(outfp,"*****************************************************");
fprintf(outfp,"**************\n");

for (i=20;i<24;++i)
 {fprintf(outfp,"%c%c* %5d  %4.1f%% | %5d  %4.1f%% | %5d  %4.1f%% |",
            NAUC[(i-20)/2 +6], NAUC[i-14 -2*((i-20)/2)],
	    df[i][0], 300.*(float)df[i][0]/(float)glgth,
	    df[i][1], 300.*(float)df[i][1]/(float)glgth,
	    df[i][2], 300.*(float)df[i][2]/(float)(glgth-3));
  fprintf(outfp,"  %5d  %4.1f%% * \n",
	    df[i][0]+df[i][1]+df[i][2], 
	    100.*(float)(df[i][0]+df[i][1]+df[i][2])/(float)(glgth-1)); 
  fprintf(outfp,"  *       (%4.1f%%)|       (%4.1f%%)|       (%4.1f%%)|",
	    100.*(float)df[i][0]/(float)(df[i][0]+df[i][1]+df[i][2]),
	    100.*(float)df[i][1]/(float)(df[i][0]+df[i][1]+df[i][2]),
	    100.*(float)df[i][2]/(float)(df[i][0]+df[i][1]+df[i][2]));
  fprintf(outfp,"               * \n");
  fprintf(outfp,"  *              |              |              |");
  fprintf(outfp,"               * \n");
 }
fprintf(outfp,"*****************************************************");
fprintf(outfp,"**************\n");

for (i=24;i<28;++i)
 {fprintf(outfp,"%c%c* %5d  %4.1f%% | %5d  %4.1f%% | %5d  %4.1f%% |",
            NAUC[(i-24)/2 +8], NAUC[i-16 -2*((i-24)/2)],
	    df[i][0], 300.*(float)df[i][0]/(float)glgth,
	    df[i][1], 300.*(float)df[i][1]/(float)glgth,
	    df[i][2], 300.*(float)df[i][2]/(float)(glgth-3));
  fprintf(outfp,"  %5d  %4.1f%% * \n",
	    df[i][0]+df[i][1]+df[i][2], 
	    100.*(float)(df[i][0]+df[i][1]+df[i][2])/(float)(glgth-1)); 
  fprintf(outfp,"  *       (%4.1f%%)|       (%4.1f%%)|       (%4.1f%%)|",
	    100.*(float)df[i][0]/(float)(df[i][0]+df[i][1]+df[i][2]),
	    100.*(float)df[i][1]/(float)(df[i][0]+df[i][1]+df[i][2]),
	    100.*(float)df[i][2]/(float)(df[i][0]+df[i][1]+df[i][2]));
  fprintf(outfp,"               * \n");
  fprintf(outfp,"  *              |              |              |");
  fprintf(outfp,"               * \n");
 }
fprintf(outfp,"*****************************************************");
fprintf(outfp,"**************\n");

} /* end print_dinusage() */



print_resusage(outfp,rf,cf,glgth)
FILE *outfp;
int rf[23], cf[4][4][4], glgth;
{

fprintf(outfp,"\n\nAMINO ACID USAGE\n\n");
fprintf(outfp,"amino acid, count, frequency/gene:  ");
fprintf(outfp,"(codon, count, frequency/amino acid)\n\n");

fprintf(outfp,"A=Ala %5d", rf[1]);
if (rf[1])
 {fprintf(outfp," %4.1f%%: GCU %3d %5.1f%% GCC %3d %5.1f%% ",
		300.*(float)rf[1]/(float)glgth, 
		cf[3][1][0], 100.*(float)cf[3][1][0]/(float)rf[1],
		cf[3][1][1], 100.*(float)cf[3][1][1]/(float)rf[1]);
  fprintf(outfp,"GCA %3d %5.1f%% GCG %3d %5.1f%%",
		cf[3][1][2], 100.*(float)cf[3][1][2]/(float)rf[1],
		cf[3][1][3], 100.*(float)cf[3][1][3]/(float)rf[1]);
 }
fprintf(outfp,"\n");
fprintf(outfp,"C=Cys %5d", rf[18]);
if (rf[18])
 {fprintf(outfp," %4.1f%%: UGU %3d %5.1f%% UGC %3d %5.1f%% ",
		300.*(float)rf[18]/(float)glgth, 
		cf[0][3][0], 100.*(float)cf[0][3][0]/(float)rf[18],
		cf[0][3][1], 100.*(float)cf[0][3][1]/(float)rf[18]);
 }
fprintf(outfp,"\n");
fprintf(outfp,"D=Asp %5d", rf[8]);
if (rf[8])
 {fprintf(outfp," %4.1f%%: GAU %3d %5.1f%% GAC %3d %5.1f%% ",
		300.*(float)rf[8]/(float)glgth, 
		cf[3][2][0], 100.*(float)cf[3][2][0]/(float)rf[8],
		cf[3][2][1], 100.*(float)cf[3][2][1]/(float)rf[8]);
 }
fprintf(outfp,"\n");
fprintf(outfp,"E=Glu %5d", rf[6]);
if (rf[6])
 {fprintf(outfp," %4.1f%%: GAA %3d %5.1f%% GAG %3d %5.1f%% ",
		300.*(float)rf[6]/(float)glgth, 
		cf[3][2][2], 100.*(float)cf[3][2][2]/(float)rf[6],
		cf[3][2][3], 100.*(float)cf[3][2][3]/(float)rf[6]);
 }
fprintf(outfp,"\n");
fprintf(outfp,"F=Phe %5d", rf[13]);
if (rf[13])
 {fprintf(outfp," %4.1f%%: UUU %3d %5.1f%% UUC %3d %5.1f%%",
		300.*(float)rf[13]/(float)glgth, 
		cf[0][0][0], 100.*(float)cf[0][0][0]/(float)rf[13],
		cf[0][0][1], 100.*(float)cf[0][0][1]/(float)rf[13]);
 }
fprintf(outfp,"\n");
fprintf(outfp,"G=Gly %5d", rf[2]);
if (rf[2])
 {fprintf(outfp," %4.1f%%: GGU %3d %5.1f%% GGC %3d %5.1f%% ",
		300.*(float)rf[2]/(float)glgth, 
		cf[3][3][0], 100.*(float)cf[3][3][0]/(float)rf[2],
		cf[3][3][1], 100.*(float)cf[3][3][1]/(float)rf[2]);
  fprintf(outfp,"GGA %3d %5.1f%% GGG %3d %5.1f%%",
		cf[3][3][2], 100.*(float)cf[3][3][2]/(float)rf[2],
		cf[3][3][3], 100.*(float)cf[3][3][3]/(float)rf[2]);
 }
fprintf(outfp,"\n");
fprintf(outfp,"H=His %5d", rf[16]);
if (rf[16])
 {fprintf(outfp," %4.1f%%: CAU %3d %5.1f%% CAC %3d %5.1f%%",
		300.*(float)rf[16]/(float)glgth, 
		cf[1][2][0], 100.*(float)cf[1][2][0]/(float)rf[16],
		cf[1][2][1], 100.*(float)cf[1][2][1]/(float)rf[16]);
 }
fprintf(outfp,"\n");
fprintf(outfp,"I=Ile %5d", rf[9]);
if (rf[9])
 {fprintf(outfp," %4.1f%%: AUU %3d %5.1f%% AUC %3d %5.1f%% ",
		300.*(float)rf[9]/(float)glgth, 
		cf[2][0][0], 100.*(float)cf[2][0][0]/(float)rf[9],
		cf[2][0][1], 100.*(float)cf[2][0][1]/(float)rf[9]);
  fprintf(outfp,"AUA %3d %5.1f%%",
		cf[2][0][2], 100.*(float)cf[2][0][2]/(float)rf[9]);
 }
fprintf(outfp,"\n");
fprintf(outfp,"K=Lys %5d", rf[5]);
if (rf[5])
 {fprintf(outfp," %4.1f%%: AAA %3d %5.1f%% AAG %3d %5.1f%%",
		300.*(float)rf[5]/(float)glgth, 
		cf[2][2][2], 100.*(float)cf[2][2][2]/(float)rf[5],
		cf[2][2][3], 100.*(float)cf[2][2][3]/(float)rf[5]);
 }
fprintf(outfp,"\n");
fprintf(outfp,"L=Leu %5d", rf[0]);
if (rf[0])
 {fprintf(outfp," %4.1f%%: UUA %3d %5.1f%% UUG %3d %5.1f%%\n",
		300.*(float)rf[0]/(float)glgth, 
		cf[0][0][2], 100.*(float)cf[0][0][2]/(float)rf[0],
		cf[0][0][3], 100.*(float)cf[0][0][3]/(float)rf[0]);
  fprintf(outfp,"                   CUU %3d %5.1f%% CUC %3d %5.1f%% ",
		cf[1][0][0], 100.*(float)cf[1][0][0]/(float)rf[0],
		cf[1][0][1], 100.*(float)cf[1][0][1]/(float)rf[0]);
  fprintf(outfp,"CUA %3d %5.1f%% CUG %3d %5.1f%%",
		cf[1][0][2], 100.*(float)cf[1][0][2]/(float)rf[0],
		cf[1][0][3], 100.*(float)cf[1][0][3]/(float)rf[0]);
 }
fprintf(outfp,"\n");
fprintf(outfp,"M=Met %5d", rf[17]);
if (rf[17])
 {fprintf(outfp," %4.1f%%: AUG %3d",
		300.*(float)rf[17]/(float)glgth, cf[2][0][3]);
 }
fprintf(outfp,"\n");
fprintf(outfp,"N=Asn %5d", rf[12]);
if (rf[12])
 {fprintf(outfp," %4.1f%%: AAU %3d %5.1f%% AAC %3d %5.1f%%",
		300.*(float)rf[12]/(float)glgth, 
		cf[2][2][0], 100.*(float)cf[2][2][0]/(float)rf[12],
		cf[2][2][1], 100.*(float)cf[2][2][1]/(float)rf[12]);
 }
fprintf(outfp,"\n");
fprintf(outfp,"P=Pro %5d", rf[11]);
if (rf[11])
 {fprintf(outfp," %4.1f%%: CCU %3d %5.1f%% CCC %3d %5.1f%% ",
		300.*(float)rf[11]/(float)glgth, 
		cf[1][1][0], 100.*(float)cf[1][1][0]/(float)rf[11],
		cf[1][1][1], 100.*(float)cf[1][1][1]/(float)rf[11]);
  fprintf(outfp,"CCA %3d %5.1f%% CCG %3d %5.1f%%",
		cf[1][1][2], 100.*(float)cf[1][1][2]/(float)rf[11],
		cf[1][1][3], 100.*(float)cf[1][1][3]/(float)rf[11]);
 }
fprintf(outfp,"\n");
fprintf(outfp,"Q=Gln %5d", rf[14]);
if (rf[14])
 {fprintf(outfp," %4.1f%%: CAA %3d %5.1f%% CAG %3d %5.1f%%",
		300.*(float)rf[14]/(float)glgth, 
		cf[1][2][2], 100.*(float)cf[1][2][2]/(float)rf[14],
		cf[1][2][3], 100.*(float)cf[1][2][3]/(float)rf[14]);
 }
fprintf(outfp,"\n");
fprintf(outfp,"R=Arg %5d", rf[10]);
if (rf[10])
 {fprintf(outfp," %4.1f%%: CGU %3d %5.1f%% CGC %3d %5.1f%% ",
		300.*(float)rf[10]/(float)glgth, 
		cf[1][3][0], 100.*(float)cf[1][3][0]/(float)rf[10],
		cf[1][3][1], 100.*(float)cf[1][3][1]/(float)rf[10]);
  fprintf(outfp,"CGA %3d %5.1f%% CGG %3d %5.1f%%\n",
		cf[1][3][2], 100.*(float)cf[1][3][2]/(float)rf[10],
		cf[1][3][3], 100.*(float)cf[1][3][3]/(float)rf[10]);
  fprintf(outfp,"                   AGA %3d %5.1f%% AGG %3d %5.1f%%",
		cf[2][3][2], 100.*(float)cf[2][3][2]/(float)rf[10],
		cf[2][3][3], 100.*(float)cf[2][3][3]/(float)rf[10]);
 }
fprintf(outfp,"\n");
fprintf(outfp,"S=Ser %5d", rf[3]);
if (rf[3])
 {fprintf(outfp," %4.1f%%: UCU %3d %5.1f%% UCC %3d %5.1f%% ",
		300.*(float)rf[3]/(float)glgth, 
		cf[0][1][0], 100.*(float)cf[0][1][0]/(float)rf[3],
		cf[0][1][1], 100.*(float)cf[0][1][1]/(float)rf[3]);
  fprintf(outfp,"UCA %3d %5.1f%% UCG %3d %5.1f%%\n",
		cf[0][1][2], 100.*(float)cf[0][1][2]/(float)rf[3],
		cf[0][1][3], 100.*(float)cf[0][1][3]/(float)rf[3]);
  fprintf(outfp,"                   AGU %3d %5.1f%% AGC %3d %5.1f%%",
		cf[2][3][0], 100.*(float)cf[2][3][0]/(float)rf[3],
		cf[2][3][1], 100.*(float)cf[2][3][1]/(float)rf[3]);
 }
fprintf(outfp,"\n");
fprintf(outfp,"T=Thr %5d", rf[7]);
if (rf[7])
 {fprintf(outfp," %4.1f%%: ACU %3d %5.1f%% ACC %3d %5.1f%% ",
		300.*(float)rf[7]/(float)glgth, 
		cf[2][1][0], 100.*(float)cf[2][1][0]/(float)rf[7],
		cf[2][1][1], 100.*(float)cf[2][1][1]/(float)rf[7]);
  fprintf(outfp,"ACA %3d %5.1f%% ACG %3d %5.1f%%",
		cf[2][1][2], 100.*(float)cf[2][1][2]/(float)rf[7],
		cf[2][1][3], 100.*(float)cf[2][1][3]/(float)rf[7]);
 }
fprintf(outfp,"\n");
fprintf(outfp,"V=Val %5d", rf[4]);
if (rf[4])
 {fprintf(outfp," %4.1f%%: GUU %3d %5.1f%% GUC %3d %5.1f%% ",
		300.*(float)rf[4]/(float)glgth, 
		cf[3][0][0], 100.*(float)cf[3][0][0]/(float)rf[4],
		cf[3][0][1], 100.*(float)cf[3][0][1]/(float)rf[4]);
  fprintf(outfp,"GUA %3d %5.1f%% GUG %3d %5.1f%%",
		cf[3][0][2], 100.*(float)cf[3][0][2]/(float)rf[4],
		cf[3][0][3], 100.*(float)cf[3][0][3]/(float)rf[4]);
 }
fprintf(outfp,"\n");
fprintf(outfp,"W=Trp %5d", rf[19]);
if (rf[19])
  fprintf(outfp," %4.1f%%: UGG %3d",
		300.*(float)rf[19]/(float)glgth, cf[0][3][3]);
fprintf(outfp,"\n");
fprintf(outfp,"Y=Tyr %5d", rf[15]);
if (rf[15])
  fprintf(outfp," %4.1f%%: UAU %3d %5.1f%% UAC %3d %5.1f%%",
		300.*(float)rf[15]/(float)glgth, 
		cf[0][2][0], 100.*(float)cf[0][2][0]/(float)rf[15],
		cf[0][2][1], 100.*(float)cf[0][2][1]/(float)rf[15]);
fprintf(outfp,"\n");

} /* end print_resusage() */



print_syncodonusage(outfp,rf,cf,aaf,glgth)
FILE *outfp;
int rf[20], cf[4][4][4], aaf[4][4][4], glgth;
{
int i,j,k;

fprintf(outfp,"\n  *        T       *        C       *        A       *");
fprintf(outfp,"        G       *  \n");
fprintf(outfp,"*************************************************");
fprintf(outfp,"*****************************\n");

for (i=0; i<=3; ++i)
 {for (k=0; k<=3; ++k)
   {fprintf(outfp,"%c *  ", NAUC[i]);
    for (j=0; j<=3; ++j)
     {fprintf(outfp,"%c      ", AAUC[codtoaa[i][j][k]]);
      if (aaf[i][j][k]==0 || AAUC[codtoaa[i][j][k]]=='-')
		fprintf(outfp,"       |  ");
      else
		fprintf(outfp,"%5.1f%% |  ",100.*(float)cf[i][j][k]/(float)aaf[i][j][k]);
     }
    fprintf(outfp,"  %c\n", NAUC[k]);
   }
  fprintf(outfp,"\n");
 }

} /* end print_syncodonusage() */
