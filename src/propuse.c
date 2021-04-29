/* PROPUSE.C;                                 Last update: December 19, 2010. */
/*   - a subroutine to establish property usage for a protein sequence.       */
/* Dependencies:   called by sgap.c                                           */
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
extern int af[23];


propuse(numaa,Tflag)
int numaa, Tflag;
			/* Tflag==0: print property usage numbers
			   Tflag==1: only tally, do not print property usage
				      numbers
			   Tflag==2: print tally of property usage (Tcpf etc.)
			*/
{
int i;
int cpf[8], fpf[4], spf[3], dgf[5];
static int Taf[23], Tcpf[8], Tfpf[4], Tspf[3], Tdgf[5];

if (Tflag!=2)
 {for (i=0;i<20;++i) Taf[i]+= af[i];
  cpf[0] = af[8] + af[6]; Tcpf[0]+= cpf[0];
  cpf[1] = af[16] + af[5] + af[10]; Tcpf[1]+= cpf[1];
  cpf[2] = af[3] + af[7]; Tcpf[2]+= cpf[2];
  cpf[3] = af[11]; Tcpf[3]+= cpf[3];
  cpf[4] = af[1] + af[2] + af[9] + af[0] + af[4]; Tcpf[4]+= cpf[4];
  cpf[5] = af[12] + af[14]; Tcpf[5]+= cpf[5];
  cpf[6] = af[13] + af[19] + af[15]; Tcpf[6]+= cpf[6];
  cpf[7] = af[18] + af[17]; Tcpf[7]+= cpf[7];

  fpf[0] = af[1] + af[13] + af[9] + af[0] + af[17] + af[11] + af[4] + af[19];
	Tfpf[0]+= fpf[0];
  fpf[1] = cpf[0]; Tfpf[1]+= fpf[1];
  fpf[2] = cpf[1]; Tfpf[2]+= fpf[2];
  fpf[3] = af[18] + af[2] + cpf[5] + cpf[2] + af[15]; Tfpf[3]+= fpf[3];

  spf[0] = af[1] + af[18] + af[2] + af[11] + cpf[2] + af[19] + af[15];
	Tspf[0]+= spf[0];
  spf[1] = cpf[0] + cpf[1] + cpf[5]; Tspf[1]+= spf[1];
  spf[2] = af[13] + af[9] + af[0] + af[17] + af[4]; Tspf[2]+= spf[2];

  dgf[0] = af[17] + af[19]; Tdgf[0]+= dgf[0];
  dgf[1] = af[18] + cpf[0] + af[13] + af[16] + af[5] + cpf[5] + af[15]; 
	Tdgf[1]+= dgf[1];
  dgf[2] = af[9]; Tdgf[2]+= dgf[2];
  dgf[3] = af[1] + af[2] + af[11] + af[7] + af[4]; Tdgf[3]+= dgf[3];
  dgf[4] = af[0] + af[10] + af[3]; Tdgf[4]+= dgf[4];
 }

if (Tflag==0) print_propusage(numaa,af,cpf,fpf,spf,dgf);
if (Tflag==2) print_propusage(numaa,Taf,Tcpf,Tfpf,Tspf,Tdgf);

} /* end propuse() */



print_propusage(slgth,rf,cpf,fpf,spf,dgf)
int rf[23], cpf[8], fpf[4], spf[3], dgf[5];
{

fprintf(outfp,"\nCHEMICAL PROPERTY USAGE\n\n" );

fprintf(outfp,"A=acidic     %5d %5.1f%% :  D %5d %5.1f%%  E %5d %5.1f%%\n",
		cpf[0], 100.*(float)cpf[0]/(float)slgth,
		rf[8], 100.*(float)rf[8]/(float)cpf[0],
		rf[6], 100.*(float)rf[6]/(float)cpf[0] );
fprintf(outfp,"B=basic      %5d %5.1f%% :  H %5d %5.1f%%  K %5d %5.1f%%  ",
		cpf[1], 100.*(float)cpf[1]/(float)slgth,
		rf[16], 100.*(float)rf[16]/(float)cpf[1],
		rf[5], 100.*(float)rf[5]/(float)cpf[1] );
fprintf(outfp,"R %5d %5.1f%%\n",
		rf[10], 100.*(float)rf[10]/(float)cpf[1] );
fprintf(outfp,"H=hydroxyl   %5d %5.1f%% :  S %5d %5.1f%%  T %5d %5.1f%%\n",
		cpf[2], 100.*(float)cpf[2]/(float)slgth,
		rf[3], 100.*(float)rf[3]/(float)cpf[2],
		rf[7], 100.*(float)rf[7]/(float)cpf[2] );
fprintf(outfp,"I=imino      %5d %5.1f%% :  P %5d\n",
		cpf[3], 100.*(float)cpf[3]/(float)slgth, rf[11] );
fprintf(outfp,"L=aliphatic  %5d %5.1f%% :  A %5d %5.1f%%  G %5d %5.1f%%  ",
		cpf[4], 100.*(float)cpf[4]/(float)slgth,
		rf[1], 100.*(float)rf[1]/(float)cpf[4],
		rf[2], 100.*(float)rf[2]/(float)cpf[4] );
fprintf(outfp,"I %5d %5.1f%%\n",
		rf[9], 100.*(float)rf[9]/(float)cpf[4] );
fprintf(outfp,"                          :  L %5d %5.1f%%  V %5d %5.1f%%\n",
		rf[0], 100.*(float)rf[0]/(float)cpf[4],
		rf[4], 100.*(float)rf[4]/(float)cpf[4] );
fprintf(outfp,"M=amide      %5d %5.1f%% :  N %5d %5.1f%%  Q %5d %5.1f%%\n",
		cpf[5], 100.*(float)cpf[5]/(float)slgth,
		rf[12], 100.*(float)rf[12]/(float)cpf[5],
		rf[14], 100.*(float)rf[14]/(float)cpf[5] );
fprintf(outfp,"R=aromatic   %5d %5.1f%% :  F %5d %5.1f%%  W %5d %5.1f%%  ",
		cpf[6], 100.*(float)cpf[6]/(float)slgth,
		rf[13], 100.*(float)rf[13]/(float)cpf[6],
		rf[19], 100.*(float)rf[19]/(float)cpf[6] );
fprintf(outfp,"Y %5d %5.1f%%\n",
		rf[15], 100.*(float)rf[15]/(float)cpf[6] );
fprintf(outfp,"S=sulfur     %5d %5.1f%% :  C %5d %5.1f%%  M %5d %5.1f%%\n\n",
		cpf[7], 100.*(float)cpf[7]/(float)slgth,
		rf[18], 100.*(float)rf[18]/(float)cpf[7],
		rf[17], 100.*(float)rf[17]/(float)cpf[7] );


fprintf(outfp,"\nFUNCTIONAL PROPERTY USAGE\n\n" );

fprintf(outfp,"H=nonpolar   %5d %5.1f%% :  A %5d %5.1f%%  F %5d %5.1f%%  ",
		fpf[0], 100.*(float)fpf[0]/(float)slgth,
		rf[1], 100.*(float)rf[1]/(float)fpf[0],
		rf[13], 100.*(float)rf[13]/(float)fpf[0] );
fprintf(outfp,"I %5d %5.1f%%\n",
		rf[9], 100.*(float)rf[9]/(float)fpf[0] );
fprintf(outfp,"                          :  L %5d %5.1f%%  M %5d %5.1f%%  ",
		rf[0], 100.*(float)rf[0]/(float)fpf[0],
		rf[17], 100.*(float)rf[17]/(float)fpf[0] );
fprintf(outfp,"P %5d %5.1f%%\n",
		rf[11], 100.*(float)rf[11]/(float)fpf[0] );
fprintf(outfp,"                          :  V %5d %5.1f%%  W %5d %5.1f%%\n",
		rf[4], 100.*(float)rf[4]/(float)fpf[0],
		rf[19], 100.*(float)rf[19]/(float)fpf[0] );
fprintf(outfp,"N=negative   %5d %5.1f%% :  D %5d %5.1f%%  E %5d %5.1f%%\n",
		fpf[1], 100.*(float)fpf[1]/(float)slgth,
		rf[8], 100.*(float)rf[8]/(float)fpf[1],
		rf[6], 100.*(float)rf[6]/(float)fpf[1] );
fprintf(outfp,"P=positive   %5d %5.1f%% :  H %5d %5.1f%%  K %5d %5.1f%%  ",
		fpf[2], 100.*(float)fpf[2]/(float)slgth,
		rf[16], 100.*(float)rf[16]/(float)fpf[2],
		rf[5], 100.*(float)rf[5]/(float)fpf[2] );
fprintf(outfp,"R %5d %5.1f%%\n",
		rf[10], 100.*(float)rf[10]/(float)fpf[2] );
fprintf(outfp,"U=polar      %5d %5.1f%% :  C %5d %5.1f%%  G %5d %5.1f%%  ",
		fpf[3], 100.*(float)fpf[3]/(float)slgth,
		rf[18], 100.*(float)rf[18]/(float)fpf[3],
		rf[2], 100.*(float)rf[2]/(float)fpf[3] );
fprintf(outfp,"N %5d %5.1f%%\n",
		rf[12], 100.*(float)rf[12]/(float)fpf[3] );
fprintf(outfp,"                          :  Q %5d %5.1f%%  S %5d %5.1f%%  ",
		rf[14], 100.*(float)rf[14]/(float)fpf[3],
		rf[3], 100.*(float)rf[3]/(float)fpf[3] );
fprintf(outfp,"T %5d %5.1f%%\n",
		rf[7], 100.*(float)rf[7]/(float)fpf[3] );
fprintf(outfp,"                          :  Y %5d %5.1f%%\n\n",
		rf[15], 100.*(float)rf[15]/(float)fpf[3] );


fprintf(outfp,"\nSTRUCTURAL PROPERTY USAGE\n\n" );

fprintf(outfp,"A=ambivalent %5d %5.1f%% :  A %5d %5.1f%%  C %5d %5.1f%%  ",
		spf[0], 100.*(float)spf[0]/(float)slgth,
		rf[1], 100.*(float)rf[1]/(float)spf[0],
		rf[18], 100.*(float)rf[18]/(float)spf[0] );
fprintf(outfp,"G %5d %5.1f%%\n",
		rf[2], 100.*(float)rf[2]/(float)spf[0] );
fprintf(outfp,"                          :  P %5d %5.1f%%  S %5d %5.1f%%  ",
		rf[11], 100.*(float)rf[11]/(float)spf[0],
		rf[3], 100.*(float)rf[3]/(float)spf[0] );
fprintf(outfp,"T %5d %5.1f%%\n",
		rf[7], 100.*(float)rf[7]/(float)spf[0] );
fprintf(outfp,"                          :  W %5d %5.1f%%  Y %5d %5.1f%%\n",
		rf[19], 100.*(float)rf[19]/(float)spf[0],
		rf[15], 100.*(float)rf[15]/(float)spf[0] );
fprintf(outfp,"E=external   %5d %5.1f%% :  D %5d %5.1f%%  E %5d %5.1f%%  ",
		spf[1], 100.*(float)spf[1]/(float)slgth,
		rf[8], 100.*(float)rf[8]/(float)spf[1],
		rf[6], 100.*(float)rf[6]/(float)spf[1] );
fprintf(outfp,"H %5d %5.1f%%\n",
		rf[16], 100.*(float)rf[16]/(float)spf[1] );
fprintf(outfp,"                          :  K %5d %5.1f%%  N %5d %5.1f%%  ",
		rf[5], 100.*(float)rf[5]/(float)spf[1],
		rf[12], 100.*(float)rf[12]/(float)spf[1] );
fprintf(outfp,"Q %5d %5.1f%%\n",
		rf[14], 100.*(float)rf[14]/(float)spf[1] );
fprintf(outfp,"                          :  R %5d %5.1f%%\n",
		rf[10], 100.*(float)rf[10]/(float)spf[1] );
fprintf(outfp,"I=internal   %5d %5.1f%% :  F %5d %5.1f%%  I %5d %5.1f%%  ",
		spf[2], 100.*(float)spf[2]/(float)slgth,
		rf[13], 100.*(float)rf[13]/(float)spf[2],
		rf[9], 100.*(float)rf[9]/(float)spf[2] );
fprintf(outfp,"L %5d %5.1f%%\n",
		rf[0], 100.*(float)rf[0]/(float)spf[2] );
fprintf(outfp,"                          :  M %5d %5.1f%%  V %5d %5.1f%%\n\n",
		rf[17], 100.*(float)rf[17]/(float)spf[2],
		rf[4], 100.*(float)rf[4]/(float)spf[2] );


fprintf(outfp,"\nAMINO ACID USAGE BY DEGREE OF CODON DEGENERACY\n\n" );

fprintf(outfp,"1 degeneracy %5d %5.1f%% :  M %5d %5.1f%%  W %5d %5.1f%%\n",
		dgf[0], 100.*(float)dgf[0]/(float)slgth,
		rf[17], 100.*(float)rf[17]/(float)dgf[0],
		rf[19], 100.*(float)rf[19]/(float)dgf[0] );
fprintf(outfp,"2 degeneracy %5d %5.1f%% :  C %5d %5.1f%%  D %5d %5.1f%%  ",
		dgf[1], 100.*(float)dgf[1]/(float)slgth,
		rf[18], 100.*(float)rf[18]/(float)dgf[1],
		rf[8], 100.*(float)rf[8]/(float)dgf[1] );
fprintf(outfp,"E %5d %5.1f%%\n",
		rf[6], 100.*(float)rf[6]/(float)dgf[1] );
fprintf(outfp,"                          :  F %5d %5.1f%%  H %5d %5.1f%%  ",
		rf[13], 100.*(float)rf[13]/(float)dgf[1],
		rf[16], 100.*(float)rf[16]/(float)dgf[1] );
fprintf(outfp,"K %5d %5.1f%%\n",
		rf[5], 100.*(float)rf[5]/(float)dgf[1] );
fprintf(outfp,"                          :  N %5d %5.1f%%  Q %5d %5.1f%%  ",
		rf[12], 100.*(float)rf[12]/(float)dgf[1],
		rf[14], 100.*(float)rf[14]/(float)dgf[1] );
fprintf(outfp,"Y %5d %5.1f%%\n",
		rf[15], 100.*(float)rf[15]/(float)dgf[1] );
fprintf(outfp,"3 degeneracy %5d %5.1f%% :  I %5d \n",
		dgf[2], 100.*(float)dgf[2]/(float)slgth, rf[9] );
fprintf(outfp,"4 degeneracy %5d %5.1f%% :  A %5d %5.1f%%  G %5d %5.1f%%  ",
		dgf[3], 100.*(float)dgf[3]/(float)slgth,
		rf[1], 100.*(float)rf[1]/(float)dgf[3],
		rf[2], 100.*(float)rf[2]/(float)dgf[3] );
fprintf(outfp,"P %5d %5.1f%%\n",
		rf[11], 100.*(float)rf[11]/(float)dgf[3] );
fprintf(outfp,"                          :  T %5d %5.1f%%  V %5d %5.1f%%\n",
		rf[7], 100.*(float)rf[7]/(float)dgf[3],
		rf[4], 100.*(float)rf[4]/(float)dgf[3] );
fprintf(outfp,"6 degeneracy %5d %5.1f%% :  L %5d %5.1f%%  R %5d %5.1f%%  ",
		dgf[4], 100.*(float)dgf[4]/(float)slgth,
		rf[0], 100.*(float)rf[0]/(float)dgf[4],
		rf[10], 100.*(float)rf[10]/(float)dgf[4] );
fprintf(outfp,"S %5d %5.1f%%\n",
		rf[3], 100.*(float)rf[3]/(float)dgf[4] );

} /* end print_propusage() */
