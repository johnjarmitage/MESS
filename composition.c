#include <stdio.h>
#include <math.h>
#include <omp.h>
#include <stdlib.h>
#include "main.h"
#include "composition.h"
#include "momentummarkerfuncs_par.h"

/* calculates the rare earth composition using parameterisation and equations from McKenzie and
   O'Nions (1991) J. of Petrology vol 32, 1021, and partition coefficeints Gibson and Geist, EPSL
   2010.

   i=0, La | i=1, Ce | i=2, Pr | i=3, Nd | i=4, Sm | i=5, Eu | i=6, Gd | i=7, Tb | i=8, Dy |
   i=9, Ho | i= 10, Er | i=11, Tm | i=12, Yb | i=13, Lu */

void mocompinst(){

  /* Gibson & Geist Supp. Matt., Earth and Planetary Science Letters, 2010 (+ McK & O'N for plag and spinel)*/
  double Dol[14] = {0.0005, 0.0005, 0.0008, 0.00042, 0.0011, 0.0016, 0.0011, 0.0015, 0.0027, 0.0016, 0.013, 0.0015, 0.02, 0.02};
  double Dopx[14] = {0.0031, 0.004, 0.0048, 0.012, 0.02, 0.013, 0.0065, 0.019, 0.011, 0.026, 0.045, 0.04, 0.08, 0.12};
  double Dcpx[14] = {0.049, 0.08, 0.126, 0.178, 0.293, 0.335, 0.35, 0.403, 0.400, 0.427, 0.420, 0.40968, 0.400, 0.376};
  double Dplag[14] = {0.27, 0.20, 0.17, 0.14, 0.11, 0.73, 0.066, 0.06, 0.055, 0.048, 0.041, 0.036, 0.031, 0.025};
  double Dspinel[14] = {0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01};
  double Dgarnet[14] = {0.001, 0.005, 0.014, 0.052, 0.250, 0.496, 0.848, 1.477, 2.200, 3.315, 4.400, 5.495, 6.600, 7.100};
  double Damph[14] = {0.17, 0.26, 0.35, 0.44, 0.76, 0.88, 0.86, 0.83, 0.78, 0.73, 0.68, 0.64, 0.59, 0.51};

  /* La, Sm and Yb replaced with McDade et al. PEPI, 2003 -- MPY-90 data is interpolated */
  //double Dol[14] =     {8.3e-6, 0.0005, 0.0008, 0.00042, 8.6e-4, 0.0016, 0.0011, 0.0015, 0.0027, 0.0016, 0.013, 0.0015,  0.036, 0.02};
  //double Dopx[14] =    {0.003,  0.004,  0.0048, 0.012,   0.068,  0.013,  0.0065, 0.019,  0.011,  0.026,  0.045, 0.04,    0.25,  0.12};
  //double Dcpx[14] =    {0.089,  0.08,   0.126,  0.178,   0.67,   0.335,  0.35,   0.403,  0.400,  0.427,  0.420, 0.40968, 1.43,  0.376};
  //double Dplag[14] =   {0.27,   0.20,   0.17,   0.14,    0.11,   0.73,   0.066,  0.06,   0.055,  0.048,  0.041, 0.036,   0.031, 0.025};
  //double Dspinel[14] = {0.01,   0.01,   0.01,   0.01,    0.01,   0.01,   0.01,   0.01,   0.01,   0.01,   0.01,  0.01,    0.01,  0.01};
  //double Dgarnet[14] = {0.001,  0.005,  0.014,  0.052,   0.210,  0.496,  0.848,  1.477,  2.200,  3.315,  4.400, 5.495,   4.18,  7.100};
  //double Damph[14] =   {0.17,   0.26,   0.35,   0.44,    0.76,   0.88,   0.86,   0.83,   0.78,   0.73,   0.68,  0.64,    0.59,  0.51};
  
  /*                  Ol     Opx    Cpx       */
  double Plag[4] = {0.636, 0.263, 0.012, 0.089};
  double Spinel[4] = {0.578, 0.270, 0.119, 0.033};
  double Garnet[4] = {0.598, 0.211, 0.079, 0.115};
  double Amph[4] = {0.599, 0.247, 0.038, 0.116};
  int comp;
  float noz;
  double p,phi,D,maxD,Da,Db,m,c,Pspinel_out,Pgarnet_in;
  
  double F,CompF,*gridX,*gridXComp,*weightX,*weightXComp,integratedX,integratedXComp,maxX,maxXComp;
  maxD = 0;

  int e,t; /* to stop warnings - but does it matter? */
#pragma omp parallel private(e,p,Pspinel_out,Pgarnet_in,comp,D,Da,Db,m,c)
  {

#pragma omp for schedule(static) reduction(max : maxD)
    for (e=0;e<N;e++) {     

      p = 1e-9*mP[e]; /* in Gpa */
      F = 1 - 1/mX[e];
      
      /* and what about spinel-out and garnet-in, in GPa */
      Pspinel_out = (mTemp[e]+400)/666.7;
      Pgarnet_in = (mTemp[e]+533)/666.7;
      
      /* if melt is present what is the partition coefficient and composition */
      if (mdm[e] > 0) {
	for (comp=0;comp<num_comp;comp++) {

	  /* shallow first:
	     if < 25 km (0.81 GPa) deep we are in the plagioclase stability
	     field
	     Which I have converted to pressure assuming a density of 3300
	     kgm3, as we have an upper air layer in the model! */
	  if (p <= 0.81)
	    D = Dol[comp]*Plag[0] + Dopx[comp]*Plag[1] + Dcpx[comp]*Plag[2] + Dplag[comp]*Plag[3];	  
	  /* if between stability, linear transition... */
	  if (p > 0.81 && p <= 1.13) {
	    Da = Dol[comp]*Plag[0] + Dopx[comp]*Plag[1] + Dcpx[comp]*Plag[2] + Dplag[comp]*Plag[3];
	    Db = Dol[comp]*Spinel[0] + Dopx[comp]*Spinel[1] + Dcpx[comp]*Spinel[2] + Dspinel[comp]*Spinel[3];
	    m = (Db-Da)/0.32;
	    c = Da-m*1.13;
	    D = m*p+c;
	  }
	  
	  /* if >35 km (1.13 GPa) and < spinel-out
	     Which I have converted to pressure assuming a density of 3300
	     kgm3, as we have an upper air layer in the model! */
	  if (p > 1.13 && p <= Pgarnet_in)
	    D = Dol[comp]*Spinel[0] + Dopx[comp]*Spinel[1] + Dcpx[comp]*Spinel[2] + Dspinel[comp]*Spinel[3];
	  
	  /* again between stability, linear transition... */
	  if (p > Pspinel_out && p <= Pgarnet_in) {
	    Da = Dol[comp]*Spinel[0] + Dopx[comp]*Spinel[1] + Dcpx[comp]*Spinel[2] + Dspinel[comp]*Spinel[3];
	    Db = Dol[comp]*Garnet[0] + Dopx[comp]*Garnet[1] + Dcpx[comp]*Garnet[2] + Dgarnet[comp]*Garnet[3];
	    m = (Db-Da)/(Pgarnet_in-Pspinel_out);
	    c = Da-m*Pspinel_out;
	    D = m*p+c;
	  }
	  
	  /* if > garnet_in */
	  if (p > Pgarnet_in)
	    D = Dol[comp]*Garnet[0] + Dopx[comp]*Garnet[1] + Dcpx[comp]*Garnet[2] + Dgarnet[comp]*Garnet[3];

	  if (comp == 0 && D > maxD) maxD = D;
	  
	  if (D > 0) { /* just to check that we have sensible values for D */
	    
	    /* Cl = Cs/D */
	    mXComp[comp][e] = mRComp[comp][e]/D;
	    
	    /* Cs_new*(1-F) + F*Cl = Cs_old -- mass balance batch melting? */
	    mRComp[comp][e] = mRComp[comp][e]/(1 + mdm[e]*(1/D - 1));
	    
	    /* Cs cannot be less than zero! */
	    if (mRComp[comp][e] <= 0) mRComp[comp][e] = 0;
	    
	  }
	}
      }
      if (mdm[e] <= 0) {
	for (comp=0;comp<num_comp;comp++)
	  mXComp[comp][e] = 0;
      }
    }
  }

  // Here marker values of dXComp*mX and mX are interpolated to a grid (without dividing correcting for weights yet).
  // Each marker point updates the value of the 4 grid points that corner the cell of the marker (the indices for the cell are mcs[e],mct[e])

  if (maxD > 0) {

    gridX=cx0; weightX=cz0; resetbasegrid(gridX); resetbasegrid(weightX);
#pragma omp parallel private(e)
    {
      int Ncells, Nt, Nblocks, threadmustworkatthisblock, localblock, odd;
      Ncells=Sbase-1; Nt=omp_get_num_threads(); Nblocks=2*Nt;
      for (odd=0;odd<2;odd++) {
#pragma omp barrier
	threadmustworkatthisblock=2*omp_get_thread_num()+odd;
	for (e=0;e<N;e++)    {
	  F = 0;	  
	  if (mdm[e] > 0) F = mdm[e]; /* small degree melting is not pooled ? */
	  localblock=(mcs[e]*Nblocks)/Ncells;
	  if (localblock==threadmustworkatthisblock) {
	    markertorhogrid(e, F, gridX, weightX);
	  }
	}
      }
    }
    // Here weights are divided by, finishing the interpolation to the grid
#pragma omp parallel for schedule(static) private(s, t, st)
    for (t=0;t<Tbase;t++) {
      for (s=0;s<Sbase;s++) {
	st=IN0(s, t);
	if (weightX[st]>0) gridX[st]/=weightX[st]; else { gridX[st]=0; }
      }
    }
    // And here the grids are integraded
    integratedX=0;
#pragma omp parallel for schedule(static) private(s, t, st) reduction(+:integratedX)
    for (t=0;t<Tbase-1;t++) {
      for (s=0;s<Sbase-1;s++) {
	integratedX+=0.25*(gridX[IN0(s, t)]+gridX[IN0(s+1, t)]+gridX[IN0(s, t+1)]+gridX[IN0(s+1, t+1)])*dx[s]*dz[t];
      }
    }

    for (comp=0;comp<num_comp;comp++) {
  
      gridXComp=cx0; weightXComp=cz0; resetbasegrid(gridXComp); resetbasegrid(weightXComp);
#pragma omp parallel private(e)
      {
	int Ncells, Nt, Nblocks, threadmustworkatthisblock, localblock, odd;
	Ncells=Sbase-1; Nt=omp_get_num_threads(); Nblocks=2*Nt;
	for (odd=0;odd<2;odd++) {
#pragma omp barrier
	  threadmustworkatthisblock=2*omp_get_thread_num()+odd;
	  for (e=0;e<N;e++)    {
	    CompF = 0;	  
	    if (mdm[e] > 0)  CompF = mXComp[comp][e]*mdm[e]; /* small degree melting is not pooled ? */
	    localblock=(mcs[e]*Nblocks)/Ncells;
	    if (localblock==threadmustworkatthisblock) {
	      markertorhogrid(e, CompF, gridXComp, weightXComp);
	    }
	  }
	}
      }
      // Here weights are divided by, finishing the interpolation to the grid
#pragma omp parallel for schedule(static) private(s, t, st)
      for (t=0;t<Tbase;t++) {
	for (s=0;s<Sbase;s++) {
	  st=IN0(s, t);
	  if (weightXComp[st]>0) gridXComp[st]/=weightXComp[st]; else { gridXComp[st]=0; }
	}
      }
      // And here the grids are integraded
      integratedXComp=0;
#pragma omp parallel for schedule(static) private(s, t, st) reduction(+:integratedXComp)
      for (t=0;t<Tbase-1;t++) {
	for (s=0;s<Sbase-1;s++) {
	  integratedXComp+=0.25*(gridXComp[IN0(s, t)]+gridXComp[IN0(s+1, t)]+gridXComp[IN0(s, t+1)]+gridXComp[IN0(s+1, t+1)])*dx[s]*dz[t];
	}
      }
      
      bulkcomp[comp] = integratedXComp/integratedX;
      if (comp == 0)
	fprintf(statusfile,"\nMelt Partition: max DLa %g, CX*F %g, meanF %g.\n",maxD,integratedXComp,integratedX);
    }
  }
  else {
    fprintf(statusfile,"\nMelt Partition should be zero: max DLa %g.\n",maxD);
    for (comp=0;comp<num_comp;comp++)
      bulkcomp[comp] = 0;
   }
}

void outputmeltcomp() {
  FILE *file0;
  char *writemode;
  if (iterationnumber==0) writemode="wb"; else writemode="ab";
  
  file0 = fopen("outfiles/bulk_comp.txt", writemode);
  
  fwrite(bulkcomp, sizeof(double), num_comp, file0);
    
  fclose(file0);
}
