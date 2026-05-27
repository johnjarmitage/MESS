#ifndef MAIN_H
#define MAIN_H


//should some integers be unsigned
extern int N,Sbase,Tbase;
extern int s,t,st,e,S,T,G,multigridlevels;
extern int *mcs,*mct,*mtype;
extern int leftbound,rightbound,upbound,downbound; //0: no slip, 1: free slip 2 free slip/free external slip
extern double Tp,Ttop,Tbottom;
extern double leftboundvel,rightboundvel,upboundvel,downboundvel;
extern double vold,voldold,Iold,Ioldold,stressintegral;
extern int *firstelementgrid,*firstelementdx,*firstelementdz,*pow2;
extern int numberofyieldmarkers,firstmarkerupdate;
extern int iterationnumber;
extern double time,dt,dtnom,dtmaxwell;
extern double *dx,*dz,*x,*z,mincellsize,xmin,xmax,zmin,zmax;
extern double *mx,*mz,*mTemp,*mstressxx,*mstressxz,*mstrainxx,*mstrainxz;
extern double *nZn,*nZs,*inZn,*inZs; //ekstra arrays for varying viscosity contrasts
extern double *weight,*rightP,*rightx,*rightz,*resP,*resx,*resz,*cx0,*cx1,*cx2,*cx3,*cx4,*cx5,*cx6,*cx7,*cx8,*cx9,*cx10,*cz0,*cz1,*cz2,*cz3,*cz4,*cz5,*cz6,*cz7,*cz8,*cz9,*cz10;
extern double *vx0,*vz0,*P0,*rightx0,*rightz0,*rightP0; //only needed for multigrid
extern double *P,*vx,*vz;
extern double *gridviscn,*gridviscs,*gridrho,*gridoldstressxx,*gridoldstressxz,*gridmyn,*gridmys; //some of these grids  do not get allocated independent memory, but is use memory allocated to the cx coeffs which is not used simultaneously
extern double *gridHs; //global variable since used in both stressstrainupdate and updatetemp
extern double *mrho,*mvisc,*mmy,*mP,*mfinitestrain;
extern double gx,gz;
extern double minvisc,maxvisc,alphaP,alphamom,alphacorrect,stressprecision;
extern FILE *statusfile,*surfacefile;
extern int Nsurfaces,sur;
extern double *surfaces;
extern double *sedipool;
extern double sealevel;
//update: Lookup table of S in IN. Nthreads, max number of threads, temperature arrays global
extern int *Stable;
extern int Nthreads;
extern double *mC, *mk, *mHr;

//variable for age tracking
extern double *mage;
extern int convertcenterpoints;

// melting
extern int num_comp,comp;
extern double *mX,*mdm,*mWater,meltproduction,meltproductionVOL,maxmeltpos;
extern double **mXComp,**mRComp,*bulkcomp;

// ice
extern int glacier,deglaciate,markstart,markend;
extern double icethick,deglaciation_start,deglaciation_end;


// Functions

int IN(int s,int t,int G);
int IN0(int s,int t);
int INdx(int s, int G);
int INdz(int t, int G);
int INsur(int s, int sur);
int INsur(int s, int sur);

int calculatemarkerparams(int yield);

#endif
