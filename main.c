#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <omp.h>
#include "main.h"
#include "momentummarkerfuncs_par.h"
#include "multigridfuncs2.h"
#include "tempfuncs_par.h"
#include "surfacefuncs.h"
#include "composition.h"


//should some integers be unsigned
int N,Sbase,Tbase;
int s,t,st,e,S,T,G,multigridlevels;
int *mcs,*mct,*mtype;
int leftbound,rightbound,upbound,downbound; //0: no slip, 1: free slip 2 free slip/free external slip
double Tp,Ttop,Tbottom;
double leftboundvel,rightboundvel,upboundvel,downboundvel;
double vold,voldold,Iold,Ioldold,stressintegral;
int *firstelementgrid,*firstelementdx,*firstelementdz,*pow2;
int numberofyieldmarkers,firstmarkerupdate;
int iterationnumber;
double time,dt,dtnom,dtmaxwell;
double *dx,*dz,*x,*z,mincellsize,xmin,xmax,zmin,zmax;
double *mx,*mz,*mTemp,*mstressxx,*mstressxz,*mstrainxx,*mstrainxz;
double *nZn,*nZs,*inZn,*inZs; //ekstra arrays for varying viscosity contrasts
double *weight,*rightP,*rightx,*rightz,*resP,*resx,*resz,*cx0,*cx1,*cx2,*cx3,*cx4,*cx5,*cx6,*cx7,*cx8,*cx9,*cx10,*cz0,*cz1,*cz2,*cz3,*cz4,*cz5,*cz6,*cz7,*cz8,*cz9,*cz10;
double *vx0,*vz0,*P0,*rightx0,*rightz0,*rightP0; //only needed for multigrid
double *P,*vx,*vz;
double *gridviscn,*gridviscs,*gridrho,*gridoldstressxx,*gridoldstressxz,*gridmyn,*gridmys; //some of these grids  do not get allocated independent memory, but is use memory allocated to the cx coeffs which is not used simultaneously
double *gridHs; //global variable since used in both stressstrainupdate and updatetemp
double *mrho,*mvisc,*mmy,*mP,*mfinitestrain;
double gx,gz;
double minvisc,maxvisc,alphaP,alphamom,alphacorrect,stressprecision;
FILE *statusfile,*surfacefile;
int Nsurfaces,sur;
double *surfaces;
double *sedipool;
double sealevel;
//update: Lookup table of S in IN. Nthreads, max number of threads, temperature arrays global
int *Stable;
int Nthreads;
double *mC, *mk, *mHr;

//variable for age tracking
double *mage;
int convertcenterpoints;

// melting
int num_comp,comp;
double *mX,*mdm,*mWater,meltproduction,meltproductionVOL,maxmeltpos;
double **mXComp,**mRComp,*bulkcomp;

// ice
int glacier,deglaciate,markstart,markend;
double icethick,deglaciation_start,deglaciation_end;

double riftvelo = 0;

/*
#include "momentmummarkerfuncs_par.c"
#include "multigridfuncs2.c"
#include "tempfuncs_par.c"
#include "surfacefuncs.c"
#include "composition.c"
//use other memory fro gridoldstress?
*/

int IN(int s,int t,int G) {
    return t*Stable[G]+s+firstelementgrid[G];
    //return t*((Sbase-1)/pow2[G]+1)+s+firstelementgrid[G];
}
int IN0(int s,int t) {
    return t*Sbase+s;
}
int INdx(int s, int G) {
    return s+firstelementdx[G];
}
int INdz(int t, int G) {
    return t+firstelementdz[G];
}
int INsur(int s, int sur) {
    return s+sur*Sbase;
}
/*
//input: 0: boundary has 0 perpendicular flow. >0: flow out of the model domain. <0: flow into the model domain
void velocityboundarychange(int left,int right,int up,int down,double inwardsvelo) {
    double totalinwardsflux,totaloutwardslength; //absolute values of fluxes are so far considered
    totalinwardsflux=0; totaloutwardslength=0;
    leftboundvel=0; rightboundvel=0; upboundvel=0; downboundvel=0;
    if (left<0) { leftboundvel=inwardsvelo; totalinwardsflux+=inwardsvelo*(zmax-zmin); } else if (left>0) totaloutwardslength+=zmax-zmin;
    if (right<0) { rightboundvel=-inwardsvelo; totalinwardsflux+=inwardsvelo*(zmax-zmin); } else if (right>0) totaloutwardslength+=zmax-zmin;
    if (up<0) { upboundvel=-inwardsvelo; totalinwardsflux+=inwardsvelo*(xmax-xmin); } else if (up>0) totaloutwardslength+=xmax-xmin;
    if (down<0) { downboundvel=inwardsvelo; totalinwardsflux+=inwardsvelo*(xmax-xmin); } else if (down>0) totaloutwardslength+=xmax-xmin;
    //constant outvel: totaloutwardslength*outvel=totalinwardsflux
    if (left>0) leftboundvel=-totalinwardsflux/totaloutwardslength;
    if (right>0) rightboundvel=totalinwardsflux/totaloutwardslength;
    if (up>0) upboundvel=-totalinwardsflux/totaloutwardslength;
    if (down>0) downboundvel=totalinwardsflux/totaloutwardslength;
    mexPrintf("leftboundvel,rightboundvel,upboundvel,downboundvel, %g %g %g %g\n",leftboundvel,rightboundvel,upboundvel,downboundvel);
    mexPrintf("inwards %g net flow %g\n",totalinwardsflux,rightboundvel*(zmax-zmin)-leftboundvel*(zmax-zmin)+downboundvel*(xmax-xmin)-upboundvel*(xmax-xmin));
}
*/
//input: 0: boundary has 0 perpendicular flow. >0: flow out of the model domain. <0: flow into the model domain
void velocityboundarychange(int left,int right,int up,int down,double inwardsvelo) {
    double totalinwardsflux,totaloutwardslength; //absolute values of fluxes are so far considered
    totalinwardsflux=0; totaloutwardslength=0;
    leftboundvel=0; rightboundvel=0; upboundvel=0; downboundvel=0;
    if (left<0) { leftboundvel=inwardsvelo; totalinwardsflux+=inwardsvelo*(zmax-zmin); } else if (left>0) totaloutwardslength+=zmax-zmin;
    if (right<0) { rightboundvel=-inwardsvelo; totalinwardsflux+=inwardsvelo*(zmax-zmin); } else if (right>0) totaloutwardslength+=zmax-zmin;
    if (up<0) { upboundvel=inwardsvelo; totalinwardsflux+=inwardsvelo*(xmax-xmin); } else if (up>0) totaloutwardslength+=xmax-xmin;
    if (down<0) { downboundvel=-inwardsvelo; totalinwardsflux+=inwardsvelo*(xmax-xmin); } else if (down>0) totaloutwardslength+=xmax-xmin;
    //constant outvel: totaloutwardslength*outvel=totalinwardsflux
/*    if (left>0) leftboundvel=-(zmax-zmin)/totaloutwardslength*totalinwardsflux/totaloutwardslength;
    if (right>0) rightboundvel=(zmax-zmin)/totaloutwardslength*totalinwardsflux/totaloutwardslength;
    if (up>0) upboundvel=-(xmax-xmin)/totaloutwardslength*totalinwardsflux/totaloutwardslength;
    if (down>0) downboundvel=(xmax-xmin)/totaloutwardslength*totalinwardsflux/totaloutwardslength;*/
    if (left>0) leftboundvel=-totalinwardsflux/totaloutwardslength;
    if (right>0) rightboundvel=totalinwardsflux/totaloutwardslength;
    if (up>0) upboundvel=-totalinwardsflux/totaloutwardslength;
    if (down>0) downboundvel=totalinwardsflux/totaloutwardslength;
    fprintf(statusfile,"leftboundvel,rightboundvel,upboundvel,downboundvel, %g %g %g %g\n",leftboundvel,rightboundvel,upboundvel,downboundvel);
    fprintf(statusfile,"inwards %g net flow %g\n",totalinwardsflux,rightboundvel*(zmax-zmin)-leftboundvel*(zmax-zmin)+downboundvel*(xmax-xmin)-upboundvel*(xmax-xmin));
}

void velocityboundarychange1(int left,int right,int up,int down,double outwardsvelo) {
    double totaloutwardsflux,totalinwardslength; //absolute values of fluxes are so far considered
    totaloutwardsflux=0; totalinwardslength=0;
    leftboundvel=0; rightboundvel=0; upboundvel=0; downboundvel=0;
    if (left>0) { leftboundvel=-outwardsvelo; totaloutwardsflux+=outwardsvelo*(zmax-zmin); } else if (left<0) totalinwardslength+=zmax-zmin;
    if (right>0) { rightboundvel=outwardsvelo; totaloutwardsflux+=outwardsvelo*(zmax-zmin); } else if (right<0) totalinwardslength+=zmax-zmin;
    if (up>0) { upboundvel=-outwardsvelo; totaloutwardsflux+=outwardsvelo*(xmax-xmin); } else if (up<0) totalinwardslength+=xmax-xmin;
    if (down>0) { downboundvel=outwardsvelo; totaloutwardsflux+=outwardsvelo*(xmax-xmin); } else if (down<0) totalinwardslength+=xmax-xmin;

    if (left<0) leftboundvel=totaloutwardsflux/totalinwardslength;
    if (right<0) rightboundvel=-totaloutwardsflux/totalinwardslength;
    if (up<0) upboundvel=totaloutwardsflux/totalinwardslength;
    if (down<0) downboundvel=-totaloutwardsflux/totalinwardslength;
    fprintf(statusfile,"leftboundvel,rightboundvel,upboundvel,downboundvel, %g %g %g %g\n",leftboundvel,rightboundvel,upboundvel,downboundvel);
    fprintf(statusfile,"outwards %g net flow %g\n",totaloutwardsflux,rightboundvel*(zmax-zmin)-leftboundvel*(zmax-zmin)+downboundvel*(xmax-xmin)-upboundvel*(xmax-xmin));
}

void initialize() {
    FILE *file;
    file=fopen("0.txt", "rb");
    //Scalar data from beginning of file
    fread(&N,sizeof(int),1,file);
    fread(&Sbase,sizeof(int),1,file);
    fread(&Tbase,sizeof(int),1,file);
    fread(&multigridlevels,sizeof(int),1,file);
    fread(&leftbound,sizeof(int),1,file);
    fread(&rightbound,sizeof(int),1,file);
    fread(&upbound,sizeof(int),1,file);
    fread(&downbound,sizeof(int),1,file);
    fread(&time,sizeof(double),1,file);
    fread(&dtnom,sizeof(double),1,file); dtmaxwell=dtnom;
    //allocation and calculation of basic indexing values
    firstelementgrid=calloc(multigridlevels+2, sizeof(int)); //these show the index of the first element in a grid or dx/dz-vector. Last value [multigridlevels+1] is the total numbers of values
    firstelementdx=calloc(multigridlevels+2, sizeof(int));
    firstelementdz=calloc(multigridlevels+2, sizeof(int));
    pow2=calloc(multigridlevels+1, sizeof(int));
    Stable=calloc(multigridlevels+1, sizeof(int));
    S=Sbase;
    T=Tbase;
    G=0;
    firstelementgrid[0]=0;
    firstelementdx[0]=0;
    firstelementdz[0]=0;
    firstelementgrid[1]=Sbase*Tbase;
    firstelementdx[1]=Sbase-1;
    firstelementdz[1]=Tbase-1;
    for (G=2;G<=multigridlevels+1;G++){
        if ((S-1)%2!=0 || (T-1)%2!=0 || S==1 || T==1) {
            fprintf(statusfile,"Error: multigridlevel too high or bad array size\n"); }
        S=(S-1)/2+1;
        T=(T-1)/2+1;
        firstelementgrid[G]=firstelementgrid[G-1]+S*T;
        firstelementdx[G]=firstelementdx[G-1]+S-1;
        firstelementdz[G]=firstelementdz[G-1]+T-1;
    }
    for (G=0;G<=multigridlevels+1;G++){
        fprintf(statusfile,"%d %d %d\n",firstelementgrid[G],firstelementdx[G],firstelementdz[G]); }
    pow2[0]=1; //Is used in the function IN() since pow() takes many ressources
    for (G=1;G<=multigridlevels;G++) pow2[G]=pow2[G-1]*2;
    for (G=0;G<=multigridlevels;G++) Stable[G]=(Sbase-1)/pow2[G]+1;
    //allocation of grid memory
    dx=calloc(firstelementdx[multigridlevels+1],sizeof(double));
    dz=calloc(firstelementdz[multigridlevels+1],sizeof(double));
    P=calloc(firstelementgrid[multigridlevels+1],sizeof(double));
    vx=calloc(firstelementgrid[multigridlevels+1],sizeof(double));
    vz=calloc(firstelementgrid[multigridlevels+1],sizeof(double));
    nZs=calloc(firstelementgrid[multigridlevels+1],sizeof(double));
    nZn=calloc(firstelementgrid[multigridlevels+1],sizeof(double));
    rightP=calloc(firstelementgrid[multigridlevels+1],sizeof(double));
    rightx=calloc(firstelementgrid[multigridlevels+1],sizeof(double));
    rightz=calloc(firstelementgrid[multigridlevels+1],sizeof(double));
    resP=calloc(firstelementgrid[multigridlevels+1],sizeof(double));
    resx=calloc(firstelementgrid[multigridlevels+1],sizeof(double));
    resz=calloc(firstelementgrid[multigridlevels+1],sizeof(double));
    cx0=calloc(firstelementgrid[multigridlevels+1],sizeof(double));
    cx1=calloc(firstelementgrid[multigridlevels+1],sizeof(double));
    cx2=calloc(firstelementgrid[multigridlevels+1],sizeof(double));
    cx3=calloc(firstelementgrid[multigridlevels+1],sizeof(double));
    cx4=calloc(firstelementgrid[multigridlevels+1],sizeof(double));
    cx5=calloc(firstelementgrid[multigridlevels+1],sizeof(double));
    cx6=calloc(firstelementgrid[multigridlevels+1],sizeof(double));
    cx7=calloc(firstelementgrid[multigridlevels+1],sizeof(double));
    cx8=calloc(firstelementgrid[multigridlevels+1],sizeof(double));
    cx9=calloc(firstelementgrid[multigridlevels+1],sizeof(double));
    cx10=calloc(firstelementgrid[multigridlevels+1],sizeof(double));
    cz0=calloc(firstelementgrid[multigridlevels+1],sizeof(double));
    cz1=calloc(firstelementgrid[multigridlevels+1],sizeof(double));
    cz2=calloc(firstelementgrid[multigridlevels+1],sizeof(double));
    cz3=calloc(firstelementgrid[multigridlevels+1],sizeof(double));
    cz4=calloc(firstelementgrid[multigridlevels+1],sizeof(double));
    cz5=calloc(firstelementgrid[multigridlevels+1],sizeof(double));
    cz6=calloc(firstelementgrid[multigridlevels+1],sizeof(double));
    cz7=calloc(firstelementgrid[multigridlevels+1],sizeof(double));
    cz8=calloc(firstelementgrid[multigridlevels+1],sizeof(double));
    cz9=calloc(firstelementgrid[multigridlevels+1],sizeof(double));
    cz10=calloc(firstelementgrid[multigridlevels+1],sizeof(double));
    P0=calloc(firstelementgrid[1],sizeof(double));
    vx0=calloc(firstelementgrid[1],sizeof(double));
    vz0=calloc(firstelementgrid[1],sizeof(double));
    rightx0=calloc(firstelementgrid[1],sizeof(double));
    rightz0=calloc(firstelementgrid[1],sizeof(double));
    rightP0=calloc(firstelementgrid[1],sizeof(double));
    weight=calloc(firstelementgrid[multigridlevels+1],sizeof(double));
    P0=calloc(firstelementgrid[1],sizeof(double));
    vx0=calloc(firstelementgrid[1],sizeof(double));
    vz0=calloc(firstelementgrid[1],sizeof(double));
    rightx0=calloc(firstelementgrid[1],sizeof(double));
    rightz0=calloc(firstelementgrid[1],sizeof(double));
    rightP0=calloc(firstelementgrid[1],sizeof(double));
    gridrho=calloc(firstelementgrid[1],sizeof(double));
    gridoldstressxx=calloc(firstelementgrid[1],sizeof(double));
    gridoldstressxz=calloc(firstelementgrid[1],sizeof(double));
    gridviscn=calloc(firstelementgrid[1],sizeof(double));
    gridviscs=calloc(firstelementgrid[1],sizeof(double));
    gridmyn=calloc(firstelementgrid[1],sizeof(double));
    gridmys=calloc(firstelementgrid[1],sizeof(double));
    inZn=calloc(firstelementgrid[1],sizeof(double));
    inZs=calloc(firstelementgrid[1],sizeof(double));
    //allocation and input of grid and marker data
/*    for (s=0;s<Sbase-1;s++) fscanf(file,"%lg",&dx[s]);
    for (t=0;t<Tbase-1;t++) fscanf(file,"%lg",&dz[t]);
    mx=calloc(N,sizeof(double)); for (e=0;e<N;e++) fscanf(file,"%lg",&mx[e]);
    mz=calloc(N,sizeof(double)); for (e=0;e<N;e++) fscanf(file,"%lg",&mz[e]);
    mcs=calloc(N,sizeof(int)); for (e=0;e<N;e++) fscanf(file,"%d",&mcs[e]);
    mct=calloc(N,sizeof(int)); for (e=0;e<N;e++) fscanf(file,"%d",&mct[e]);
    mtype=calloc(N,sizeof(int)); for (e=0;e<N;e++) fscanf(file,"%d",&mtype[e]);
    mTemp=calloc(N,sizeof(double)); for (e=0;e<N;e++) fscanf(file,"%lg",&mTemp[e]);
    mstressxx=calloc(N,sizeof(double)); for (e=0;e<N;e++) fscanf(file,"%lg",&mstressxx[e]);
    mstressxz=calloc(N,sizeof(double)); for (e=0;e<N;e++) fscanf(file,"%lg",&mstressxz[e]);
    mstrainxx=calloc(N,sizeof(double)); for (e=0;e<N;e++) fscanf(file,"%lg",&mstrainxx[e]);
    mstrainxz=calloc(N,sizeof(double)); for (e=0;e<N;e++) fscanf(file,"%lg",&mstrainxz[e]);
    mfinitestrain=calloc(N,sizeof(double)); for (e=0;e<N;e++) fscanf(file,"%lg",&mfinitestrain[e]);*/
    //allocation and input of grid and marker data
    fread(dx,sizeof(double),Sbase-1,file);
    fread(dz,sizeof(double),Tbase-1,file);
    mx=calloc(N,sizeof(double)); fread(mx,sizeof(double),N,file);
    mz=calloc(N,sizeof(double)); fread(mz,sizeof(double),N,file);
    mcs=calloc(N,sizeof(int)); fread(mcs,sizeof(int),N,file);
    mct=calloc(N,sizeof(int)); fread(mct,sizeof(int),N,file);
    mtype=calloc(N,sizeof(int)); fread(mtype,sizeof(int),N,file);
    mTemp=calloc(N,sizeof(double)); fread(mTemp,sizeof(double),N,file);
    mstressxx=calloc(N,sizeof(double)); fread(mstressxx,sizeof(double),N,file);
    mstressxz=calloc(N,sizeof(double)); fread(mstressxz,sizeof(double),N,file);
    mstrainxx=calloc(N,sizeof(double)); fread(mstrainxx,sizeof(double),N,file);
    mstrainxz=calloc(N,sizeof(double)); fread(mstrainxz,sizeof(double),N,file);
    mfinitestrain=calloc(N,sizeof(double)); fread(mfinitestrain,sizeof(double),N,file);
    //allocation of additional marker memory
    mvisc=calloc(N,sizeof(double));
    mmy=calloc(N,sizeof(double));
    mrho=calloc(N,sizeof(double));
    mP=calloc(N,sizeof(double));
    mC=calloc(N, sizeof(double));
    mk=calloc(N, sizeof(double));
    mHr=calloc(N, sizeof(double));
    //input of grid data
    fread(vx,sizeof(double),Sbase*Tbase,file);
    fread(vz,sizeof(double),Sbase*Tbase,file);
    fread(P,sizeof(double),Sbase*Tbase,file);
    
    //input of surfaces
    fread(&Nsurfaces,sizeof(int),1,file);
    surfaces=calloc(Nsurfaces*Sbase,sizeof(double));
    for (sur=0;sur<Nsurfaces;sur++) fread(surfaces+sur*Sbase,sizeof(double),Sbase,file);
    fprintf(statusfile,"Nsurfaces %d\n",Nsurfaces);
    fclose(file);
    
    //age tracking initiation
    mage=calloc(N, sizeof(double)); for (e=0;e<N;e++) mage[e]=-1.0;
    
    //depletion tracking
    mX=calloc(N, sizeof(double));
    mdm=calloc(N, sizeof(double));
    for (e=0;e<N;e++) {
        mX[e]=1;
	mdm[e]=0;
	// if (mtype[e]==3) mX[e]=1.2;
	// if (mz[e]>=40e3 && mz[e]<175e3+40e3) mX[e]=1.0/(0.5+2.8571e-6*(mz[e]-40e3)); else mX[e]=1;
    }

    // water content tracking
    mWater=calloc(N, sizeof(double));
    for (e=0;e<N;e++) {
        mWater[e]=100;
    }

    //melt compostion tracking
    num_comp = 14;
    mXComp=(double **)calloc(num_comp, sizeof(double *));
    mRComp=(double **)calloc(num_comp, sizeof(double *));
    bulkcomp=calloc(num_comp, sizeof(double));
    for (comp=0;comp<num_comp;comp++){
      mXComp[comp] = (double *)calloc(N, sizeof(double));
      mRComp[comp] = (double *)calloc(N, sizeof(double));
    }
    for (e=0;e<N;e++) {
      for (comp=0;comp<num_comp;comp++)
	mXComp[comp][e]=0;
       /* i=0, La | i=1,   Ce | i=2,  Pr | i=3,  Nd | i=4,  Sm | i=5, Eu | i=6, Gd | i=7, Tb | i=8, Dy |
	  i=9, Ho | i= 10, Er | i=11, Tm | i=12, Yb | i=13, Lu */
       mRComp[0][e]=0.648;
       mRComp[1][e]=1.675;
       mRComp[2][e]=0.254;
       mRComp[3][e]=1.250;
       mRComp[4][e]=0.406;
       mRComp[5][e]=0.154;
       mRComp[6][e]=0.544;
       mRComp[7][e]=0.099;
       mRComp[8][e]=0.674;
       mRComp[9][e]=0.149;
       mRComp[10][e]=0.438;
       mRComp[11][e]=0.068;
       mRComp[12][e]=0.441;
       mRComp[13][e]=0.0675;
    }
    
    // ice
    glacier = 1; // 0 = off, 1 = on.
    icethick = 2e3; // ice thickness of 2 km
    deglaciate = 1; // do some deglaciation...
    markstart = 1; // if equal 0 increase the start deglaciation time
    markend = 1; // if equal 0 increase the end deglaciation time
    deglaciation_start = 9e6; // start deglaciation time (Myr)
    deglaciation_end = deglaciation_start+1e6; // stop deglaciation time (Myr)

    sedipool=calloc(Sbase,sizeof(double)); //sedipool is reset at each beginning and not saved during each save
    //for (s=0;s<Sbase;s++) sedipool[s]=200;
    
    //rho nodes position on the grid

    x=calloc(Sbase,sizeof(double)); x[0]=0; for (s=0;s<Sbase-1;s++) x[s+1]=x[s]+dx[s];
    z=calloc(Tbase,sizeof(double)); z[0]=0; for (t=0;t<Tbase-1;t++) z[t+1]=z[t]+dz[t];
    xmin=x[0]; xmax=x[Sbase-1]; zmin=z[0]; zmax=z[Tbase-1];
    //dx and dz on coarse scale are calculated
    for (G=1;G<=multigridlevels;G++){
        S=(Sbase-1)/pow2[G]+1;
        T=(Tbase-1)/pow2[G]+1;
        for (s=0;s<S-1;s++) {
            dx[INdx(s,G)]=dx[INdx(2*s,G-1)]+dx[INdx(2*s+1,G-1)]; }
        for (t=0;t<T-1;t++) {
            dz[INdz(t,G)]=dz[INdz(2*t,G-1)]+dz[INdz(2*t+1,G-1)]; }
    }
    mincellsize=x[Sbase-1]+z[Tbase-1]; //this value will always be reduced in the following loops
    for (s=0;s<Sbase-1;s++) {if (mincellsize>dx[s]) {mincellsize=dx[s]; }}
    for (t=0;t<Tbase-1;t++) {if (mincellsize>dz[t]) {mincellsize=dz[t]; }}
    //assignment of constants
    leftboundvel=0; rightboundvel=0; upboundvel=0; downboundvel=0;
    alphaP=0.1; alphamom=0.5;
    gx=0; gz=9.82;
}



double getmarkerdevstressinv2(int e,double visc,double my) { //second stress invariant squared
    double stressxx,stressxz;
    stressxx=2*visc*mstrainxx[e]*dt*my/(dt*my+visc)+mstressxx[e]*(visc/(dt*my+visc));
    stressxz=2*visc*mstrainxz[e]*dt*my/(dt*my+visc)+mstressxz[e]*visc/(dt*my+visc);
    return secdevinv2(stressxx,stressxz);
}

double nonlinvisc(int e,double n,double Fn,double Ad,double E,double V,double my) { //Fn is F^n
    double P,T;
    double R;
    double visc,sigmaII,lowsigma,highsigma,sigmaold;
    int i;
    P=mP[e]; T=mTemp[e]+273;
    R=8.314472;
    sigmaII=sqrt(secdevinv2(mstressxx[e],mstressxz[e]));
    lowsigma=sigmaII; highsigma=sigmaII;
    
    visc=pow(2/sigmaII,(n-1)/1)*Fn/Ad*exp((E+P*V)/(R*T));
    if (isinf(visc)) visc=1e100; //if the stress is very low, the viscosity may become infinite
    sigmaII=sqrt(getmarkerdevstressinv2(e,visc,my));
    if (sigmaII>highsigma) highsigma=sigmaII;
    if (sigmaII<lowsigma) lowsigma=sigmaII;
    visc=pow(2/sigmaII,(n-1)/1)*Fn/Ad*exp((E+P*V)/(R*T));
    if (isinf(visc)) visc=1e100; //if the stress is very low, the viscosity may become infinite
    sigmaII=sqrt(getmarkerdevstressinv2(e,visc,my));
    if (sigmaII>highsigma) highsigma=sigmaII;
    if (sigmaII<lowsigma) lowsigma=sigmaII;
    
    i=0;
    while (highsigma-lowsigma>0.1*stressprecision && i<100) {
        i++;
        sigmaII=0.5*(highsigma+lowsigma); sigmaold=sigmaII;
        visc=pow(2/sigmaII,(n-1)/1)*Fn/Ad*exp((E+P*V)/(R*T));
        sigmaII=sqrt(getmarkerdevstressinv2(e,visc,my));
        if (sigmaII>sigmaold) lowsigma=sigmaold; else highsigma=sigmaold;
    }
    return visc;
}

double nonlinviscdiff(int e,double n,double Fn,double Ad,double E,double V,double my,double viscdiff) { //Fn is F^n
    double P,T;
    double R;
    double visc,sigmaII,lowsigma,highsigma,sigmaold;
    int i;
    P=mP[e]; T=mTemp[e]+273;
    R=8.314472;
    sigmaII=sqrt(secdevinv2(mstressxx[e],mstressxz[e]));
    lowsigma=sigmaII; highsigma=sigmaII;
    
    visc=1/(1/(pow(2/sigmaII,(n-1)/1)*Fn/Ad*exp((E+P*V)/(R*T)))+1/viscdiff);
    if (isinf(visc)) visc=1e100; //if the stress is very low, the viscosity may become infinite
    sigmaII=sqrt(getmarkerdevstressinv2(e,visc,my));
    if (sigmaII>highsigma) highsigma=sigmaII;
    if (sigmaII<lowsigma) lowsigma=sigmaII;
    visc=1/(1/(pow(2/sigmaII,(n-1)/1)*Fn/Ad*exp((E+P*V)/(R*T)))+1/viscdiff);
    if (isinf(visc)) visc=1e100; //if the stress is very low, the viscosity may become infinite
    sigmaII=sqrt(getmarkerdevstressinv2(e,visc,my));
    if (sigmaII>highsigma) highsigma=sigmaII;
    if (sigmaII<lowsigma) lowsigma=sigmaII;
    
    i=0;
    while (highsigma-lowsigma>0.1*stressprecision && i<100) {
        i++;
        sigmaII=0.5*(highsigma+lowsigma); sigmaold=sigmaII;
        visc=1/(1/(pow(2/sigmaII,(n-1)/1)*Fn/Ad*exp((E+P*V)/(R*T)))+1/viscdiff);
        sigmaII=sqrt(getmarkerdevstressinv2(e,visc,my));
        if (sigmaII>sigmaold) lowsigma=sigmaold; else highsigma=sigmaold;
    }
    return visc;
}

double frictionsoftening(double phi_peak,double phi_stable,double totalstrain) {
    double emax,emin,output;
    // emax=0.5; emin=1e-5;
    emax=0.1; emin=1e-5;
    if (totalstrain>emin) {
        if (totalstrain<emax) output=phi_peak+(totalstrain-emin)*(phi_stable-phi_peak)/(emax-emin);
        else output=phi_stable;
    }
    else
        output=phi_peak;
    return output;
}

double peierl(double T, double inputstrain) {
    double R,sig0,H,A,sig,strain,b,sigmin,sigmax;
    R=8.314472; H=540e3; A=6.3096e-5; sig0=9.1e+9; b=-H/(R*T); //T is now in Kelvin
    sigmin=1; sigmax=0.999*sig0;
    while (sigmax/sigmin>1.001) {
        sig=pow(sigmax*sigmin,0.5);
        strain=A*sig*sig*exp(b*(1-sig/sig0)*(1-sig/sig0));
        if (strain<inputstrain) sigmin=sig; else sigmax=sig;
    }
    if (T>1000 || sig<1e7) sig=1e100;
    return sig;
}


int calculatemarkerparams(int yield) { //if yield is 1, effective yielding viscosities will be calculated
    int e;
    int numberoffirsttimeyieldmarkers=0;
    numberofyieldmarkers=0;
    #pragma omp parallel private(e) 
    {
        double visc, my, rho, rho0, viscmin, viscmax;
        double Ad, E, V, sigmayield, n; //let R be global. implement with sigmayield later
        double factor, highvisc, lowvisc, sigmaII2;
        double sigmaP;
        double depletionhardening,F;
        int plastyield, peierlyield;
        viscmin=1e-15; viscmax=1e80; //change!
        
        #pragma omp for schedule(static) reduction(+:numberoffirsttimeyieldmarkers,numberofyieldmarkers)
        for (e=0;e<N;e++) {
            sigmayield=1e100; sigmaP=1e100;
            switch(mtype[e]) {
                case 0:
                    visc=1e+21;
                    my=1e+10;
                    rho=1;
                    sigmaP=1e5;
                    break;
	        case 16: /* ice */
	            visc=1e+25;
	            my=1e+10;
	            rho=1000;
	            sigmaP=1e5;
	            break;
                case 1: case 11:
                    visc=1e+24;
                    my=2.5e+10;
		    // visc=nonlinvisc(e, 4, pow(0.5, 4), 1.00e-028, 223000, 0, my); //upper crust Miguel
		    visc=nonlinvisc(e, 3.3, pow(0.5, 3.3), 3.16e-026, 190000, 0, my); //upper crust (Schmalholz-2009).
                    // visc=nonlinvisc(e, 3.2, pow(0.5, 3.2), 2.0822e-023, 238000, 0, my);
		    // visc=nonlinvisc(e, 2.3, pow(2, 1-2*2.3), 5.0119e-018, 154000, 0, my);//ucrust0
                    rho0=2800;
		    // rho=2800;
                    rho=rho0-3.28e-5*mTemp[e]*rho0;
		    sigmayield=10.0e+6+sin(0.0175*frictionsoftening(30, 15, mfinitestrain[e]))*mP[e];
                    // sigmayield=21.6e+6+sin(0.0175*frictionsoftening(36, 0, mfinitestrain[e]))*mP[e];
                    if (sigmayield<0) sigmayield=10.0e+6;
                    break;
                case 2: case 12:
                    //visc=1e+24;
                    my=2.5e+10;
		    // visc=nonlinvisc(e, 4.2, pow(0.5, 4.2), 8.9125e-022, 445000, 0, my); // lower crust Miguel
		    visc=nonlinvisc(e, 3.0, pow(0.5, 3.0), 3.20e-020, 276000, 0, my); //weak lower crust (Schmalholz-2009).
		    // visc=nonlinvisc(e, 4.7, pow(0.5, 4.7), 1.20e-026, 485000, 0, my); //strong lower crust (Schmalholz-2009).
                    // visc=nonlinvisc(e, 3.2, pow(0.5, 3.2), 2.0822e-023, 238000, 0, my);//lcrust
                    rho0=2800;
                    // rho=2800;
		    rho = rho0-3.28e-5*mTemp[e]*rho0;
                    sigmayield=10.0e+6+sin(0.0175*frictionsoftening(30, 15, mfinitestrain[e]))*mP[e];
		    // sigmayield=21.6e+6+sin(0.0175*frictionsoftening(36, 0, mfinitestrain[e]))*mP[e]; 
                    if (sigmayield<0) sigmayield=10.0e+6;
                    break;
                case 3: case 4: case 7: case 8:
                    my=6.7e+10;
                    F=1.0-1.0/mX[e];
		    depletionhardening=(1.0+29.0*erf(F/0.05))*exp(8.0265*F);//31/30*1.0/(exp(-68.0239*F)+1.0/30*exp(-8.0265*F));
		    // visc=2.238e-009*0.5*exp(375000/(8.314472*(mTemp[e]+273))); // Miguel's diffusion creep 
		    visc=2.6009e+010*exp((300e3+5e-6*mP[e])/(8.314472*(mTemp[e]+273))); // original
		    if (mX[e]>1) visc*=depletionhardening;
                    // visc=nonlinviscdiff(e, 3.5, pow(0.5, 3.5), 2.7542e-016, 530e3, 0, my, visc);
		    visc=nonlinviscdiff(e, 3.5, pow(0.5, 3.5), 2.4169e-016, 540e3, 15e-6, my, visc);
                    rho0=3333;
                    rho=rho0*exp(-2.77e-5*mTemp[e]-0.5*0.97e-8*((mTemp[e]+273)*(mTemp[e]+273)-273*273)-0.32*(1.0/(mTemp[e]+273)-1.0/273));
                    if (rho<2000) rho=2000;
		    // sigmayield=10.0e+6+sin(0.0175*frictionsoftening(30, 15, mfinitestrain[e]))*mP[e];
                    sigmayield=21.6e+6+sin(0.0175*frictionsoftening(36, 0, mfinitestrain[e]))*mP[e];
                    if (sigmayield<0) sigmayield=21.6e+6;
                    rho-=0.0591*(1-1/mX[e])*rho0;
//                     sigmaP=50e6;
                    break;
                case 5: case 15:
                    my=2.5e+10;
                    visc=nonlinvisc(e, 2.3, pow(2, 1-2*2.3), 5.0119e-018, 154000, 0, my);//sedi
                    rho0=2300;
                    rho=rho0-3.28e-5*mTemp[e]*rho0;
                    sigmayield=21.6e+6+sin(0.0175*frictionsoftening(36, 0, mfinitestrain[e]))*mP[e];
                    if (sigmayield<0) sigmayield=21.6e+6;
                    break;
                case 6: case 10:
                default:
                    visc=1e+20;
                    my=1e+10;
                    rho=2800;
            }
            peierlyield=0; plastyield=0;
            sigmaII2=getmarkerdevstressinv2(e, visc, my);
            if (sigmayield*sigmayield<sigmaII2) plastyield=1;
            if (sigmaP*sigmaP<sigmaII2 && sigmaP*sigmaP<sigmayield*sigmayield){ peierlyield=1; plastyield=0; }
            
            if (!plastyield) {
                mfinitestrain[e]=-fabs(mfinitestrain[e]); } //no yielding. finite strain becomes negative to indicate this
            else {
                numberofyieldmarkers++;
                if  (mfinitestrain[e]<=0) {
                    numberoffirsttimeyieldmarkers++;
                    if (yield)          {
                        mfinitestrain[e]=fabs(mfinitestrain[e]); } //the number of first time yielding markers is counted. if yield is activated in the function call, mfinitestrain turns positive
                } }
            if (plastyield && yield) { // yielding if marker is above yield and yield is activated from the function call
                factor=sigmayield/sqrt(secdevinv2(mstressxx[e], mstressxz[e]));
                visc=sigmayield/(2*sqrt(secdevinv2(mstrainxx[e], mstrainxz[e])));
                mstressxx[e]*=factor; mstressxz[e]*=factor;
                
            }
            
            if (peierlyield) {
                factor=sigmaP/sqrt(secdevinv2(mstressxx[e], mstressxz[e]));
                visc=sigmaP/(2*sqrt(secdevinv2(mstrainxx[e], mstrainxz[e])));
                mstressxx[e]*=factor; mstressxz[e]*=factor;
            }
            
            if (visc<viscmin) visc=viscmin;
            if (visc>viscmax) visc=viscmax;
            mvisc[e]=visc; mmy[e]=my; mrho[e]=rho;
        }
    }
    return numberoffirsttimeyieldmarkers;
}

double calcnewoutvelocity() {
    double vnew,slope,intersect;
/*    if (stressintegral<0) {
        if (vold==voldold || Iold==Ioldold) { //change the velocity to make a gradient in the next time step
            if (vold==0) vnew=0.0001*3.17e-8; else vnew=0;
        }
        else {
            slope=(Ioldold-Iold)/(voldold-vold); intersect=Iold-slope*vold; vnew=-intersect/slope; //extrapolating to find a velocity that yields a 0 integral
        }
    }
    else vnew=0; //if tension, no pulling*/
    if (stressintegral<0) vnew=0.0001*3.17e-8; else vnew=0;
    return vnew;
}

void timestep() {
    int accetable_vel,convergence;
    int timestepsbetweenprint,tempstepsbetweensmallprint;
    double olddtnom,dtnew,timestart,timesolve,timeother;
    timestart=omp_get_wtime();
    statusfile = fopen("status.txt","a");
    fprintf(statusfile,"ITERATION: %d\n",iterationnumber); 
    timestepsbetweenprint=2000; tempstepsbetweensmallprint=10;
    // if (iterationnumber%timestepsbetweenprint==0) printstate(iterationnumber/timestepsbetweenprint);
    // if (iterationnumber%tempstepsbetweensmallprint==0) outputsurfaces();
    
    Tp=1400; // Potential temperature in Celcius for bottom boundary condition.
    Tbottom=Tp+273-(710e3-zmax)*0.5/1000;
    olddtnom=dtnom;
    voldold=vold; vold=leftboundvel; Ioldold=Iold; Iold=stressintegral;

    // change ice thickness at these times
    if (time>=deglaciation_start*60*60*24*365 && time<=(deglaciation_start + 1e5)*60*60*24*365) {
      markstart = 0;
      markend = 0;
      timestepsbetweenprint=1;
    }
    if (time>=deglaciation_end*60*60*24*365 && markstart == 0) {
      deglaciation_start = deglaciation_start + 1e6;
      markstart = 1;
    }
    if (time>=(deglaciation_end + 1e5)*60*60*24*365 && markend == 0) {
      deglaciation_end = deglaciation_start + 1e6;
      markend = 1;
    }
    
    if (iterationnumber%timestepsbetweenprint==0) printstate(iterationnumber/timestepsbetweenprint);
    if (iterationnumber%tempstepsbetweensmallprint==0) outputsurfaces();
        
    riftvelo=0.010; //0.015; // m/yr
    if (xmax<3300e3 && time>1e5*60*60*24*365) {
        velocityboundarychange1(0,1,0,-1,riftvelo*3.17e-8);
//         if (riftvelo>0.005) dtnom*=0.005/riftvelo;
    }
    else {
        velocityboundarychange1(0,1,0,-1,0);
    }

    accetable_vel=0;
    while (accetable_vel==0) {
        timeother=omp_get_wtime();
        markermomentumparamstogrid(); fprintf(statusfile,"momentum interp time %g\n",omp_get_wtime()-timeother);
        fprintf(statusfile,"time:%g Ma. dt: %g Ma. dtnom %g Ma dtmaxwell %g Ma\n",time/3.1536e13,dt/3.1536e13,dtnom/3.1536e13,dtmaxwell/3.1536e13); 
        timesolve=omp_get_wtime();
        convergence=solve(stressprecision);  fprintf(statusfile,"solve time %g\n",omp_get_wtime()-timesolve);
//         fprintf(statusfile,"max vel %g m/s, dt %g Ma, mincellsize %g m\n",maxvelocity(),dt,mincellsize);
        if (dt<5*mincellsize/maxvelocity()) accetable_vel=1;
        if (!convergence) {accetable_vel=0; dt*=0.5; }
    }
    timeother=omp_get_wtime(); stressstrainupdate(); fprintf(statusfile,"stressstrainupdate time %g\n",omp_get_wtime()-timeother);
    timeother=omp_get_wtime(); updatetemp(); fprintf(statusfile,"updatetemp time %g ",omp_get_wtime()-timeother);
    if (iterationnumber%tempstepsbetweensmallprint==0) { outputtempprofiles(); outputmeltcomp(); }
    meltupdate();
    mocompinst();
    timeother=omp_get_wtime(); stressrotation(); fprintf(statusfile,"stressrotation time %g\n",omp_get_wtime()-timeother);
    timeother=omp_get_wtime(); movemarkers1(); fprintf(statusfile, "movemarkers time %g\n", omp_get_wtime()-timeother);
    
    timeother=omp_get_wtime();

    if (iterationnumber%1000==0) surfaceupdate(1); else surfaceupdate(0);
    fprintf(statusfile, "surfaceupdate time %g\n", omp_get_wtime()-timeother);
    
    time+=dt; iterationnumber++;
    fprintf(statusfile,"number of yield markers %d\n",numberofyieldmarkers); 
    fprintf(statusfile,"total time %g\n\n",omp_get_wtime()-timestart);
    fclose(statusfile);
    
    dtnom=olddtnom;
}

void timestep_force() {
    int accetable_vel;
    int convergence;
    int timestepsbetweenprint,tempstepsbetweensmallprint;
    double olddtnom,dtnew,timestart,timesolve,timeother;

    const double rift_step = 0.001; // m/yr, how much to increase each trial
    const double stress_threshold = 10e12; // maximum stress
    const int    max_rift_iters = 100;    // safety limit
    double stress_integral;
    double stress_tollerance = 0.01; // how close to zero the stress integral needs to be to accept the rift velocity

    timestart=omp_get_wtime();
    statusfile = fopen("status.txt","a");
    fprintf(statusfile,"ITERATION: %d\n",iterationnumber); 
    timestepsbetweenprint=2000; tempstepsbetweensmallprint=10;
    // if (iterationnumber%timestepsbetweenprint==0) printstate(iterationnumber/timestepsbetweenprint);
    // if (iterationnumber%tempstepsbetweensmallprint==0) outputsurfaces();
    
    Tp=1400; // Potential temperature in Celcius for bottom boundary condition.
    Tbottom=Tp+273-(710e3-zmax)*0.5/1000;
    olddtnom=dtnom;
    voldold=vold; vold=leftboundvel; Ioldold=Iold; Iold=stressintegral;

    // change ice thickness at these times
    if (time>=deglaciation_start*60*60*24*365 && time<=(deglaciation_start + 1e5)*60*60*24*365) {
      markstart = 0;
      markend = 0;
      timestepsbetweenprint=1;
    }
    if (time>=deglaciation_end*60*60*24*365 && markstart == 0) {
      deglaciation_start = deglaciation_start + 1e6;
      markstart = 1;
    }
    if (time>=(deglaciation_end + 1e5)*60*60*24*365 && markend == 0) {
      deglaciation_end = deglaciation_start + 1e6;
      markend = 1;
    }
    
    if (iterationnumber%timestepsbetweenprint==0) printstate(iterationnumber/timestepsbetweenprint);
    if (iterationnumber%tempstepsbetweensmallprint==0) outputsurfaces();
    
    if (iterationnumber == 0) {
        riftvelo = 0.001; // m/yr, initial guess for rift velocity
    }

    if (time <= 1e5*60*60*24*365) {
        /* Before rifting starts, just run the model with zero velocity and no adjustments */
        velocityboundarychange1(0, 1, 0, -1, 0);

        accetable_vel = 0;
        while (accetable_vel==0) {
            timeother = omp_get_wtime();
            markermomentumparamstogrid();
            fprintf(statusfile,"momentum interp time %g\n",omp_get_wtime()-timeother);
            fprintf(statusfile,"time:%g Ma. dt: %g Ma. dtnom %g Ma dtmaxwell %g Ma\n",
                time/3.1536e13,dt/3.1536e13,dtnom/3.1536e13,dtmaxwell/3.1536e13); 
            timesolve = omp_get_wtime();
            convergence = solve(stressprecision);
            fprintf(statusfile,"solve time %g\n",omp_get_wtime()-timesolve);

            if (dt<5*mincellsize/maxvelocity()) accetable_vel=1;
            if (!convergence) {
                accetable_vel=0;
                dt*=0.5; 
            }
        }
    }
    else {
        /* adapt riftvelo */
        fprintf(statusfile, "Starting rift velocity adaptation loop...\n");

        double riftvelo_trial = riftvelo;   // start from last accepted value
        int    rift_iter = 0;
        convergence = 0; // initialize convergence to false

        while (rift_iter < max_rift_iters) {

            // Update rift boundary velocity for this trial
            velocityboundarychange1(0, 1, 0, -1, riftvelo_trial * 3.17e-8);

            // Solve with this riftvelo_trial
            accetable_vel=0;
            while (accetable_vel==0) {
                timeother=omp_get_wtime();
                markermomentumparamstogrid();
                fprintf(statusfile,"momentum interp time %g\n",omp_get_wtime()-timeother);
                fprintf(statusfile,
                    "time:%g Ma. dt: %g Ma. dtnom %g Ma dtmaxwell %g Ma\n",
                    time/3.1536e13,dt/3.1536e13,dtnom/3.1536e13,dtmaxwell/3.1536e13);

                timesolve=omp_get_wtime();
                convergence=solve(stressprecision);
                fprintf(statusfile,"solve time %g\n",omp_get_wtime()-timesolve);
                // fprintf(statusfile,"max vel %g m/s, dt %g Ma, mincellsize %g m\n",maxvelocity(),dt,mincellsize);

                if (!convergence) {
                    // Solver failed: reduce dt and restart whole outer loop
                    dt *= 0.5;
                    accetable_vel = 0;
                }
                if (dt<5*mincellsize/maxvelocity()) accetable_vel=1;
            }

            // get stress integral using read only function that doesn't change the model state
            stress_integral = stressstrainupdate_readonly();
                
            // Print stress integral and check if it's acceptable
            fprintf(statusfile, "Rift trial %d: riftvelo=%g m/yr, stress_integral=%g\n",
                rift_iter, riftvelo_trial, stress_integral);
                
            // Check threshold
            if (stress_integral >= (1 - stress_tollerance) * stress_threshold && stress_integral <= (1 + stress_tollerance) * stress_threshold) {
                // Accept this riftvelo_trial
                riftvelo = riftvelo_trial;
                break;
            }
            else if (stress_integral > (1 + stress_tollerance) * stress_threshold) {
                // Stress too high, decrease riftvelo
                riftvelo_trial -= rift_step;
            }
            else {
                // Stress too low, increase riftvelo
                riftvelo_trial += rift_step;
            }

            // Stop the model if the number of itterations is too high, to avoid infinite loops
            if (rift_iter == max_rift_iters - 1) {
                fprintf(statusfile, "Warning: Maximum rift velocity iterations reached without \
                    achieving the threshold. Continuing with final force = %g.\n", stress_integral);
                riftvelo = riftvelo_trial;
            }

            ++rift_iter;
        }
    }
    
    /* update stresses and strains */
    timeother=omp_get_wtime();
    stressstrainupdate();
    fprintf(statusfile,"stressstrainupdate time %g\n",omp_get_wtime()-timeother);

    timeother=omp_get_wtime();
    updatetemp();
    fprintf(statusfile,"updatetemp time %g ",omp_get_wtime()-timeother);
    if (iterationnumber%tempstepsbetweensmallprint==0) { 
        outputtempprofiles();
        outputmeltcomp();
    }
    meltupdate();
    mocompinst();
    timeother=omp_get_wtime();
    stressrotation();
    fprintf(statusfile,"stressrotation time %g\n",omp_get_wtime()-timeother);
    timeother=omp_get_wtime();
    movemarkers1();
    fprintf(statusfile, "movemarkers time %g\n", omp_get_wtime()-timeother);
    
    timeother=omp_get_wtime();
    if (iterationnumber%1000==0) surfaceupdate(1);
    else surfaceupdate(0);
    fprintf(statusfile, "surfaceupdate time %g\n", omp_get_wtime()-timeother);
    
    time+=dt; iterationnumber++;
    fprintf(statusfile,"number of yield markers %d\n",numberofyieldmarkers); 
    fprintf(statusfile,"total time %g\n\n",omp_get_wtime()-timestart);
    fclose(statusfile);
    
    dtnom=olddtnom;
}

void firsttimesteps() { //number of short timesteps to make elastic stress buildups
    double olddtnom;
    int a;
    olddtnom=dtnom; //backup of nominal timestep
     convertcenterpoints=1;
    for (dtnom=0.0000001*olddtnom;dtnom<olddtnom;dtnom*=2){
        for (a=0;a<1;a++) {
            timestep();
            //dtnom=olddtnom;
        }
    }
    dtnom=olddtnom;
}

int main(int argc, char *argv[] ) {
    int a,gravityon;
    double lengthx,lengthz;
    int riftinghasbegun,plottask,NNN;
    double riftendtime;
    char mod_type;

    // Check that we have enough arguments
    if (argc < 3) {
        fprintf(stderr, "Usage: %s <num_threads> <V|F>\n", argv[0]);
        return 1;
    }

    sscanf(argv[1], "%d", &Nthreads);
    mod_type = argv[2][0];  // Read first character of second argument

    // Validate the mode type
    if (mod_type != 'V' && mod_type != 'F') {
        fprintf(stderr, "Error: mod_type must be 'V' or 'F', got '%c'\n", mod_type);
        return 1;
    }

    omp_set_num_threads(Nthreads);

    statusfile = fopen("status.txt","w");
    surfacefile = fopen("outfiles/surfaces.txt","w");
    fclose(surfacefile);
        
    initialize();
    stressprecision=1e4;
    
    #pragma omp parallel
    {
    NNN=omp_get_num_threads();
    printf("NNN: %d\n",NNN);
    }

    Ttop=0;
    
    iterationnumber=0;
    gz=9.82;
    //gx=9.82;
    //dtnom=1e+3*60*60*24*365;
    fclose(statusfile);
    
//     e=1e5;
//     printf("e: %d, mxd=%f mzd=%f\n",e,mx[e],mz[e]);
//     printf("nearest marker of %d is %d\n",e,getnearestmarker(mx[e],mz[e]));
    
    dtnom=1e2*60*60*24*365; 
    firsttimesteps();
    riftinghasbegun=0; riftendtime=1e80; plottask=0;
    for (a=0;time<40e6*60*60*24*365;a++) {
    //for (a=0;a=2;a++) { // speed testing -> 3 time steps
        if (iterationnumber>100000000) break;
        if (mod_type == 'V') {
            timestep();
        } else if (mod_type == 'F') {
            timestep_force();
        }
    }
    printstate(20000);
    return 0;
}
