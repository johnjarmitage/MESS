#include <stdio.h>
#include <math.h>
#include <omp.h>
#include <stdlib.h>
#include "main.h"
#include "tempfuncs_par.h"

void outputsurfaces() {
    FILE *file0,*file1,*file2,*file3,*file4;
    char *writemode;
    if (iterationnumber==0) writemode="wb"; else writemode="ab";
    
    file0 = fopen("outfiles/xmax.txt",writemode);
    file1 = fopen("outfiles/times.txt",writemode);
    file2 = fopen("outfiles/surface.txt",writemode);
    file3 = fopen("outfiles/basement.txt",writemode);
    file4 = fopen("outfiles/moho.txt",writemode);
    
//     fwrite(&time, sizeof(double), 1, file1);
//     for (s=0;s<Sbase;s++) {
//         fwrite(&surfaces[INsur(s, Nsurfaces-1)], sizeof(double), 1, file2);
//         fwrite(&surfaces[INsur(s, 3)], sizeof(double), 1, file3);
//         fwrite(&surfaces[INsur(s, 2)], sizeof(double), 1, file4);
//     }
    
//     fwrite(&time, sizeof(double), 1, file1);
//     for (s=0;s<Sbase;s++) {
//         a=(double)INsur(s, Nsurfaces-1);
//         b=(double)INsur(s, 3);
//         c=(double)INsur(s, 2);
//         fwrite(&a, sizeof(double), 1, file2);
//         fwrite(&b, sizeof(double), 1, file3);
//         fwrite(&c, sizeof(double), 1, file4);
//     }
//     
    fwrite(&x[Sbase-1], sizeof(double), 1, file0);
    fwrite(&time, sizeof(double), 1, file1);
    fwrite(surfaces+INsur(0, Nsurfaces-1), sizeof(double), Sbase, file2);
    fwrite(surfaces+INsur(0,3), sizeof(double), Sbase, file3);
    fwrite(surfaces+INsur(0,1), sizeof(double), Sbase, file4);

    fclose(file0); fclose(file1); fclose(file2); fclose(file3); fclose(file4);
        
    }
    

void landscapediff(double *h,double *old_h,double k,double dt) {
    double D;
    double *a,*b,*c,*d,*bdot,*ddot;
    
    a = calloc(Sbase,sizeof(double));
    b = calloc(Sbase,sizeof(double));
    c = calloc(Sbase,sizeof(double));
    d = calloc(Sbase,sizeof(double));
    bdot = calloc(Sbase,sizeof(double));
    ddot = calloc(Sbase,sizeof(double));
    
    for (s=1;s<Sbase-1;s++) {
        D=dx[s]+dx[s-1];
        a[s]=2*k/(D*dx[s-1]);
        b[s]=-(2*k/(D*dx[s-1])+2*k/(D*dx[s])+1/dt);
        c[s]=2*k/(D*dx[s]);
        d[s]=-old_h[s]/dt;
    }
    d[0]=0; d[Sbase-1]=0;
    b[0]=1; c[0]=-1;
    b[Sbase-1]=1; a[Sbase-1]=-1;
    tridiag(h,a,b,c,d,bdot,ddot,Sbase);
    
    free(a);
    free(b);
    free(c);
    free(d);
    free(bdot);
    free(ddot);
}
//as above, but variable diffusivity
void landscapediff1(double *h,double *old_h,double *k,double dt) {
    double alphaplus,alphaminus;
    double *a,*b,*c,*d,*bdot,*ddot;
    
    a = calloc(Sbase,sizeof(double));
    b = calloc(Sbase,sizeof(double));
    c = calloc(Sbase,sizeof(double));
    d = calloc(Sbase,sizeof(double));
    bdot = calloc(Sbase,sizeof(double));
    ddot = calloc(Sbase,sizeof(double));
    
    for (s=1;s<Sbase-1;s++) {
        alphaplus=1/(dx[s]*(dx[s-1]+dx[s]))*(k[s]+k[s+1]);
        alphaminus=1/(dx[s-1]*(dx[s-1]+dx[s]))*(k[s]+k[s-1]);
        a[s]=alphaminus;
        b[s]=-(alphaminus+alphaplus+1/dt);
        c[s]=alphaplus;
        d[s]=-old_h[s]/dt;
    }
    
    d[0]=0; d[Sbase-1]=0;
    b[0]=1; c[0]=-1;
    b[Sbase-1]=1; a[Sbase-1]=-1;
    tridiag(h,a,b,c,d,bdot,ddot,Sbase);
    
    free(a);
    free(b);
    free(c);
    free(d);
    free(bdot);
    free(ddot);
}

void sedipooldiff(double *h,double *old_h,double *newtopo,double sealevel,double kland,double ksea,double leftbasinlimit,double rightbasinlimit,double dt) {
    double kplus,kcenter,kminus,alpha;
    double *a,*b,*c,*d,*bdot,*ddot;
    
    a = calloc(Sbase,sizeof(double));
    b = calloc(Sbase,sizeof(double));
    c = calloc(Sbase,sizeof(double));
    d = calloc(Sbase,sizeof(double));
    bdot = calloc(Sbase,sizeof(double));
    ddot = calloc(Sbase,sizeof(double));
    
    for (s=1; s<Sbase-1; s++) {
        if (newtopo[s+1]>sealevel && x[s+1]<rightbasinlimit && x[s+1]>leftbasinlimit) kplus=ksea; else kplus=kland;
        if (newtopo[s]>sealevel && x[s]<rightbasinlimit && x[s]>leftbasinlimit) kcenter=ksea; else kcenter=kland;
        if (newtopo[s-1]>sealevel && x[s-1]<rightbasinlimit && x[s-1]>leftbasinlimit) kminus=ksea; else kminus=kland;
        alpha=1/(dx[s-1]+dx[s]);
        a[s]=alpha*(kminus+kcenter)/dx[s-1];
        b[s]=-(1/dt+alpha*(kminus+kcenter)/dx[s-1]+alpha*(kplus+kcenter)/dx[s]);
        c[s]=alpha*(kplus+kcenter)/dx[s];
        d[s]=-old_h[s]/dt;
        }
    d[0]=0; d[Sbase-1]=0;
    b[0]=1; c[0]=-1;
    b[Sbase-1]=1; a[Sbase-1]=-1;
    tridiag(h,a,b,c,d,bdot,ddot,Sbase);
    
    free(a);
    free(b);
    free(c);
    free(d);
    free(bdot);
    free(ddot);
}

void updatesedipool(double *sedipool,double *newtopo,double sealevel,double leftbasinlimit,double rightbasinlimit,double decaylength,double rate){ //rate is the maximum available rate in m/year
    int leftcoast,rightcoast;
    rate=rate/(60*60*24*365)*dt;
    for (s=0; s<Sbase; s++) if (x[s]<=rightbasinlimit) leftcoast=s;
    for (s=Sbase-1; s>-1; s--) if (x[s]>=leftbasinlimit) rightcoast=s;
 
    for (s=0; s<Sbase; s++) {
        if (newtopo[s]>sealevel && x[s]>leftbasinlimit && s<leftcoast) leftcoast=s;
        if (newtopo[s]>sealevel && x[s]<rightbasinlimit && s>rightcoast) rightcoast=s;
    }
    //printf("%d %d\n",leftcoast,rightcoast);
    for (s=0; s<Sbase; s++) {
        sedipool[s]=0;
        if (s>=leftcoast) sedipool[s]+=rate*exp(-(x[s]-x[leftcoast])/decaylength);
        if (s<=rightcoast) sedipool[s]+=rate*exp((x[s]-x[rightcoast])/decaylength);
        //printf("%f",sedipool[s]);
    }
    //printf("\n");
}

void printsurface() {
    double localsealevel;
    surfacefile = fopen("outfiles/surfaces.txt","a");
    fprintf(surfacefile,"%d ",Nsurfaces);
    fprintf(surfacefile,"%f ",time);
    for (s=0;s<Sbase;s++) fprintf(surfacefile,"%f ",surfaces[INsur(s,Nsurfaces-1)]-sealevel);
    fprintf(surfacefile,"\n");
    fclose(surfacefile);
}

void surfaceupdate(int newsurface) {
    double *newtopo,*oldtopo,*oldsedipool,*k_var;
    double k;
    double *oldsurfaces;
    double topo;
    double slope, interc;
    double sedipoolincrement;
    double leftbasinlimit,rightbasinlimit;
    int ice_on_top;

    ice_on_top = 2;
    
    leftbasinlimit=0; rightbasinlimit=x[S-1];
    sealevel=0.5*(surfaces[INsur(0,3)]+surfaces[INsur(Sbase-1,3)]); // sealevel is constant after rifting
    surfaces[INsur(0,2)]=sealevel; //this particular space is not in use, and saves the sealevel
    
    
    if (newsurface) { //create a new surface, allocate new memory and copy old surfaces
        oldsurfaces=surfaces;
        surfaces=calloc((Nsurfaces+1)*Sbase,sizeof(double));
        for (e=0;e<Nsurfaces*Sbase;e++) surfaces[e]=oldsurfaces[e];
        Nsurfaces++;
        free(oldsurfaces);
        for (s=0;s<Sbase;s++) surfaces[INsur(s,Nsurfaces-1)]=surfaces[INsur(s,Nsurfaces-2)]; //copying the old top surfaces to the new top surface
    }
    
    newtopo=calloc(Sbase, sizeof(double));
    oldtopo=calloc(Sbase, sizeof(double));
    oldsedipool=calloc(Sbase, sizeof(double));
    k_var=calloc(Sbase, sizeof(double));
    
    
    for (s=0;s<Sbase;s++) {
        oldtopo[s]=surfaces[INsur(s, Nsurfaces-1)];
        oldsedipool[s]=sedipool[s]; }
    
 
    
    k=6.6e-6;
    for (s=0;s<Sbase;s++) k_var[s]=k;
    landscapediff1(newtopo,oldtopo,k_var,dt);
    
//     landscapediff(newtopo,oldtopo,k,dt);
    //sedipooldiff(sedipool,oldsedipool,newtopo,sealevel,10,1e-7,leftbasinlimit,rightbasinlimit,dt);
    //updatesedipool(sedipool,newtopo,sealevel,leftbasinlimit,rightbasinlimit,10e3,0.00040);
    
    //sealevel=0.5*(surfaces[INsur(0,Nsurfaces-1)]+surfaces[INsur(Sbase-1,Nsurfaces-1)]);
    //sealevel=0.5*(surfaces[INsur(0,3)]+surfaces[INsur(Sbase-1,3)]);
    

    
    for (s=0;s<Sbase;s++) {
//         if (newtopo[s]>surfaces[INsur(s,Nsurfaces-1)] && surfaces[INsur(s,Nsurfaces-1)]<sealevel) { //erosion only above sealevel
        if (newtopo[s]>surfaces[INsur(s, Nsurfaces-1)]) { //erosion everywhere
            //sedipool[s]+=newtopo[s]-surfaces[INsur(s,Nsurfaces-1)];
            surfaces[INsur(s,Nsurfaces-1)]=newtopo[s]; //erosion criterium
            for (sur=Nsurfaces-2;sur>=0;sur--) {
                if (surfaces[INsur(s,Nsurfaces-1)]>surfaces[INsur(s,sur)]) surfaces[INsur(s,sur)]=surfaces[INsur(s,Nsurfaces-1)];
                else break;
            }
        }
        else surfaces[INsur(s,Nsurfaces-1)]=newtopo[s]; //sedimentation
    }
//     for (s=0;s<Sbase;s++) {
//         if (surfaces[INsur(s,Nsurfaces-1)]>sealevel) {  //sedimentation criterium
//             surfaces[INsur(s,Nsurfaces-1)]-=sedipool[s]; 
//             sedipool[s]=0; //emptying sedipool
//             if (surfaces[INsur(s,Nsurfaces-1)]<sealevel) { //if overflow - return surplus to sedipool
//                 sedipool[s]=sealevel-surfaces[INsur(s,Nsurfaces-1)]; 
//                 surfaces[INsur(s,Nsurfaces-1)]=sealevel; }}}
//     
    //updating markers according to the new surfaces
    for (e=0;e<N;e++) {
        s=mcs[e]; t=mct[e];
        slope=(surfaces[INsur(s+1,Nsurfaces-1)]-surfaces[INsur(s,Nsurfaces-1)])/(x[s+1]-x[s]);
        interc=surfaces[INsur(s,Nsurfaces-1)]-slope*x[s];
        topo=slope*mx[e]+interc; //local topography
	if (glacier == 0) {// no ice
	  ice_on_top = 0;
	  if (topo<mz[e]) { //possible sedimentation
            if (mtype[e]==0) { //check for air to sediment conversion
                mtype[e]=5; //WARNING: seditype may vary
                mstrainxx[e]=1e-20; mstrainxz[e]=1e-20; //decrease strain for the next time step
                //remember age
                mage[e]=time;
            }
	  }
	  else { //possible erosion
            if (mtype[e]>0) {
                mtype[e]=0; //mstressxx[e]=0; mstressxz[e]=0; no erosion!
            }
	  }
	}
	if (glacier == 1) { // ice
	  // fprintf(statusfile,"Ice layer!\n");
	  if (topo<mz[e]) { //possible sedimentation
            if (mtype[e]==0 || mtype[e]==16) { //check for air/ice to sediment conversion
                mtype[e]=5; //WARNING: seditype may vary
                mstrainxx[e]=1e-20; mstrainxz[e]=1e-20; //decrease strain for the next time step
                //remember age
                mage[e]=time;
            }
	  }
	  if (topo-icethick<mz[e] && mz[e]<=topo) { // put ice on top
	    if (deglaciate == 0){
	      mtype[e] = 16;
	      mstrainxx[e]=1e-20; mstrainxz[e]=1e-20; //decrease strain for the next time step
	      ice_on_top = 1;
	    }
	    if (deglaciate == 1){
	      if (time>=deglaciation_start*60*60*24*365) { // deglaciate
		mtype[e] = 0;
		ice_on_top = 0;
	      }
	      else { // glaciate
		mtype[e] = 16;
		mstrainxx[e]=1e-20;
		mstrainxz[e]=1e-20;
		ice_on_top = 1;
	      }
	    }
	  }
	  if (mz[e]<=topo-icethick) // then air
	    mtype[e] = 0;
	}
    }
    fprintf(statusfile,"Iteration %d deglaciate %d deglaciation_start %g deglaciation_end %g ice on top? %d\n",iterationnumber,deglaciate,deglaciation_start,deglaciation_end,ice_on_top);




    
    if (newsurface) printsurface(); //prints the state of the uppermost surface
    
    free(newtopo);
    free(oldtopo);
    free(oldsedipool);
    free(k_var);
}
