#include <stdio.h>
#include <math.h>
#include <omp.h>
#include <stdlib.h>
#include "main.h"

void iniP() {
    for (s=0;s<Sbase-1;s++){
        P[IN(s,0,0)]=0;
        for (t=1;t<Tbase-1;t++)	{
            P[IN(s,t,0)]=P[IN(s,t-1,0)]+0.5*(gridrho[IN(s,t,0)]+gridrho[IN(s+1,t,0)])*9.82*(dz[INdz(t-1,0)]+dz[INdz(t,0)])*0.5;
        } }
}



void correct(int G) {
    int S, T, st,s,t;
    S=(Sbase-1)/pow2[G]+1;
    T=(Tbase-1)/pow2[G]+1;
    //#pragma omp parallel for schedule(static,65) default(shared) private(s,t,st)
    for (t=0;t<T;t++) {
        for (s=0;s<S;s++) {
            st=IN(s, t, G);
            if (s==S-1 || t==T-1 || s==0 && (t==0 || t==T-2) || s==S-2 && (t==0 || t==T-2))	{ }
            else {
                P[st]+=resP[st]*alphacorrect; }
            if (s==0 || s==S-1 || t==0 || t==T-2 || t==T-1) { }
            else {
                vx[st]+=resx[st]*alphacorrect; }
            if (t==0 || t==T-1 || s==0 || s==S-2 || s==S-1) { }
            else {
                vz[st]+=resz[st]*alphacorrect; }
        }
    }
}

void coarsenrhogrid(double *fine,double *coarse,int G) {
    double dzN,dzS,dxV,dxE,dztot,dxtot,dw,dimlx,dimlz;
    int S,T,Snew,Tnew,st,s,t;
    S=(Sbase-1)/pow2[G]+1;
    T=(Tbase-1)/pow2[G]+1;
    Snew=(S-1)/2+1;
    Tnew=(T-1)/2+1;
    //#pragma omp parallel for schedule(static,65) default(shared) private(s,t)
    for (t=0;t<Tnew;t++) { //delete old values
        for (s=0;s<Snew;s++) {
            coarse[IN(s,t,G+1)]=0;
        } }
    //#pragma omp parallel for schedule(static,65) default(shared) private(s,t,dzN,dzS,dxV,dxE,dztot,dxtot,dw,dimlx,dimlz)
    for (t=0;t<Tnew-1;t++) {
        for (s=0;s<Snew-1;s++) { //the coarse[Snew-1] and coarse[Tnew-2] gets assigned values but these will not be used
            dxV=dx[INdx(2*s,G)];
            dxE=dx[INdx(2*s+1,G)];
            dzN=dz[INdz(2*t,G)];
            dzS=dz[INdz(2*t+1,G)];
            dxtot=dxV+dxE;
            dztot=dzN+dzS;
            
            //fine NV - only contributes to coarse NV with dw=1
            dw=1;
            coarse[IN(s,t,G+1)]+=dw*fine[IN(2*s,2*t,G)];
            weight[IN(s,t,G+1)]+=dw;
            //fine NE - only contribution to coarse NV and NE
            dimlx=dxV/dxtot;
            dimlz=0;
            dw=(1-dimlx)*(1-dimlz);
            coarse[IN(s,t,G+1)]+=dw*fine[IN(2*s+1,2*t,G)];
            weight[IN(s,t,G+1)]+=dw;
            dw=dimlx*(1-dimlz);
            coarse[IN(s+1,t,G+1)]+=dw*fine[IN(2*s+1,2*t,G)];
            weight[IN(s+1,t,G+1)]+=dw;
            //fine SV - only contributes to coarse NV and SV
            dimlx=0;
            dimlz=dzN/dztot;
            dw=(1-dimlx)*(1-dimlz);
            coarse[IN(s,t,G+1)]+=dw*fine[IN(2*s,2*t+1,G)];
            weight[IN(s,t,G+1)]+=dw;
            dw=(1-dimlx)*dimlz;
            coarse[IN(s,t+1,G+1)]+=dw*fine[IN(2*s,2*t+1,G)];
            weight[IN(s,t+1,G+1)]+=dw;
            //fine SE - contributions to all 4 coarse nodes
            dimlx=dxV/dxtot;
            dimlz=dzN/dztot;
            dw=(1-dimlx)*(1-dimlz);
            coarse[IN(s,t,G+1)]+=dw*fine[IN(2*s+1,2*t+1,G)];
            weight[IN(s,t,G+1)]+=dw;
            dw=dimlx*(1-dimlz);
            coarse[IN(s+1,t,G+1)]+=dw*fine[IN(2*s+1,2*t+1,G)];
            weight[IN(s+1,t,G+1)]+=dw;
            dw=(1-dimlx)*dimlz;
            coarse[IN(s,t+1,G+1)]+=dw*fine[IN(2*s+1,2*t+1,G)];
            weight[IN(s,t+1,G+1)]+=dw;
            dw=dimlx*dimlz;
            coarse[IN(s+1,t+1,G+1)]+=dw*fine[IN(2*s+1,2*t+1,G)];
            weight[IN(s+1,t+1,G+1)]+=dw;
        }
    }
    
    t=Tnew-1;
    for (s=0;s<Snew-1;s++) {
        dxV=dx[INdx(2*s,G)];
        dxE=dx[INdx(2*s+1,G)];
        dxtot=dxV+dxE;
        //fine NV - only contributes to coarse NV with dw=1
        dw=1;
        coarse[IN(s,t,G+1)]+=dw*fine[IN(2*s,2*t,G)];
        weight[IN(s,t,G+1)]+=dw;
        //fine NE - only contribution to coarse NV and NE
        dimlx=dxV/dxtot;
        dimlz=0;
        dw=(1-dimlx)*(1-dimlz);
        coarse[IN(s,t,G+1)]+=dw*fine[IN(2*s+1,2*t,G)];
        weight[IN(s,t,G+1)]+=dw;
        dw=dimlx*(1-dimlz);
        coarse[IN(s+1,t,G+1)]+=dw*fine[IN(2*s+1,2*t,G)];
        weight[IN(s+1,t,G+1)]+=dw; }
    
    s=Snew-1;
    for (t=0;t<Tnew-1;t++) {
        dzN=dz[INdz(2*t,G)];
        dzS=dz[INdz(2*t+1,G)];
        dztot=dzN+dzS;
        //fine NV - only contributes to coarse NV with dw=1
        dw=1;
        coarse[IN(s,t,G+1)]+=dw*fine[IN(2*s,2*t,G)];
        weight[IN(s,t,G+1)]+=dw;
        //fine SV - only contributes to coarse NV and SV
        dimlx=0;
        dimlz=dzN/dztot;
        dw=(1-dimlx)*(1-dimlz);
        coarse[IN(s,t,G+1)]+=dw*fine[IN(2*s,2*t+1,G)];
        weight[IN(s,t,G+1)]+=dw;
        dw=(1-dimlx)*dimlz;
        coarse[IN(s,t+1,G+1)]+=dw*fine[IN(2*s,2*t+1,G)];
        weight[IN(s,t+1,G+1)]+=dw; }
    
    s=Snew-1; //interpolation of most SE-point
    t=Tnew-1;
    dw=1;
    coarse[IN(s,t,G+1)]+=dw*fine[IN(2*s,2*t,G)];
    weight[IN(s,t,G+1)]+=dw;
    
    
    //#pragma omp parallel for schedule(static,65) default(shared) private(s,t,st)
    for (t=0;t<Tnew;t++) {
        for (s=0;s<Snew;s++) {
            st=IN(s,t,G+1);
            if (weight[st]>0) {
                coarse[st]=coarse[st]/weight[st];
            }
            weight[st]=0;
        }
    }
}

void finePgrid(double *coarse,double *fine,int G) {
    double dzNN,dzN,dzS,dzSS,dxVV,dxV,dxE,dxEE,dztot,dxtot,dw,dimlx,dimlz;
    int S,T,Snew,Tnew,s,t,st;
    int INfNV,INfNE,INfSV,INfSE;
    int INcNV,INcNE,INcSV,INcSE;
    S=(Sbase-1)/pow2[G]+1;
    T=(Tbase-1)/pow2[G]+1;
    Snew=(S-1)*2+1;
    Tnew=(T-1)*2+1;
    //#pragma omp parallel for schedule(static,65) default(shared) private(s,t)
    for (t=0;t<Tnew;t++) { //deleting old values
        for (s=0;s<Snew;s++) {
            fine[IN(s,t,G-1)]=0; } }
    //#pragma omp parallel for schedule(static,65) default(shared) private(s,t,dzNN,dzN,dzS,dzSS,dxVV,dxV,dxE,dxEE,dztot,dxtot,dw,dimlx,dimlz,INfNV,INfNE,INfSV,INfSE,INcNV,INcNE,INcSV,INcSE)
    for (t=0;t<T-2;t++) { //Loops skip grid[highest index] since it is only dummy memory
        for (s=0;s<S-2;s++) {
            INfNV=IN(2*s+1,2*t+1,G-1); INfNE=IN(2*s+2,2*t+1,G-1); INfSV=IN(2*s+1,2*t+2,G-1); INfSE=IN(2*s+2,2*t+2,G-1);
            INcNV=IN(s,t,G); INcNE=IN(s+1,t,G); INcSV=IN(s,t+1,G); INcSE=IN(s+1,t+1,G);
            dxVV=dx[INdx(2*s,G-1)];
            dxV=dx[INdx(2*s+1,G-1)];
            dxE=dx[INdx(2*s+2,G-1)];
            dxEE=dx[INdx(2*s+3,G-1)];
            dzNN=dz[INdz(2*t,G-1)];
            dzN=dz[INdz(2*t+1,G-1)];
            dzS=dz[INdz(2*t+2,G-1)];
            dzSS=dz[INdz(2*t+3,G-1)];
            dxtot=0.5*(dxVV+dxV+dxE+dxEE);
            dztot=0.5*(dzNN+dzN+dzS+dzSS);
            //NV
            dimlx=0.5*dxVV/dxtot;
            dimlz=0.5*dzNN/dztot;
            if (s==0 && t==0) dw=0; else dw=(1-dimlx)*(1-dimlz); //weight for NV node on coarse grid i s set to zero if in a corner
            fine[INfNV]+=coarse[INcNV]*dw;
            weight[INfNV]+=dw;
            if (s==S-3 && t==0) dw=0; else dw=dimlx*(1-dimlz); //weight for NE node on coarse grid i s set to zero if in a corner
            fine[INfNV]+=coarse[INcNE]*dw;
            weight[INfNV]+=dw;
            if (s==0 && t==T-3) dw=0; else dw=(1-dimlx)*dimlz; //weight for SV node on coarse grid i s set to zero if in a corner
            fine[INfNV]+=coarse[INcSV]*dw;
            weight[INfNV]+=dw;
            if (s==S-3 && t==T-3) dw=0; else dw=dimlx*dimlz; //weight for SE node on coarse grid i s set to zero if in a corner
            fine[INfNV]+=coarse[INcSE]*dw;
            weight[INfNV]+=dw;
            //NE
            dimlx=0.5*(dxVV+dxV+dxE)/dxtot;
            dimlz=0.5*dzNN/dztot;
            if (s==0 && t==0) dw=0; else dw=(1-dimlx)*(1-dimlz); //weight for NV node on coarse grid i s set to zero if in a corner
            fine[INfNE]+=coarse[INcNV]*dw;
            weight[INfNE]+=dw;
            if (s==S-3 && t==0) dw=0; else dw=dimlx*(1-dimlz); //weight for NE node on coarse grid i s set to zero if in a corner
            fine[INfNE]+=coarse[INcNE]*dw;
            weight[INfNE]+=dw;
            if (s==0 && t==T-3) dw=0; else dw=(1-dimlx)*dimlz; //weight for SV node on coarse grid i s set to zero if in a corner
            fine[INfNE]+=coarse[INcSV]*dw;
            weight[INfNE]+=dw;
            if (s==S-3 && t==T-3) dw=0; else dw=dimlx*dimlz; //weight for SE node on coarse grid i s set to zero if in a corner
            fine[INfNE]+=coarse[INcSE]*dw;
            weight[INfNE]+=dw;
            //SV
            dimlx=0.5*dxVV/dxtot;
            dimlz=0.5*(dzNN+dzN+dzS)/dztot;
            if (s==0 && t==0) dw=0; else dw=(1-dimlx)*(1-dimlz); //weight for NV node on coarse grid i s set to zero if in a corner
            fine[INfSV]+=coarse[INcNV]*dw;
            weight[INfSV]+=dw;
            if (s==S-3 && t==0) dw=0; else dw=dimlx*(1-dimlz); //weight for NE node on coarse grid i s set to zero if in a corner
            fine[INfSV]+=coarse[INcNE]*dw;
            weight[INfSV]+=dw;
            if (s==0 && t==T-3) dw=0; else dw=(1-dimlx)*dimlz; //weight for SV node on coarse grid i s set to zero if in a corner
            fine[INfSV]+=coarse[INcSV]*dw;
            weight[INfSV]+=dw;
            if (s==S-3 && t==T-3) dw=0; else dw=dimlx*dimlz; //weight for SE node on coarse grid i s set to zero if in a corner
            fine[INfSV]+=coarse[INcSE]*dw;
            weight[INfSV]+=dw;
            //SE
            dimlx=0.5*(dxVV+dxV+dxE)/dxtot;
            dimlz=0.5*(dzNN+dzN+dzS)/dztot;
            if (s==0 && t==0) dw=0; else dw=(1-dimlx)*(1-dimlz); //weight for NV node on coarse grid i s set to zero if in a corner
            fine[INfSE]+=coarse[INcNV]*dw;
            weight[INfSE]+=dw;
            if (s==S-3 && t==0) dw=0; else dw=dimlx*(1-dimlz); //weight for NE node on coarse grid i s set to zero if in a corner
            fine[INfSE]+=coarse[INcNE]*dw;
            weight[INfSE]+=dw;
            if (s==0 && t==T-3) dw=0; else dw=(1-dimlx)*dimlz; //weight for SV node on coarse grid i s set to zero if in a corner
            fine[INfSE]+=coarse[INcSV]*dw;
            weight[INfSE]+=dw;
            if (s==S-3 && t==T-3) dw=0; else dw=dimlx*dimlz; //weight for SE node on coarse grid i s set to zero if in a corner
            fine[INfSE]+=coarse[INcSE]*dw;
            weight[INfSE]+=dw;  } }
    
    s=-1;
    for (t=0;t<T-2;t++) {
        INcNE=IN(s+1,t,G); INcSE=IN(s+1,t+1,G);
        INfNE=IN(2*s+2,2*t+1,G-1); INfSE=IN(2*s+2,2*t+2,G-1);
        dzNN=dz[INdz(2*t,G-1)];
        dzN=dz[INdz(2*t+1,G-1)];
        dzS=dz[INdz(2*t+2,G-1)];
        dzSS=dz[INdz(2*t+3,G-1)];
        dztot=0.5*(dzNN+dzN+dzS+dzSS);
        dimlx=0.5; // put to 0.5 since each fine point should have a total weight of 0.5
        //N(E)
        dimlz=0.5*dzNN/dztot;
        if (t==0) dw=0; else dw=dimlx*(1-dimlz); //weight for NE node on coarse grid i s set to zero if in a corner
        fine[INfNE]+=coarse[INcNE]*dw;
        weight[INfNE]+=dw;
        if (t==T-3) dw=0; else dw=dimlx*dimlz; //weight for SE node on coarse grid i s set to zero if in a corner
        fine[INfNE]+=coarse[INcSE]*dw;
        weight[INfNE]+=dw;
        //S(E)
        dimlz=0.5*(dzNN+dzN+dzS)/dztot;
        if (t==0) dw=0; else dw=dimlx*(1-dimlz); //weight for NE node on coarse grid i s set to zero if in a corner
        fine[INfSE]+=coarse[INcNE]*dw;
        weight[INfSE]+=dw;
        if (t==T-3) dw=0; else dw=dimlx*dimlz; //weight for SE node on coarse grid i s set to zero if in a corner
        fine[INfSE]+=coarse[INcSE]*dw;
        weight[INfSE]+=dw;
    }
    s=S-2;
    for (t=0;t<T-2;t++) {
        INcNV=IN(s,t,G); INcSV=IN(s,t+1,G);
        INfNV=IN(2*s+1,2*t+1,G-1); INfSV=IN(2*s+1,2*t+2,G-1);
        dzNN=dz[INdz(2*t,G-1)];
        dzN=dz[INdz(2*t+1,G-1)];
        dzS=dz[INdz(2*t+2,G-1)];
        dzSS=dz[INdz(2*t+3,G-1)];
        dztot=0.5*(dzNN+dzN+dzS+dzSS);
        dimlx=0.5; // put to 0.5 since each fine point should have a total weight of 0.5
        //N(V)
        dimlz=0.5*dzNN/dztot;
        if (t==0) dw=0; else dw=(1-dimlx)*(1-dimlz); //weight for NV node on coarse grid i s set to zero if in a corner
        fine[INfNV]+=coarse[INcNV]*dw;
        weight[INfNV]+=dw;
        if (t==T-3) dw=0; else dw=(1-dimlx)*dimlz; //weight for SV node on coarse grid i s set to zero if in a corner
        fine[INfNV]+=coarse[INcSV]*dw;
        weight[INfNV]+=dw;
        //S(V)
        dimlz=0.5*(dzNN+dzN+dzS)/dztot;
        if (t==0) dw=0; else dw=(1-dimlx)*(1-dimlz); //weight for NV node on coarse grid i s set to zero if in a corner
        fine[INfSV]+=coarse[INcNV]*dw;
        weight[INfSV]+=dw;
        if (t==T-3) dw=0; else dw=(1-dimlx)*dimlz; //weight for SV node on coarse grid i s set to zero if in a corner
        fine[INfSV]+=coarse[INcSV]*dw;
        weight[INfSV]+=dw;
    }
    t=-1;
    for (s=0;s<S-2;s++) {
        INcSV=IN(s,t+1,G); INcSE=IN(s+1,t+1,G);
        INfSV=IN(2*s+1,2*t+2,G-1); INfSE=IN(2*s+2,2*t+2,G-1);
        dxVV=dx[INdx(2*s,G-1)];
        dxV=dx[INdx(2*s+1,G-1)];
        dxE=dx[INdx(2*s+2,G-1)];
        dxEE=dx[INdx(2*s+3,G-1)];
        dxtot=0.5*(dxVV+dxV+dxE+dxEE);
        dimlz=0.5;
        //(S)V
        dimlx=0.5*dxVV/dxtot;
        if (s==0) dw=0; else dw=(1-dimlx)*dimlz; //weight for SV node on coarse grid i s set to zero if in a corner
        fine[INfSV]+=coarse[INcSV]*dw;
        weight[INfSV]+=dw;
        if (s==S-3) dw=0; else dw=dimlx*dimlz; //weight for SE node on coarse grid i s set to zero if in a corner
        fine[INfSV]+=coarse[INcSE]*dw;
        weight[INfSV]+=dw;
        //(S)E
        dimlx=0.5*(dxVV+dxV+dxE)/dxtot;
        if (s==0) dw=0; else dw=(1-dimlx)*dimlz; //weight for SV node on coarse grid i s set to zero if in a corner
        fine[INfSE]+=coarse[INcSV]*dw;
        weight[INfSE]+=dw;
        if (s==S-3) dw=0; else dw=dimlx*dimlz; //weight for SE node on coarse grid i s set to zero if in a corner
        fine[INfSE]+=coarse[INcSE]*dw;
        weight[INfSE]+=dw;
    }
    t=T-2;
    for (s=0;s<S-2;s++) {
        INcNV=IN(s,t,G); INcNE=IN(s+1,t,G);
        INfNV=IN(2*s+1,2*t+1,G-1); INfNE=IN(2*s+2,2*t+1,G-1);
        dxVV=dx[INdx(2*s,G-1)];
        dxV=dx[INdx(2*s+1,G-1)];
        dxE=dx[INdx(2*s+2,G-1)];
        dxEE=dx[INdx(2*s+3,G-1)];
        dxtot=0.5*(dxVV+dxV+dxE+dxEE);
        dimlz=0.5;
        //(N)V
        dimlx=0.5*dxVV/dxtot;
        if (s==0) dw=0; else dw=(1-dimlx)*(1-dimlz); //weight for NV node on coarse grid i s set to zero if in a corner
        fine[INfNV]+=coarse[INcNV]*dw;
        weight[INfNV]+=dw;
        if (s==S-3) dw=0; else dw=dimlx*(1-dimlz); //weight for NE node on coarse grid i s set to zero if in a corner
        fine[INfNV]+=coarse[INcNE]*dw;
        weight[INfNV]+=dw;
        //(N)E
        dimlx=0.5*(dxVV+dxV+dxE)/dxtot;
        if (s==0) dw=0; else dw=(1-dimlx)*(1-dimlz); //weight for NV node on coarse grid i s set to zero if in a corner
        fine[INfNE]+=coarse[INcNV]*dw;
        weight[INfNE]+=dw;
        if (s==S-3) dw=0; else dw=dimlx*(1-dimlz); //weight for NE node on coarse grid i s set to zero if in a corner
        fine[INfNE]+=coarse[INcNE]*dw;
        weight[INfNE]+=dw;
    }      
    //#pragma omp parallel for schedule(static,65) default(shared) private(s,t,st)
    for (t=0;t<Tnew;t++) {
        for (s=0;s<Snew;s++) {
            st=IN(s,t,G-1);
            if (weight[st]>0) {
                fine[st]=fine[st]/weight[st];
            }
            weight[st]=0;
        }
    }
}

void finevxgrid(double *coarse,double *fine,int G) {
    double dzNN,dzN,dzS,dzSS,dxV,dxE,dztot,dxtot,dw,dimlx,dimlz;
    int S,T,Snew,Tnew,s,t,st;
    int INfNV,INfNE,INfSV,INfSE;
    int INcNV,INcNE,INcSV,INcSE;
    S=(Sbase-1)/pow2[G]+1;
    T=(Tbase-1)/pow2[G]+1;
    Snew=(S-1)*2+1;
    Tnew=(T-1)*2+1;
    //#pragma omp parallel for schedule(static,65) default(shared) private(s,t)
    for (t=0;t<Tnew;t++) { //deleting old values
        for (s=0;s<Snew;s++) {
            fine[IN(s,t,G-1)]=0; } }
    
    //#pragma omp parallel for schedule(static,65) default(shared) private(s,t,dzNN,dzN,dzS,dzSS,dxV,dxE,dztot,dxtot,dw,dimlx,dimlz,INfNV,INfNE,INfSV,INfSE,INcNV,INcNE,INcSV,INcSE)
    for (t=0;t<T-2;t++) {
        for (s=0;s<S-1;s++) { //the coarse grid[Snew-1] and grid[Tnew-2] gets assigned values but these will not be used
            INfNV=IN(2*s,2*t+1,G-1); INfNE=IN(2*s+1,2*t+1,G-1); INfSV=IN(2*s,2*t+2,G-1); INfSE=IN(2*s+1,2*t+2,G-1);
            INcNV=IN(s,t,G); INcNE=IN(s+1,t,G); INcSV=IN(s,t+1,G); INcSE=IN(s+1,t+1,G);
            dxV=dx[INdx(2*s,G-1)];
            dxE=dx[INdx(2*s+1,G-1)];
            dzNN=dz[INdz(2*t,G-1)];
            dzN=dz[INdz(2*t+1,G-1)];
            dzS=dz[INdz(2*t+2,G-1)];
            dzSS=dz[INdz(2*t+3,G-1)];
            dxtot=dxV+dxE;
            dztot=0.5*(dzNN+dzN+dzS+dzSS);
            //fine NV - only connections with coarse NV and SV
            dimlx=0;
            dimlz=dimlz=0.5*dzNN/dztot;
            if (s==0) dw=0; else dw=(1-dimlx)*(1-dimlz); //weight for NV node on fine grid is set to zero if s==0 because it is then a boundary point
            fine[INfNV]+=dw*coarse[INcNV];
            weight[INfNV]+=dw;
            if (s==0) dw=0; else dw=(1-dimlx)*dimlz; //weight for SV node on fine grid is set to zero if s==0 because it is then a boundary point
            fine[INfNV]+=dw*coarse[INcSV];
            weight[INfNV]+=dw;
            //fine NE - connections with all 4 coarse nodes
            dimlx=dxV/dxtot;
            dimlz=dimlz=0.5*dzNN/dztot;
            dw=(1-dimlx)*(1-dimlz);
            fine[INfNE]+=dw*coarse[INcNV];
            weight[INfNE]+=dw;
            dw=dimlx*(1-dimlz);
            fine[INfNE]+=dw*coarse[INcNE];
            weight[INfNE]+=dw;
            dw=(1-dimlx)*dimlz;
            fine[INfNE]+=dw*coarse[INcSV];
            weight[INfNE]+=dw;
            dw=dimlx*dimlz;
            fine[INfNE]+=dw*coarse[INcSE];
            weight[INfNE]+=dw;
            
            //fine SV - only connections with coarse NV and SV
            dimlx=0;
            dimlz=0.5*(dzNN+dzN+dzS)/dztot;
            if (s==0) dw=0; else dw=(1-dimlx)*(1-dimlz); //weight for NV node on fine grid is set to zero if s==0 because it is then a boundary point
            fine[INfSV]+=dw*coarse[INcNV];
            weight[INfSV]+=dw;
            if (s==0) dw=0; else dw=(1-dimlx)*dimlz; //weight for SV node on fine grid is set to zero if s==0 because it is then a boundary point
            fine[INfSV]+=dw*coarse[INcSV];
            weight[INfSV]+=dw;
            //fine SE - connections with all 4 coarse nodes
            dimlx=dxV/dxtot;
            dimlz=0.5*(dzNN+dzN+dzS)/dztot;
            dw=(1-dimlx)*(1-dimlz);
            fine[INfSE]+=dw*coarse[INcNV];
            weight[INfSE]+=dw;
            dw=dimlx*(1-dimlz);
            fine[INfSE]+=dw*coarse[INcNE];
            weight[INfSE]+=dw;
            dw=(1-dimlx)*dimlz;
            fine[INfSE]+=dw*coarse[INcSV];
            weight[INfSE]+=dw;
            dw=dimlx*dimlz;
            fine[INfSE]+=dw*coarse[INcSE];
            weight[INfSE]+=dw;
        }
    }
    //#pragma omp parallel for schedule(static,65) default(shared) private(s,t,st)
    for (t=0;t<Tnew;t++) {
        for (s=0;s<Snew;s++) {
            st=IN(s,t,G-1);
            if (weight[st]>0) {
                fine[st]=fine[st]/weight[st];
            }
            weight[st]=0;
        }
    }
}

void finevzgrid(double *coarse,double *fine,int G) {
    double dzN,dzS,dxVV,dxV,dxE,dxEE,dztot,dxtot,dw,dimlx,dimlz;
    int S,T,Snew,Tnew,s,t,st;
    int INfNV,INfNE,INfSV,INfSE;
    int INcNV,INcNE,INcSV,INcSE;
    S=(Sbase-1)/pow2[G]+1;
    T=(Tbase-1)/pow2[G]+1;
    Snew=(S-1)*2+1;
    Tnew=(T-1)*2+1;
    //#pragma omp parallel for schedule(static,65) default(shared) private(s,t)
    for (t=0;t<Tnew;t++) { //deleting old values
        for (s=0;s<Snew;s++) {
            fine[IN(s,t,G-1)]=0; } }
    
    //#pragma omp parallel for schedule(static,65) default(shared) private(s,t,dzN,dzS,dxVV,dxV,dxE,dxEE,dztot,dxtot,dw,dimlx,dimlz,INfNV,INfNE,INfSV,INfSE,INcNV,INcNE,INcSV,INcSE)
    for (t=0;t<T-1;t++) {  //the coarse grid[Tnew-1] and grid[Snew-2] gets assigned values but these will not be used
        for (s=0;s<S-2;s++) {
            INfNV=IN(2*s+1,2*t,G-1); INfNE=IN(2*s+2,2*t,G-1); INfSV=IN(2*s+1,2*t+1,G-1); INfSE=IN(2*s+2,2*t+1,G-1);
            INcNV=IN(s,t,G); INcNE=IN(s+1,t,G); INcSV=IN(s,t+1,G); INcSE=IN(s+1,t+1,G);
            dxVV=dx[INdx(2*s,G-1)];
            dxV=dx[INdx(2*s+1,G-1)];
            dxE=dx[INdx(2*s+2,G-1)];
            dxEE=dx[INdx(2*s+3,G-1)];
            dzN=dz[INdz(2*t,G-1)];
            dzS=dz[INdz(2*t+1,G-1)];
            dxtot=0.5*(dxVV+dxV+dxE+dxEE);
            dztot=dzN+dzS;
            //fine NV - only connections with coarse NV and NE
            dimlx=0.5*dxVV/dxtot;
            dimlz=0;
            if (t==0) dw=0; else dw=(1-dimlx)*(1-dimlz);
            fine[INfNV]+=coarse[INcNV]*dw;
            weight[INfNV]+=dw;
            if (t==0) dw=0; else dw=dimlx*(1-dimlz);
            fine[INfNV]+=coarse[INcNE]*dw;
            weight[INfNV]+=dw;
            //fine NE - only connections with coarse NV and NE
            dimlx=0.5*(dxVV+dxV+dxE)/dxtot;
            dimlz=0;
            if (t==0) dw=0; else dw=(1-dimlx)*(1-dimlz);
            fine[INfNE]+=coarse[INcNV]*dw;
            weight[INfNE]+=dw;
            if (t==0) dw=0; else dw=dimlx*(1-dimlz);
            fine[INfNE]+=coarse[INcNE]*dw;
            weight[INfNE]+=dw;
            
            //fine SV - connections with all 4 coarse nodes
            dimlx=0.5*dxVV/dxtot;
            dimlz=dzN/dztot;
            dw=(1-dimlx)*(1-dimlz);
            fine[INfSV]+=coarse[INcNV]*dw;
            weight[INfSV]+=dw;
            dw=dimlx*(1-dimlz);
            fine[INfSV]+=coarse[INcNE]*dw;
            weight[INfSV]+=dw;
            dw=(1-dimlx)*dimlz;
            fine[INfSV]+=coarse[INcSV]*dw;
            weight[INfSV]+=dw;
            dw=dimlx*dimlz;
            fine[INfSV]+=coarse[INcSE]*dw;
            weight[INfSV]+=dw;
            
            //fine SE - connections with all 4 coarse nodes
            dimlx=0.5*(dxVV+dxV+dxE)/dxtot;
            dimlz=dzN/dztot;
            dw=(1-dimlx)*(1-dimlz);
            fine[INfSE]+=coarse[INcNV]*dw;
            weight[INfSE]+=dw;
            dw=dimlx*(1-dimlz);
            fine[INfSE]+=coarse[INcNE]*dw;
            weight[INfSE]+=dw;
            dw=(1-dimlx)*dimlz;
            fine[INfSE]+=coarse[INcSV]*dw;
            weight[INfSE]+=dw;
            dw=dimlx*dimlz;
            fine[INfSE]+=coarse[INcSE]*dw;
            weight[INfSE]+=dw;
        }
    }
    
    //#pragma omp parallel for schedule(static,65) default(shared) private(s,t,st)
    for (t=0;t<Tnew;t++) {
        for (s=0;s<Snew;s++) {
            st=IN(s,t,G-1);
            if (weight[st]>0) {
                fine[st]=fine[st]/weight[st];
            }
            weight[st]=0;
        }
    }
}

void coarsenP(double *fine,double *coarse,int G) {
    double dzNN,dzN,dzS,dzSS,dxVV,dxV,dxE,dxEE,dztot,dxtot,dw,dimlx,dimlz;
    double fvisc,cviscNV,cviscNE,cviscSV,cviscSE;
    int S,T,Snew,Tnew,s,t,st;
    int INfNV,INfNE,INfSV,INfSE;
    int INcNV,INcNE,INcSV,INcSE;
    S=(Sbase-1)/pow2[G]+1;
    T=(Tbase-1)/pow2[G]+1;
    Snew=(S-1)/2+1;
    Tnew=(T-1)/2+1;
    //#pragma omp parallel for schedule(static,65) default(shared) private(s,t)
    for (t=0;t<Tnew;t++) { //delete old values
        for (s=0;s<Snew;s++) {
            coarse[IN(s,t,G+1)]=0;
        }
    }
    //#pragma omp parallel for schedule(static,65) default(shared) private(s,t,dzNN,dzN,dzS,dzSS,dxVV,dxV,dxE,dxEE,dztot,dxtot,dw,dimlx,dimlz,fvisc,cviscNV,cviscNE,cviscSV,cviscSE,INfNV,INfNE,INfSV,INfSE,INcNV,INcNE,INcSV,INcSE)
    for (t=0;t<Tnew-2;t++) { //Loops skip grid[highest index] since it is only dummy memory
        for (s=0;s<Snew-2;s++) {
            INfNV=IN(2*s+1,2*t+1,G); INfNE=IN(2*s+2,2*t+1,G); INfSV=IN(2*s+1,2*t+2,G); INfSE=IN(2*s+2,2*t+2,G);
            INcNV=IN(s,t,G+1); INcNE=IN(s+1,t,G+1); INcSV=IN(s,t+1,G+1); INcSE=IN(s+1,t+1,G+1);
            dxVV=dx[INdx(2*s,G)];
            dxV=dx[INdx(2*s+1,G)];
            dxE=dx[INdx(2*s+2,G)];
            dxEE=dx[INdx(2*s+3,G)];
            dzNN=dz[INdz(2*t,G)];
            dzN=dz[INdz(2*t+1,G)];
            dzS=dz[INdz(2*t+2,G)];
            dzSS=dz[INdz(2*t+3,G)];
            cviscNV=nZn[INcNV];
            cviscNE=nZn[INcNE];
            cviscSV=nZn[INcSV];
            cviscSE=nZn[INcSE];
            
            dxtot=0.5*(dxVV+dxV+dxE+dxEE);
            dztot=0.5*(dzNN+dzN+dzS+dzSS);
            //NV
            fvisc=nZn[INfNV];
            dimlx=0.5*dxVV/dxtot;
            dimlz=0.5*dzNN/dztot;
            dw=(1-dimlx)*(1-dimlz);
            coarse[INcNV]+=fine[INfNV]*dw*fvisc/cviscNV;
            weight[INcNV]+=dw;
            dw=dimlx*(1-dimlz);
            coarse[INcNE]+=fine[INfNV]*dw*fvisc/cviscNE;
            weight[INcNE]+=dw;
            dw=(1-dimlx)*dimlz;
            coarse[INcSV]+=fine[INfNV]*dw*fvisc/cviscSV;
            weight[INcSV]+=dw;
            dw=dimlx*dimlz;
            coarse[INcSE]+=fine[INfNV]*dw*fvisc/cviscSE;
            weight[INcSE]+=dw;
            //NE
            fvisc=nZn[INfNE];
            dimlx=0.5*(dxVV+dxV+dxE)/dxtot;
            dimlz=0.5*dzNN/dztot;
            dw=(1-dimlx)*(1-dimlz);
            coarse[INcNV]+=fine[INfNE]*dw*fvisc/cviscNV;
            weight[INcNV]+=dw;
            dw=dimlx*(1-dimlz);
            coarse[INcNE]+=fine[INfNE]*dw*fvisc/cviscNE;
            weight[INcNE]+=dw;
            dw=(1-dimlx)*dimlz;
            coarse[INcSV]+=fine[INfNE]*dw*fvisc/cviscSV;
            weight[INcSV]+=dw;
            dw=dimlx*dimlz;
            coarse[INcSE]+=fine[INfNE]*dw*fvisc/cviscSE;
            weight[INcSE]+=dw;
            //SV
            fvisc=nZn[INfSV];
            dimlx=0.5*dxVV/dxtot;
            dimlz=0.5*(dzNN+dzN+dzS)/dztot;
            dw=(1-dimlx)*(1-dimlz);
            coarse[INcNV]+=fine[INfSV]*dw*fvisc/cviscNV;
            weight[INcNV]+=dw;
            dw=dimlx*(1-dimlz);
            coarse[INcNE]+=fine[INfSV]*dw*fvisc/cviscNE;
            weight[INcNE]+=dw;
            dw=(1-dimlx)*dimlz;
            coarse[INcSV]+=fine[INfSV]*dw*fvisc/cviscSV;
            weight[INcSV]+=dw;
            dw=dimlx*dimlz;
            coarse[INcSE]+=fine[INfSV]*dw*fvisc/cviscSE;
            weight[INcSE]+=dw;
            //SE
            fvisc=nZn[INfSE];
            dimlx=0.5*(dxVV+dxV+dxE)/dxtot;
            dimlz=0.5*(dzNN+dzN+dzS)/dztot;
            dw=(1-dimlx)*(1-dimlz);
            coarse[INcNV]+=fine[INfSE]*dw*fvisc/cviscNV;
            weight[INcNV]+=dw;
            dw=dimlx*(1-dimlz);
            coarse[INcNE]+=fine[INfSE]*dw*fvisc/cviscNE;
            weight[INcNE]+=dw;
            dw=(1-dimlx)*dimlz;
            coarse[INcSV]+=fine[INfSE]*dw*fvisc/cviscSV;
            weight[INcSV]+=dw;
            dw=dimlx*dimlz;
            coarse[INcSE]+=fine[INfSE]*dw*fvisc/cviscSE;
            weight[INcSE]+=dw; } }
    
    s=-1;
    for (t=0;t<Tnew-2;t++) {
        INcNE=IN(s+1,t,G+1); INcSE=IN(s+1,t+1,G+1);
        INfNE=IN(2*s+2,2*t+1,G); INfSE=IN(2*s+2,2*t+2,G);
        dzNN=dz[INdz(2*t,G)];
        dzN=dz[INdz(2*t+1,G)];
        dzS=dz[INdz(2*t+2,G)];
        dzSS=dz[INdz(2*t+3,G)];
        cviscNE=nZn[INcNE];
        cviscSE=nZn[INcSE];
        dztot=0.5*(dzNN+dzN+dzS+dzSS);
        dimlx=0.5; // put to 0.5 since each fine point should have a total weight of 0.5
        //N(E)
        fvisc=nZn[INfNE];
        dimlz=0.5*dzNN/dztot;
        dw=dimlx*(1-dimlz);
        coarse[INcNE]+=fine[INfNE]*dw*fvisc/cviscNE;
        weight[INcNE]+=dw;
        dw=dimlx*dimlz;
        coarse[INcSE]+=fine[INfNE]*dw*fvisc/cviscSE;
        weight[INcSE]+=dw;
        //S(E)
        fvisc=nZn[INfSE];
        dimlz=0.5*(dzNN+dzN+dzS)/dztot;
        dw=dimlx*(1-dimlz);
        coarse[INcNE]+=fine[INfSE]*dw*fvisc/cviscNE;
        weight[INcNE]+=dw;
        dw=dimlx*dimlz;
        coarse[INcSE]+=fine[INfSE]*dw*fvisc/cviscSE;
        weight[INcSE]+=dw;
    }
    s=Snew-2;
    for (t=0;t<Tnew-2;t++) {
        INcNV=IN(s,t,G+1); INcSV=IN(s,t+1,G+1);
        INfNV=IN(2*s+1,2*t+1,G); INfSV=IN(2*s+1,2*t+2,G);
        dzNN=dz[INdz(2*t,G)];
        dzN=dz[INdz(2*t+1,G)];
        dzS=dz[INdz(2*t+2,G)];
        dzSS=dz[INdz(2*t+3,G)];
        cviscNV=nZn[INcNV];
        cviscSV=nZn[INcSV];
        dztot=0.5*(dzNN+dzN+dzS+dzSS);
        dimlx=0.5; // put to 0.5 since each fine point should have a total weight of 0.5
        //N(V)
        fvisc=nZn[INfNV];
        dimlz=0.5*dzNN/dztot;
        dw=dimlx*(1-dimlz);
        coarse[INcNV]+=fine[INfNV]*dw*fvisc/cviscNV;
        weight[INcNV]+=dw;
        dw=dimlx*dimlz;
        coarse[INcSV]+=fine[INfNV]*dw*fvisc/cviscSV;
        weight[INcSV]+=dw;
        //S(V)
        fvisc=nZn[INfSV];
        dimlz=0.5*(dzNN+dzN+dzS)/dztot;
        dw=dimlx*(1-dimlz);
        coarse[INcNV]+=fine[INfSV]*dw*fvisc/cviscNV;
        weight[INcNV]+=dw;
        dw=dimlx*dimlz;
        coarse[INcSV]+=fine[INfSV]*dw*fvisc/cviscSV;
        weight[INcSV]+=dw;
    }
    t=-1;
    for (s=0;s<Snew-2;s++) {
        INcSV=IN(s,t+1,G+1); INcSE=IN(s+1,t+1,G+1);
        INfSV=IN(2*s+1,2*t+2,G); INfSE=IN(2*s+2,2*t+2,G);
        dxVV=dx[INdx(2*s,G)];
        dxV=dx[INdx(2*s+1,G)];
        dxE=dx[INdx(2*s+2,G)];
        dxEE=dx[INdx(2*s+3,G)];
        cviscSV=nZn[INcSV];
        cviscSE=nZn[INcSE];
        dxtot=0.5*(dxVV+dxV+dxE+dxEE);
        dimlz=0.5;
        //(S)V
        fvisc=nZn[INfSV];
        dimlx=0.5*dxVV/dxtot;
        dw=(1-dimlx)*dimlz;
        coarse[INcSV]+=fine[INfSV]*dw*fvisc/cviscSV;
        weight[INcSV]+=dw;
        dw=dimlx*dimlz;
        coarse[INcSE]+=fine[INfSV]*dw*fvisc/cviscSE;
        weight[INcSE]+=dw;
        //(S)E
        fvisc=nZn[INfSE];
        dimlx=0.5*(dxVV+dxV+dxE)/dxtot;
        dw=(1-dimlx)*dimlz;
        coarse[INcSV]+=fine[INfSE]*dw*fvisc/cviscSV;
        weight[INcSV]+=dw;
        dw=dimlx*dimlz;
        coarse[INcSE]+=fine[INfSE]*dw*fvisc/cviscSE;
        weight[INcSE]+=dw;
    }
    t=Tnew-2;
    for (s=0;s<Snew-2;s++) {
        INcNV=IN(s,t,G+1); INcNE=IN(s+1,t,G+1);
        INfNV=IN(2*s+1,2*t+1,G); INfNE=IN(2*s+2,2*t+1,G);
        dxVV=dx[INdx(2*s,G)];
        dxV=dx[INdx(2*s+1,G)];
        dxE=dx[INdx(2*s+2,G)];
        dxEE=dx[INdx(2*s+3,G)];
        cviscNV=nZn[INcNV];
        cviscNE=nZn[INcNE];
        dxtot=0.5*(dxVV+dxV+dxE+dxEE);
        dimlz=0.5;
        //(N)V
        fvisc=nZn[INfNV];
        dimlx=0.5*dxVV/dxtot;
        dw=(1-dimlx)*(1-dimlz);
        coarse[INcNV]+=fine[INfNV]*dw*fvisc/cviscNV;
        weight[INcNV]+=dw;
        dw=dimlx*(1-dimlz);
        coarse[INcNE]+=fine[INfNV]*dw*fvisc/cviscNE;
        weight[INcNE]+=dw;
        //(N)E
        fvisc=nZn[INfNE];
        dimlx=0.5*(dxVV+dxV+dxE)/dxtot;
        dw=(1-dimlx)*(1-dimlz);
        coarse[INcNV]+=fine[INfNE]*dw*fvisc/cviscNV;
        weight[INcNV]+=dw;
        dw=dimlx*(1-dimlz);
        coarse[INcNE]+=fine[INfNE]*dw*fvisc/cviscNE;
        weight[INcNE]+=dw;
    }
    //#pragma omp parallel for schedule(static,65) default(shared) private(s,t,st)
    for (t=0;t<Tnew-1;t++) {
        for (s=0;s<Snew-1;s++) {
            st=IN(s,t,G+1);
            if (weight[st]>0) {
                coarse[st]=coarse[st]/weight[st];
            }
            weight[st]=0;
        }
    }
}

void coarsenPgrid(double *fine,double *coarse,int G) {
    double dzNN,dzN,dzS,dzSS,dxVV,dxV,dxE,dxEE,dztot,dxtot,dw,dimlx,dimlz;
    int S,T,Snew,Tnew,s,t,st;
    int INfNV,INfNE,INfSV,INfSE;
    int INcNV,INcNE,INcSV,INcSE;
    S=(Sbase-1)/pow2[G]+1;
    T=(Tbase-1)/pow2[G]+1;
    Snew=(S-1)/2+1;
    Tnew=(T-1)/2+1;
    //#pragma omp parallel for schedule(static,65) default(shared) private(s,t)
    for (t=0;t<Tnew;t++) { //delete old values
        for (s=0;s<Snew;s++) {
            coarse[IN(s,t,G+1)]=0;
        }
    }
    //#pragma omp parallel for schedule(static,65) default(shared) private(s,t,dzNN,dzN,dzS,dzSS,dxVV,dxV,dxE,dxEE,dztot,dxtot,dw,dimlx,dimlz,INfNV,INfNE,INfSV,INfSE,INcNV,INcNE,INcSV,INcSE)
    for (t=0;t<Tnew-2;t++) { //Loops skip grid[highest index] since it is only dummy memory
        for (s=0;s<Snew-2;s++) {
            INfNV=IN(2*s+1,2*t+1,G); INfNE=IN(2*s+2,2*t+1,G); INfSV=IN(2*s+1,2*t+2,G); INfSE=IN(2*s+2,2*t+2,G);
            INcNV=IN(s,t,G+1); INcNE=IN(s+1,t,G+1); INcSV=IN(s,t+1,G+1); INcSE=IN(s+1,t+1,G+1);
            dxVV=dx[INdx(2*s,G)];
            dxV=dx[INdx(2*s+1,G)];
            dxE=dx[INdx(2*s+2,G)];
            dxEE=dx[INdx(2*s+3,G)];
            dzNN=dz[INdz(2*t,G)];
            dzN=dz[INdz(2*t+1,G)];
            dzS=dz[INdz(2*t+2,G)];
            dzSS=dz[INdz(2*t+3,G)];
            dxtot=0.5*(dxVV+dxV+dxE+dxEE);
            dztot=0.5*(dzNN+dzN+dzS+dzSS);
            //NV
            dimlx=0.5*dxVV/dxtot;
            dimlz=0.5*dzNN/dztot;
            dw=(1-dimlx)*(1-dimlz);
            coarse[INcNV]+=fine[INfNV]*dw;
            weight[INcNV]+=dw;
            dw=dimlx*(1-dimlz);
            coarse[INcNE]+=fine[INfNV]*dw;
            weight[INcNE]+=dw;
            dw=(1-dimlx)*dimlz;
            coarse[INcSV]+=fine[INfNV]*dw;
            weight[INcSV]+=dw;
            dw=dimlx*dimlz;
            coarse[INcSE]+=fine[INfNV]*dw;
            weight[INcSE]+=dw;
            //NE
            dimlx=0.5*(dxVV+dxV+dxE)/dxtot;
            dimlz=0.5*dzNN/dztot;
            dw=(1-dimlx)*(1-dimlz);
            coarse[INcNV]+=fine[INfNE]*dw;
            weight[INcNV]+=dw;
            dw=dimlx*(1-dimlz);
            coarse[INcNE]+=fine[INfNE]*dw;
            weight[INcNE]+=dw;
            dw=(1-dimlx)*dimlz;
            coarse[INcSV]+=fine[INfNE]*dw;
            weight[INcSV]+=dw;
            dw=dimlx*dimlz;
            coarse[INcSE]+=fine[INfNE]*dw;
            weight[INcSE]+=dw;
            //SV
            dimlx=0.5*dxVV/dxtot;
            dimlz=0.5*(dzNN+dzN+dzS)/dztot;
            dw=(1-dimlx)*(1-dimlz);
            coarse[INcNV]+=fine[INfSV]*dw;
            weight[INcNV]+=dw;
            dw=dimlx*(1-dimlz);
            coarse[INcNE]+=fine[INfSV]*dw;
            weight[INcNE]+=dw;
            dw=(1-dimlx)*dimlz;
            coarse[INcSV]+=fine[INfSV]*dw;
            weight[INcSV]+=dw;
            dw=dimlx*dimlz;
            coarse[INcSE]+=fine[INfSV]*dw;
            weight[INcSE]+=dw;
            //SE
            dimlx=0.5*(dxVV+dxV+dxE)/dxtot;
            dimlz=0.5*(dzNN+dzN+dzS)/dztot;
            dw=(1-dimlx)*(1-dimlz);
            coarse[INcNV]+=fine[INfSE]*dw;
            weight[INcNV]+=dw;
            dw=dimlx*(1-dimlz);
            coarse[INcNE]+=fine[INfSE]*dw;
            weight[INcNE]+=dw;
            dw=(1-dimlx)*dimlz;
            coarse[INcSV]+=fine[INfSE]*dw;
            weight[INcSV]+=dw;
            dw=dimlx*dimlz;
            coarse[INcSE]+=fine[INfSE]*dw;
            weight[INcSE]+=dw; } }
    
    s=-1;
    for (t=0;t<Tnew-2;t++) {
        INcNE=IN(s+1,t,G+1); INcSE=IN(s+1,t+1,G+1);
        INfNE=IN(2*s+2,2*t+1,G); INfSE=IN(2*s+2,2*t+2,G);
        dzNN=dz[INdz(2*t,G)];
        dzN=dz[INdz(2*t+1,G)];
        dzS=dz[INdz(2*t+2,G)];
        dzSS=dz[INdz(2*t+3,G)];
        dztot=0.5*(dzNN+dzN+dzS+dzSS);
        dimlx=0.5; // put to 0.5 since each fine point should have a total weight of 0.5
        //N(E)
        dimlz=0.5*dzNN/dztot;
        dw=dimlx*(1-dimlz);
        coarse[INcNE]+=fine[INfNE]*dw;
        weight[INcNE]+=dw;
        dw=dimlx*dimlz;
        coarse[INcSE]+=fine[INfNE]*dw;
        weight[INcSE]+=dw;
        //S(E)
        dimlz=0.5*(dzNN+dzN+dzS)/dztot;
        dw=dimlx*(1-dimlz);
        coarse[INcNE]+=fine[INfSE]*dw;
        weight[INcNE]+=dw;
        dw=dimlx*dimlz;
        coarse[INcSE]+=fine[INfSE]*dw;
        weight[INcSE]+=dw;
    }
    s=Snew-2;
    for (t=0;t<Tnew-2;t++) {
        INcNV=IN(s,t,G+1); INcSV=IN(s,t+1,G+1);
        INfNV=IN(2*s+1,2*t+1,G); INfSV=IN(2*s+1,2*t+2,G);
        dzNN=dz[INdz(2*t,G)];
        dzN=dz[INdz(2*t+1,G)];
        dzS=dz[INdz(2*t+2,G)];
        dzSS=dz[INdz(2*t+3,G)];
        dztot=0.5*(dzNN+dzN+dzS+dzSS);
        dimlx=0.5; // put to 0.5 since each fine point should have a total weight of 0.5
        //N(V)
        dimlz=0.5*dzNN/dztot;
        dw=dimlx*(1-dimlz);
        coarse[INcNV]+=fine[INfNV]*dw;
        weight[INcNV]+=dw;
        dw=dimlx*dimlz;
        coarse[INcSV]+=fine[INfNV]*dw;
        weight[INcSV]+=dw;
        //S(V)
        dimlz=0.5*(dzNN+dzN+dzS)/dztot;
        dw=dimlx*(1-dimlz);
        coarse[INcNV]+=fine[INfSV]*dw;
        weight[INcNV]+=dw;
        dw=dimlx*dimlz;
        coarse[INcSV]+=fine[INfSV]*dw;
        weight[INcSV]+=dw;
    }
    t=-1;
    for (s=0;s<Snew-2;s++) {
        INcSV=IN(s,t+1,G+1); INcSE=IN(s+1,t+1,G+1);
        INfSV=IN(2*s+1,2*t+2,G); INfSE=IN(2*s+2,2*t+2,G);
        dxVV=dx[INdx(2*s,G)];
        dxV=dx[INdx(2*s+1,G)];
        dxE=dx[INdx(2*s+2,G)];
        dxEE=dx[INdx(2*s+3,G)];
        dxtot=0.5*(dxVV+dxV+dxE+dxEE);
        dimlz=0.5;
        //(S)V
        dimlx=0.5*dxVV/dxtot;
        dw=(1-dimlx)*dimlz;
        coarse[INcSV]+=fine[INfSV]*dw;
        weight[INcSV]+=dw;
        dw=dimlx*dimlz;
        coarse[INcSE]+=fine[INfSV]*dw;
        weight[INcSE]+=dw;
        //(S)E
        dimlx=0.5*(dxVV+dxV+dxE)/dxtot;
        dw=(1-dimlx)*dimlz;
        coarse[INcSV]+=fine[INfSE]*dw;
        weight[INcSV]+=dw;
        dw=dimlx*dimlz;
        coarse[INcSE]+=fine[INfSE]*dw;
        weight[INcSE]+=dw;
    }
    t=Tnew-2;
    for (s=0;s<Snew-2;s++) {
        INcNV=IN(s,t,G+1); INcNE=IN(s+1,t,G+1);
        INfNV=IN(2*s+1,2*t+1,G); INfNE=IN(2*s+2,2*t+1,G);
        dxVV=dx[INdx(2*s,G)];
        dxV=dx[INdx(2*s+1,G)];
        dxE=dx[INdx(2*s+2,G)];
        dxEE=dx[INdx(2*s+3,G)];
        dxtot=0.5*(dxVV+dxV+dxE+dxEE);
        dimlz=0.5;
        //(N)V
        dimlx=0.5*dxVV/dxtot;
        dw=(1-dimlx)*(1-dimlz);
        coarse[INcNV]+=fine[INfNV]*dw;
        weight[INcNV]+=dw;
        dw=dimlx*(1-dimlz);
        coarse[INcNE]+=fine[INfNV]*dw;
        weight[INcNE]+=dw;
        //(N)E
        dimlx=0.5*(dxVV+dxV+dxE)/dxtot;
        dw=(1-dimlx)*(1-dimlz);
        coarse[INcNV]+=fine[INfNE]*dw;
        weight[INcNV]+=dw;
        dw=dimlx*(1-dimlz);
        coarse[INcNE]+=fine[INfNE]*dw;
        weight[INcNE]+=dw;
    }
    //#pragma omp parallel for schedule(static,65) default(shared) private(s,t,st)
    for (t=0;t<Tnew-1;t++) {
        for (s=0;s<Snew-1;s++) {
            st=IN(s,t,G+1);
            if (weight[st]>0) {
                coarse[st]=coarse[st]/weight[st];
            }
            weight[st]=0;
        }
    }
}


void coarsenvxgrid(double *fine,double *coarse,int G) {
    double dzNN,dzN,dzS,dzSS,dxV,dxE,dztot,dxtot,dw,dimlx,dimlz;
    int S,T,Snew,Tnew,s,t,st;
    int INfNV,INfNE,INfSV,INfSE;
    int INcNV,INcNE,INcSV,INcSE;
    S=(Sbase-1)/pow2[G]+1;
    T=(Tbase-1)/pow2[G]+1;
    Snew=(S-1)/2+1;
    Tnew=(T-1)/2+1;
    //#pragma omp parallel for schedule(static,65) default(shared) private(s,t)
    for (t=0;t<Tnew;t++) { //delete old values
        for (s=0;s<Snew;s++) {
            coarse[IN(s,t,G+1)]=0;
        } }
    //#pragma omp parallel for schedule(static,65) default(shared) private(s,t,dzNN,dzN,dzS,dzSS,dxV,dxE,dztot,dxtot,dw,dimlx,dimlz,INfNV,INfNE,INfSV,INfSE,INcNV,INcNE,INcSV,INcSE)
    for (t=0;t<Tnew-2;t++) {
        for (s=0;s<Snew-1;s++) { //the coarse[Snew-1] and coarse[Tnew-2] gets assigned values but these will not be used
            INfNV=IN(2*s,2*t+1,G); INfNE=IN(2*s+1,2*t+1,G); INfSV=IN(2*s,2*t+2,G); INfSE=IN(2*s+1,2*t+2,G);
            INcNV=IN(s,t,G+1); INcNE=IN(s+1,t,G+1); INcSV=IN(s,t+1,G+1); INcSE=IN(s+1,t+1,G+1);
            dxV=dx[INdx(2*s,G)];
            dxE=dx[INdx(2*s+1,G)];
            dzNN=dz[INdz(2*t,G)];
            dzN=dz[INdz(2*t+1,G)];
            dzS=dz[INdz(2*t+2,G)];
            dzSS=dz[INdz(2*t+3,G)];
            dxtot=dxV+dxE;
            dztot=0.5*(dzNN+dzN+dzS+dzSS);
            
            //fine NV - only connections with coarse NV and SV
            dimlx=0;
            dimlz=dimlz=0.5*dzNN/dztot;
            if (s==0) dw=0; else dw=(1-dimlx)*(1-dimlz); //weight for NV node on fine grid is set to zero if s==0 because it is then a boundary point
            coarse[INcNV]+=dw*fine[INfNV];
            weight[INcNV]+=dw;
            if (s==0) dw=0; else dw=(1-dimlx)*dimlz; //weight for SV node on fine grid is set to zero if s==0 because it is then a boundary point
            coarse[INcSV]+=dw*fine[INfNV];
            weight[INcSV]+=dw;
            //fine NE - connections with all 4 coarse nodes
            dimlx=dxV/dxtot;
            dimlz=dimlz=0.5*dzNN/dztot;
            dw=(1-dimlx)*(1-dimlz);
            coarse[INcNV]+=dw*fine[INfNE];
            weight[INcNV]+=dw;
            dw=dimlx*(1-dimlz);
            coarse[INcNE]+=dw*fine[INfNE];
            weight[INcNE]+=dw;
            dw=(1-dimlx)*dimlz;
            coarse[INcSV]+=dw*fine[INfNE];
            weight[INcSV]+=dw;
            dw=dimlx*dimlz;
            coarse[INcSE]+=dw*fine[INfNE];
            weight[INcSE]+=dw;
            //fine SV - only connections with coarse NV and SV
            dimlx=0;
            dimlz=0.5*(dzNN+dzN+dzS)/dztot;
            if (s==0) dw=0; else dw=(1-dimlx)*(1-dimlz); //weight for NV node on fine grid is set to zero if s==0 because it is then a boundary point
            coarse[INcNV]+=dw*fine[INfSV];
            weight[INcNV]+=dw;
            if (s==0) dw=0; else dw=(1-dimlx)*dimlz; //weight for SV node on fine grid is set to zero if s==0 because it is then a boundary point
            coarse[INcSV]+=dw*fine[INfSV];
            weight[INcSV]+=dw;
            //fine SE - connections with all 4 coarse nodes
            dimlx=dxV/dxtot;
            dimlz=0.5*(dzNN+dzN+dzS)/dztot;
            dw=(1-dimlx)*(1-dimlz);
            coarse[INcNV]+=dw*fine[INfSE];
            weight[INcNV]+=dw;
            dw=dimlx*(1-dimlz);
            coarse[INcNE]+=dw*fine[INfSE];
            weight[INcNE]+=dw;
            dw=(1-dimlx)*dimlz;
            coarse[INcSV]+=dw*fine[INfSE];
            weight[INcSV]+=dw;
            dw=dimlx*dimlz;
            coarse[INcSE]+=dw*fine[INfSE];
            weight[INcSE]+=dw;
            
        }
    }
    //#pragma omp parallel for schedule(static,65) default(shared) private(s,t,st)
    for (t=0;t<Tnew;t++) {
        for (s=0;s<Snew;s++) {
            st=IN(s,t,G+1);
            if (weight[st]>0) {
                coarse[st]=coarse[st]/weight[st];
            }
            weight[st]=0;
        }
    }
}

void coarsenvzgrid(double *fine,double *coarse,int G) {
    double dzN,dzS,dxVV,dxV,dxE,dxEE,dztot,dxtot,dw,dimlx,dimlz;
    int S,T,Snew,Tnew,s,t,st;
    int INfNV,INfNE,INfSV,INfSE;
    int INcNV,INcNE,INcSV,INcSE;
    S=(Sbase-1)/pow2[G]+1;
    T=(Tbase-1)/pow2[G]+1;
    Snew=(S-1)/2+1;
    Tnew=(T-1)/2+1;
    //#pragma omp parallel for schedule(static,65) default(shared) private(s,t)
    for (t=0;t<Tnew;t++) { //delete old values
        for (s=0;s<Snew;s++) {
            coarse[IN(s,t,G+1)]=0;
        } }
    //#pragma omp parallel for schedule(static,65) default(shared) private(s,t,dzN,dzS,dxVV,dxV,dxE,dxEE,dztot,dxtot,dw,dimlx,dimlz,INfNV,INfNE,INfSV,INfSE,INcNV,INcNE,INcSV,INcSE)
    for (t=0;t<Tnew-1;t++) {  //the coarse grid[Tnew-1] and grid[Snew-2] gets assigned values but these will not be used
        for (s=0;s<Snew-2;s++) {
            INfNV=IN(2*s+1,2*t,G); INfNE=IN(2*s+2,2*t,G); INfSV=IN(2*s+1,2*t+1,G); INfSE=IN(2*s+2,2*t+1,G);
            INcNV=IN(s,t,G+1); INcNE=IN(s+1,t,G+1); INcSV=IN(s,t+1,G+1); INcSE=IN(s+1,t+1,G+1);
            dxVV=dx[INdx(2*s,G)];
            dxV=dx[INdx(2*s+1,G)];
            dxE=dx[INdx(2*s+2,G)];
            dxEE=dx[INdx(2*s+3,G)];
            dzN=dz[INdz(2*t,G)];
            dzS=dz[INdz(2*t+1,G)];
            dxtot=0.5*(dxVV+dxV+dxE+dxEE);
            dztot=dzN+dzS;
            //fine NV - only connections with coarse NV and NE
            dimlx=0.5*dxVV/dxtot;
            dimlz=0;
            if (t==0) dw=0; else dw=(1-dimlx)*(1-dimlz);
            coarse[INcNV]+=fine[INfNV]*dw;
            weight[INcNV]+=dw;
            if (t==0) dw=0; else dw=dimlx*(1-dimlz);
            coarse[INcNE]+=fine[INfNV]*dw;
            weight[INcNE]+=dw;
            
            //fine NE - only connections with coarse NV and NE
            dimlx=0.5*(dxVV+dxV+dxE)/dxtot;
            dimlz=0;
            if (t==0) dw=0; else dw=(1-dimlx)*(1-dimlz);
            coarse[INcNV]+=fine[INfNE]*dw;
            weight[INcNV]+=dw;
            if (t==0) dw=0; else dw=dimlx*(1-dimlz);
            coarse[INcNE]+=fine[INfNE]*dw;
            weight[INcNE]+=dw;
            
            //fine SV - connections with all 4 coarse nodes
            dimlx=0.5*dxVV/dxtot;
            dimlz=dzN/dztot;
            dw=(1-dimlx)*(1-dimlz);
            coarse[INcNV]+=fine[INfSV]*dw;
            weight[INcNV]+=dw;
            dw=dimlx*(1-dimlz);
            coarse[INcNE]+=fine[INfSV]*dw;
            weight[INcNE]+=dw;
            dw=(1-dimlx)*dimlz;
            coarse[INcSV]+=fine[INfSV]*dw;
            weight[INcSV]+=dw;
            dw=dimlx*dimlz;
            coarse[INcSE]+=fine[INfSV]*dw;
            weight[INcSE]+=dw;
            
            //fine SE - connections with all 4 coarse nodes
            dimlx=0.5*(dxVV+dxV+dxE)/dxtot;
            dimlz=dzN/dztot;
            dw=(1-dimlx)*(1-dimlz);
            coarse[INcNV]+=fine[INfSE]*dw;
            weight[INcNV]+=dw;
            dw=dimlx*(1-dimlz);
            coarse[INcNE]+=fine[INfSE]*dw;
            weight[INcNE]+=dw;
            dw=(1-dimlx)*dimlz;
            coarse[INcSV]+=fine[INfSE]*dw;
            weight[INcSV]+=dw;
            dw=dimlx*dimlz;
            coarse[INcSE]+=fine[INfSE]*dw;
            weight[INcSE]+=dw;
        }
    }
    //#pragma omp parallel for schedule(static,65) default(shared) private(s,t,st)
    for (t=0;t<Tnew;t++) {
        for (s=0;s<Snew;s++) {
            st=IN(s,t,G+1);
            if (weight[st]>0) {
                coarse[st]=coarse[st]/weight[st];
            }
            weight[st]=0;
        }
    }
}

void resetgrid(double *grid,int G) {
    int S,T,s,t;
    S=(Sbase-1)/pow2[G]+1;
    T=(Tbase-1)/pow2[G]+1;
    //#pragma omp parallel for schedule(static,65) default(shared) private(s,t)
    for (t=0;t<T;t++) {
        for (s=0;s<S;s++) {
            grid[IN(s,t,G)]=0; } }
}

void restoright(int G) { // let the coarsened residuals become the right hand side of the equations
    int S,T,s,t,st;
    S=(Sbase-1)/pow2[G]+1;
    T=(Tbase-1)/pow2[G]+1;
    //#pragma omp parallel for schedule(static,65) default(shared) private(s,t,st)
    for (t=0;t<T;t++) {
        for (s=0;s<S;s++) {
            st=IN(s,t,G);
            if (s==S-1 || t==T-1 || s==0 && (t==0 || t==T-2) || s==S-2 && (t==0 || t==T-2))	{
                rightP[st]=0; }
            else {
                rightP[st]=resP[st]; }
            if (s==0 || s==S-1 || t==0 || t==T-2 || t==T-1) {
                rightx[st]=0; }
            else {
                rightx[st]=resx[st]; }
            if (t==0 || t==T-1 || s==0 || s==S-2 || s==S-1) {
                rightz[st]=0; }
            else {
                rightz[st]=resz[st]; }
        }
    }
}

void pointupdate(int s,int t,int G,int S,int T) {
    int st,s_1t,sp1t,stp1,st_1;
    st=IN(s,t,G);
    s_1t=IN(s-1,t,G);
    sp1t=IN(s+1,t,G);
    stp1=IN(s,t+1,G);
    st_1=IN(s,t-1,G);
    // IN(s-1,t+1,G) IN(s+1,t-1,G) used only once
    if (s==S-1 || t==T-1 || s==0 && (t==0 || t==T-2) || s==S-2 && (t==0 || t==T-2))	{
        if (s==S-2 && (t==0 || t==T-2)) { //gradient-rand i nord�st og syd�st
            P[st]=P[s_1t]; }
        if (s==0 && (t==0 || t==T-2)) { //gradient-rand i nordvest og sydvest
            P[st]=P[sp1t];
        }
    }
    else { //P=Pold-div(v)*local_visc or if rightP>0 P=Pold+(rightP-div(v))*local_visc=Pold+*local_visc
        resP[st]=rightP[st]-((vx[sp1t]-vx[st])/dx[INdx(s,G)]+(vz[stp1]-vz[st])/dz[INdz(t,G)]);
        P[st]+=resP[st]*nZn[st]*alphaP;
    }
    //randbetingelser for vx
    if (s==0 || s==S-1 || t==0 || t==T-2 || t==T-1) {
        if (s==0) {
            if (leftbound==0) vx[st]=0;
            if (leftbound==1) vx[st]=0;
            if (G==0) vx[st]=leftboundvel;
        }
        if (s==S-1) {
            if (rightbound==0) vx[st]=0;
            if (rightbound==1) vx[st]=0;
            if (G==0) vx[st]=rightboundvel;
        }
        if (t==0) {
            if (upbound==0) vx[st]=0;
            if (upbound==1) vx[st]=vx[stp1];
        }
        if (t==T-2) {
            if (downbound==0) vx[st]=0;
            if (downbound==1) vx[st]=vx[st_1];
        } }
    else {
        resx[st]=rightx[st]-(cx0[st]*P[s_1t]+cx1[st]*vx[s_1t]+cx2[st]*vz[s_1t]+cx3[st]*vz[IN(s-1,t+1,G)]+cx4[st]*vx[st_1]+cx5[st]*P[st]+cx6[st]*vx[st]+cx7[st]*vz[st]+cx8[st]*vx[stp1]+cx9[st]*vz[stp1]+cx10[st]*vx[sp1t]);
        vx[st]+=alphamom*resx[st]/cx6[st];
    }
    //randbetingelser for vz
    if (t==0 || t==T-1 || s==0 || s==S-2 || s==S-1) {
        if (s==0) {
            if (leftbound==0) vz[st]=0;
            if (leftbound==1) vz[st]=vz[sp1t];
        }
        if (s==S-2) {
            if (rightbound==0) vz[st]=0;
            if (rightbound==1) vz[st]=vz[s_1t];
        }
        if (t==0) {
            if (upbound==0) vz[st]=0;
            if (upbound==1) vz[st]=0;
            if (G==0) vz[st]=upboundvel;
        }
        if (t==T-1) {
            if (downbound==0) vz[st]=0;
            if (downbound==1) vz[st]=0;
            if (G==0) vz[st]=downboundvel;
        } }
    else {
        resz[st]=rightz[st]-(cz0[st]*vz[s_1t]+cz1[st]*P[st_1]+cz2[st]*vx[st_1]+cz3[st]*vz[st_1]+cz4[st]*P[st]+cz5[st]*vx[st]+cz6[st]*vz[st]+cz7[st]*vz[stp1]+cz8[st]*vx[IN(s+1,t-1,G)]+cz9[st]*vx[sp1t]+cz10[st]*vz[sp1t]);
        vz[st]+=alphamom*resz[st]/cz6[st];
    }
}

void iterate(int times,int up,int basevisclevel,int G) { //times: number of times to iterate. If up is true, the residuals corresponding to the solution will be calculated in the end
        double Perr;
        int count,S,T,Nb,b,s,t,st;
        S=(Sbase-1)/pow2[G]+1;
        T=(Tbase-1)/pow2[G]+1;

        //if ((S-1)/Nb>8) {}
        //if (S>=129) {}
        //else Nb=1;
        //omp_set_num_threads(3);
        //if (S==513) omp_set_num_threads(2);
        //if (S==257) omp_set_num_threads(3);
        //if (S==129) omp_set_num_threads(1);
        //if (S==65) omp_set_num_threads(1);
        //if (S<65) omp_set_num_threads(1);
/*
        int Nthreads=omp_get_max_threads();
        int start_t[Nthreads+1]; //last start_t is T, number of Ts
        start_t[0]=0;
        for (e=0;e<Nthreads;e++) start_t[e+1]=start_t[e]+T/Nthreads;
        start_t[Nthreads]=T;
        
        //#pragma omp parallel for private(s,t,st,count)
        for (count=0;count<times;count++) {
            int ID= omp_get_thread_num();
            //printf("ID: %d start: %d end: %d Nthreads: %d\n",ID,start_t[ID],start_t[ID+1],Nthreads);
            for (t=start_t[ID];t<start_t[ID+1];t++)     {
                for (s=0;s<S;s++)    pointupdate(s, t, G, S, T);
                //for (t=0;t<T;t++)    pointupdate(s,t,G,S,T);
                //for (t=0;t<T;t++)    pointupdate(s,t,G,S,T);
                //for (t=0;t<T;t++)    pointupdate(s,t,G,S,T);
            }
        }*/
        
   /*     int minT=2048;
        if (T<minT) {
            for (count=0;count<times;count++) {
                //#pragma omp parallel for schedule(static) default(shared) private(s, t, st)
                for (t=0;t<T;t++)     {
                    for (s=0;s<S;s++)    pointupdate(s, t, G, S, T);
                } } }
        else {
            for (count=0;count<times;count++) {
                #pragma omp parallel for schedule(static) default(shared) private(s, t, st)
                for (t=0;t<T;t=t+2)     {
                    for (s=0;s<S;s++)    pointupdate(s, t, G, S, T);
                }
                #pragma omp parallel for schedule(static) default(shared) private(s, t, st)
                for (t=1;t<T;t=t+2)     {
                    for (s=0;s<S;s++)    pointupdate(s, t, G, S, T);
                }
            }
        }*/
        for (count=0;count<times;count++) {
            #pragma omp parallel for schedule(static) default(shared) private(s, t)
            for (t=0;t<T;t++)     {
                for (s=0;s<S;s++)    pointupdate(s, t, G, S, T);
            }
        }
    /*    for (count=0;count<times;count++) {
            #pragma omp parallel for schedule(static) default(shared) private(s, t, st)
            for (st=firstelementgrid[G];st<firstelementgrid[G+1];st++)     {
                t=(st-firstelementgrid[G])/S; s=st-firstelementgrid[G]-t*S;
                for (s=0;s<S;s++)    pointupdate(s, t, G, S, T);
            }
        } */
        
        //IDEA: Parallelize this loop!
        
  /*      if (G==0)
            for (count=0;count<times;count++) {
            #pragma omp parallel for schedule(static) default(shared) private(s, t, st)
            for (t=0;t<T;t++)     {
                for (s=0;s<S;s++)    pointupdate(s, t, G, S, T);
            }
            }
        
        else
            for (count=0;count<times;count++) {
            //#pragma omp parallel for schedule(static) default(shared) private(s, t, st)
            for (t=0;t<T;t++)     {
                for (s=0;s<S;s++)    pointupdate(s, t, G, S, T);
            }
            }*/

        Perr=P[IN(1,0,G)]; // fejl i P, der opst�r i hvert iterativt trin. Alle P korrigeres herfor. V�rdien i P(1,0) s�ttes til at skulle v�re 0 (randbetingelse for P)
        //#pragma omp parallel for schedule(static,65) default(shared) private(s,t,st)
        for (t=0;t<T-1;t++)	{
            for (s=0;s<S;s++){
                st=IN(s,t,G);
                P[st]-=Perr; }  } //principielt rammes hj�rnepunkter, men det korrigeres der for senere med en gradientrand
        //omp_set_num_threads(1);
}


double maxerr_mom_press(int G) {
    int S,T,s,t,st;
    double length,localerror,maxerror2;
    double error2;
    maxerror2=0;
    S=(Sbase-1)/pow2[G]+1;
    T=(Tbase-1)/pow2[G]+1;
    for (t=1;t<T-2;t++)	{ //resx
        for (s=1;s<S-1;s++) {
            st=IN(s,t,G);
            length=(0.5*(dx[INdx(s,G)]+dx[INdx(s-1,G)])+dz[INdz(t,G)])/2;
            localerror=resx[st]*length;
            error2=localerror*localerror;
            if (error2>maxerror2) maxerror2=error2;
            if (isnan(error2)) maxerror2=1e100;
        } }
    for (t=1;t<T-1;t++)	{
        for (s=1;s<S-2;s++) { //resz
            st=IN(s,t,G);
            length=(0.5*(dz[INdz(t,G)]+dz[INdz(t-1,G)])+dx[INdx(s,G)])/2;
            localerror=resz[st]*length;
            error2=localerror*localerror;
            if (error2>maxerror2) maxerror2=error2;
            if (isnan(error2)) maxerror2=1e100;
        } }
    return sqrt(maxerror2);
}

double maxerr_cont_press(int G) {
    int S,T,s,t,st;
    double localerror,error2;
    double maxerror2=0;
    S=(Sbase-1)/pow2[G]+1;
    T=(Tbase-1)/pow2[G]+1;
    for (t=0;t<T-1;t++)	{
        for (s=0;s<S-1;s++) {
            if (s==0 && (t==0 || t==T-2) || s==S-2 && (t==0 || t==T-2))	{ }
            else {
                st=IN(s,t,G);
                localerror=resP[st]*nZn[st];
                error2=localerror*localerror;
                if (error2>maxerror2) maxerror2=error2;
                if (isnan(error2)) maxerror2=1e100;
            } } }
    return sqrt(maxerror2);
}

void Vcycle(int G,int baseiterations) {
    iterate(baseiterations,0,1,G);
    if (G==multigridlevels) {
        iterate(2*baseiterations,0,1,G);
    }
    else {
        iterate(1,1,1,G);
        coarsenP(resP,resP,G);
        coarsenvxgrid(resx,resx,G);
        coarsenvzgrid(resz,resz,G);
        restoright(G+1);
        resetgrid(P,G+1);
        resetgrid(vx,G+1);
        resetgrid(vz,G+1);
        Vcycle(G+1,baseiterations*2);
        finePgrid(P,resP,G+1); //res is used as dummy memory to rember the interpolated increments from coarse to a finer level
        finevxgrid(vx,resx,G+1);
        finevzgrid(vz,resz,G+1);
        correct(G);
    }
    iterate(baseiterations,0,1,G);
}

void Vcycle1(int G,int baseiterations) {
    iterate(baseiterations,0,0,G);
    if (G==multigridlevels) {
        iterate(2*baseiterations,0,0,G);
    }
    else {
        iterate(1,1,0,G);
        coarsenP(resP,resP,G);
        coarsenvxgrid(resx,resx,G);
        coarsenvzgrid(resz,resz,G);
        restoright(G+1);
        resetgrid(P,G+1);
        resetgrid(vx,G+1);
        resetgrid(vz,G+1);
        Vcycle(G+1,baseiterations*2);
        finePgrid(P,resP,G+1); //res is used as dummy memory to rember the interpolated increments from coarse to a finer level
        finevxgrid(vx,resx,G+1);
        finevzgrid(vz,resz,G+1);
        correct(G);
    }
    iterate(baseiterations,0,0,G);
}

void Fcycle(int G,int baseiterations) {
    iterate(baseiterations,0,1,G);
    if (G==multigridlevels) {
        iterate(2*baseiterations,0,1,G);
    }
    else {
        iterate(1,1,1,G);
        coarsenP(resP,resP,G);
        coarsenvxgrid(resx,resx,G);
        coarsenvzgrid(resz,resz,G);
        restoright(G+1);
        resetgrid(P,G+1);
        resetgrid(vx,G+1);
        resetgrid(vz,G+1);
        Fcycle(G+1,baseiterations*2);
        Vcycle(G+1,baseiterations*2);
        finePgrid(P,resP,G+1); //res is used as dummy memory to rember the interpolated increments from coarse to a finer level
        finevxgrid(vx,resx,G+1);
        finevzgrid(vz,resz,G+1);
        correct(G);
    }
    iterate(baseiterations,0,1,G);
}

void Fcycle1(int G,int baseiterations) {
    iterate(baseiterations,0,0,G);
    if (G==multigridlevels) {
        iterate(2*baseiterations,0,0,G);
    }
    else {
        iterate(1,1,0,G);
        coarsenP(resP,resP,G);
        coarsenvxgrid(resx,resx,G);
        coarsenvzgrid(resz,resz,G);
        restoright(G+1);
        resetgrid(P,G+1);
        resetgrid(vx,G+1);
        resetgrid(vz,G+1);
        Fcycle(G+1,baseiterations*2);
        Vcycle(G+1,baseiterations*2);
        finePgrid(P,resP,G+1); //res is used as dummy memory to rember the interpolated increments from coarse to a finer level
        finevxgrid(vx,resx,G+1);
        finevzgrid(vz,resz,G+1);
        correct(G);
    }
    iterate(baseiterations,0,0,G);
}

int multigridVcycle(double tolerance,int baseiterations,int maxiterations,int *isfinished) {
    int count=0;
    *isfinished=1;
    //fprintf(statusfile,"Vcycles: ");
    while (maxerr_mom_press(0)>tolerance || maxerr_cont_press(0)>tolerance) {
        fprintf(statusfile,"%d ",count); 
        Vcycle(0,baseiterations);
        count++;
        if (count>maxiterations) { *isfinished=0; break; }
    }
    fprintf(statusfile,"\n");
    return count;
}

int multigridFcycle(double tolerance,int baseiterations,int maxiterations,int *isfinished) {
    int count=0;
    *isfinished=1;
    fprintf(statusfile,"Fcycles: ");
    while (maxerr_mom_press(0)>tolerance || maxerr_cont_press(0)>tolerance) {
        fprintf(statusfile,"%d ",count); 
        Fcycle(0,baseiterations);
        count++;
        if (count>maxiterations) { *isfinished=0; break; }
    }
    fprintf(statusfile,"\n");
    return count;
}


void updatevisc(double alpha) {
    G=0;
    S=(Sbase-1)/pow2[G]+1;
    T=(Tbase-1)/pow2[G]+1;
    if (alpha==1) {
        for (t=0;t<T;t++)	{
            for (s=0;s<S;s++){
                nZn[IN(s,t,G)]=inZn[IN(s,t,G)];
                nZs[IN(s,t,G)]=inZs[IN(s,t,G)]; } } }
    else {
        for (t=0;t<T;t++)	{
            for (s=0;s<S;s++){
                nZn[IN(s,t,G)]=minvisc*pow(inZn[IN(s,t,G)]/minvisc,alpha);
                nZs[IN(s,t,G)]=minvisc*pow(inZs[IN(s,t,G)]/minvisc,alpha); } } }
}

double absmax(double a,double b) {
    if (fabs(a)>fabs(b)) return fabs(a);
    else return fabs(b);
}
double absmin(double a,double b) {
    if (fabs(a)<fabs(b)) return fabs(a);
    else return fabs(b);
}

void calcextremevisc() {
    int s,t,st;
    maxvisc=0;
    for (t=0;t<Tbase;t++)	{
        for (s=0;s<Sbase;s++){       
            st=IN(s,t,0);
            //maxvisc=absmax(maxvisc,inZn[st]);
            maxvisc=absmax(maxvisc,inZs[st]);
        }
    }
    minvisc=maxvisc;
    for (t=0;t<T;t++)	{
        for (s=0;s<Sbase;s++){
            st=IN(s,t,0);
            //minvisc=absmin(minvisc,inZn[st]);
            minvisc=absmin(minvisc,inZs[st]);
        }
    }
}

void makecoeff() {
    int s,t,G,s_1t,stp1,st_1,sp1t,idxs,idxs_1,idxsp1,idxt,idxt_1,idxtp1,st;
    S=Sbase;
    T=Tbase;
    for (G=0;G<multigridlevels;G++){
        //coarsenvisc(nZs,nZn,nZs,nZn,G);
        coarsenPgrid(nZn,nZn,G);
        coarsenrhogrid(nZs,nZs,G);
        //coarsenrhogrid(rho,rho,G);
        S=(S-1)/2+1;
        T=(T-1)/2+1;
    }
    
    S=Sbase;
    T=Tbase;
    for (G=0;G<=multigridlevels;G++){
        for (t=0;t<T;t++)	{
            for (s=0;s<S;s++){      
                st=IN(s,t,G); s_1t=IN(s-1,t,G); stp1=IN(s,t+1,G); st_1=IN(s,t-1,G); sp1t=IN(s+1,t,G); idxs=INdx(s,G); idxs_1=INdx(s-1,G); idxsp1=INdx(s+1,G); idxt=INdz(t,G); idxt_1=INdz(t-1,G); idxtp1=INdz(t+1,G);
                if (s==0 || s==S-1 || t==0 || t==T-2 || t==T-1) { }
                else {
                    cx0[st]=2/(dx[idxs_1]+dx[idxs]); //P(s-1,t)
                    cx1[st]=4*nZn[s_1t]/(dx[idxs_1]*(dx[idxs_1]+dx[idxs])); //vx(s-1,t)
                    cx2[st]=2*nZs[st]/(dz[idxt]*(dx[idxs_1]+dx[idxs])); //vz(s-1,t)
                    cx3[st]=-2*nZs[stp1]/(dz[idxt]*(dx[idxs_1]+dx[idxs])); //vz(s-1,t+1)
                    cx4[st]=2*nZs[st]/(dz[idxt]*(dz[idxt_1]+dz[idxt]));  //vx(s,t-1)
                    cx5[st]=-2/(dx[idxs_1]+dx[idxs]); //P(s,t)
                    cx6[st]=-(4*nZn[st]/(dx[idxs]*(dx[idxs_1]+dx[idxs]))+4*nZn[s_1t]/(dx[idxs_1]*(dx[idxs_1]+dx[idxs]))+2*nZs[stp1]/(dz[idxt]*(dz[idxt]+dz[idxtp1]))+2*nZs[st]/(dz[idxt]*(dz[idxt_1]+dz[idxt])));   //vx(s,t)
                    cx7[st]=-2*nZs[st]/(dz[idxt]*(dx[idxs_1]+dx[idxs])); //vz(s,t)
                    cx8[st]=2*nZs[stp1]/(dz[idxt]*(dz[idxt]+dz[idxtp1])); //vx(s,t+1)
                    cx9[st]=2*nZs[stp1]/(dz[idxt]*(dx[idxs_1]+dx[idxs])); //vz(s,t+1)
                    cx10[st]=4*nZn[st]/(dx[idxs]*(dx[idxs_1]+dx[idxs])); //vx(s+1,t)
                }
                
                if (t==0 || t==T-1 || s==0 || s==S-2 || s==S-1) { }
                else {
                    cz0[st]=2*nZs[st]/(dx[idxs]*(dx[idxs_1]+dx[idxs]));  //vz(s-1,t)
                    cz1[st]=2/(dz[idxt_1]+dz[idxt]); //P(s,t-1)
                    cz2[st]=2*nZs[st]/(dx[idxs]*(dz[idxt_1]+dz[idxt]));  //vx(s,t-1)
                    cz3[st]=4*nZn[st_1]/(dz[idxt_1]*(dz[idxt_1]+dz[idxt]));  //vz(s,t-1)
                    cz4[st]=-2/(dz[idxt_1]+dz[idxt]);//P(s,t)
                    cz5[st]=-2*nZs[st]/(dx[idxs]*(dz[idxt_1]+dz[idxt]));    //vx(s,t)
                    cz6[st]=-(4*nZn[st]/(dz[idxt]*(dz[idxt_1]+dz[idxt]))+4*nZn[st_1]/(dz[idxt_1]*(dz[idxt_1]+dz[idxt]))+2*nZs[sp1t]/(dx[idxs]*(dx[idxs]+dx[idxsp1]))+2*nZs[st]/(dx[idxs]*(dx[idxs_1]+dx[idxs])));//vz(s,t)
                    cz7[st]=4*nZn[st]/(dz[idxt]*(dz[idxt_1]+dz[idxt]));  //vz(s,t+1)
                    cz8[st]=-2*nZs[sp1t]/(dx[idxs]*(dz[idxt_1]+dz[idxt]));//vx(s+1,t-1)
                    cz9[st]=2*nZs[sp1t]/(dx[idxs]*(dz[idxt_1]+dz[idxt]));  //vx(s+1,t)
                    cz10[st]=2*nZs[sp1t]/(dx[idxs]*(dx[idxs]+dx[idxsp1]));  //vz(s+1,t)
                } }}
        S=(S-1)/2+1;
        T=(T-1)/2+1; }
}

int solve(double tolerance) {
    int mgcycle;
    double currentmaxvisc,alpha,range;
    double momerror,conterror,conterror0,momerror0;
    double *vxsave,*vzsave,*Psave,*rightxsave,*rightzsave,*rightPsave;
    int multimulticycles;
    int Vcycles,convergence;
    double olderror;
    
    alphacorrect=0.5; //alphamom=0.5; alphaP=0.1;
 
    Vcycles=0; multimulticycles=0;
    
    updatevisc(1); makecoeff(); calcextremevisc();
    range=log10(maxvisc/minvisc);
    iterate(2,1,1,0); //perform a single iteration to calculate the residual
    momerror0=maxerr_mom_press(0); conterror0=maxerr_cont_press(0);
    olderror=momerror0;
    if (isinf(momerror0) || isnan(momerror0) || isinf(conterror0) || isnan(conterror0) || conterror0>1e49 || momerror0>1e49 ) {
        fprintf(statusfile,"Resetting grids because non finite error\n");
        resetgrid(vx,0); resetgrid(vz,0); resetgrid(P,0); iterate(2,1,1,0); momerror0=maxerr_mom_press(0); conterror0=maxerr_cont_press(0); olderror=momerror0;
    }
    
    fprintf(statusfile,"Initial error (momentum, continuity): %3.2g %3.2g\n",momerror0,conterror0);
    if (momerror0>tolerance || conterror0>tolerance) { //iterations are needed
        if (range<0.1) { //regular Vcycles - determined by emperical limit of visc contrast
            fprintf(statusfile,"Viscosity range is %g. Performing regular Vcycles\n",range);
            while(momerror0>tolerance || conterror0>tolerance) {
                if (Vcycles%10==0)  {fprintf(statusfile,"Cycle %d error %3.2g.   ",Vcycles,momerror0);  }
                Vcycles++; Vcycle(0,4);
                if (Vcycles>1000) { fprintf(statusfile,"\nError: convergence was not obtained with V-cycles Trying with F-cycles\n");  break; }
                if (Vcycles%10==0) {momerror0=maxerr_mom_press(0); conterror0=maxerr_cont_press(0);}
                if(momerror0>1e3*olderror) {break; }
            } }

        if (momerror0>tolerance || conterror0>tolerance) {
            fprintf(statusfile,"Viscosity range is %g. Performing regular Fcycles\n",range);
            while(momerror0>tolerance || conterror0>tolerance) {
                if (Vcycles%10==0)  {fprintf(statusfile,"Cycle %d error %3.2g.   ",Vcycles,momerror0);  }
                Vcycles++; Fcycle(0,5);
                if (Vcycles>2000) { fprintf(statusfile,"\nError: convergence was not obtained with F-cycles\n"); break; }
                if (Vcycles%10==0) {momerror0=maxerr_mom_press(0); conterror0=maxerr_cont_press(0);}
                if(momerror0>1e4*olderror) {break; }
            } } }
    
    if (momerror0>tolerance || conterror0>tolerance) {
        fprintf(statusfile,"\nError: convergence was not obtained with normal cycles. Tring with secure cycles\n");
        
//         printgrid(rightx, "rightx", iterationnumber);
//         printgrid(rightz, "rightz", iterationnumber);
//         printgrid(gridoldstressxx, "gridoldstressxx", iterationnumber);
//         printgrid(gridoldstressxz, "gridoldstressxz", iterationnumber);
//         printgrid(nZn, "nZn", iterationnumber);
//         printgrid(nZs, "nZs", iterationnumber);
//         printgrid(gridviscn, "gridviscn", iterationnumber);
//         printgrid(gridviscs, "gridviscs", iterationnumber);
//         printgrid(gridrho, "gridrho", iterationnumber);
//         printgrid(resx, "resx", iterationnumber);
//         printgrid(resz, "resz", iterationnumber);
//         printgrid(resP, "resP", iterationnumber);
//         
        resetgrid(vx,0); resetgrid(vz,0); resetgrid(P,0); iterate(2,1,1,0); momerror0=maxerr_mom_press(0); conterror0=maxerr_cont_press(0); olderror=momerror0;
        alphacorrect=0.3;
        while(momerror0>tolerance || conterror0>tolerance) {
            if (Vcycles%10==0)  {fprintf(statusfile,"Cycle %d error %3.2g.   ",Vcycles,momerror0);  }
            Vcycles++; Fcycle(0,20);
            if (Vcycles>4000) { fprintf(statusfile,"\nWARNING: convergence was not obtained!\n"); break; }
            if (Vcycles%10==0) {momerror0=maxerr_mom_press(0); conterror0=maxerr_cont_press(0);}
            if (momerror0>1e6*olderror) {
                fprintf(statusfile,"\nError: Divergence - reducing relaxation coefficients \n");
                resetgrid(vx,0); resetgrid(vz,0); resetgrid(P,0); iterate(2,1,1,0); momerror0=maxerr_mom_press(0); conterror0=maxerr_cont_press(0); olderror=momerror0;
                alphacorrect*=0.5; 
            }
        } }
    if (momerror0<=tolerance && conterror0<=tolerance) { fprintf(statusfile,"\nConvergence was  obtained, # multigrid cycles: %d\n",Vcycles); convergence=1;}
    else { fprintf(statusfile,"\nConvergence was NOT  obtained. solve returns 0 - see call to solve in timestep\n"); convergence=0; }
    return convergence;
    //printf("Multigrid cycles: %d\n",Vcycles);
}
