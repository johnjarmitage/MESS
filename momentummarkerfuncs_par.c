#include <stdio.h>
#include <math.h>
#include <omp.h>
#include <stdlib.h>
#include "main.h"
#include "momentummarkerfuncs_par.h"

void printgrid(double *grid,char *name, int number) {
    FILE *file;
    char a[100];
    sprintf(a,"outfiles/%s%d.txt",name,number);
    file = fopen(a, "w");

    for (t=0;t<Tbase;t++) {
        for (s=0;s<Sbase-1;s++) {
            fprintf(file, "%g ", grid[IN(s, t, 0)]);
        }
        fprintf(file, "%g\n",grid[IN(s, t, 0)]);
    }
    fclose(file);
}


double vxtopos(double xpos, double zpos) {
    int INNV,INNE,INSV,INSE;
    double dlx,dlz;
    int isnorth;
    int s,t;
    s=get_s(xpos); t=get_t(zpos);
    isnorth=zpos<z[t]+dz[t]/2;
    dlx=(xpos-x[s])/dx[s];
    if (isnorth) {
        INSV=IN0(s,t); INSE=IN0(s+1,t);
        if (t==0) { INNV=INSV; INNE=INSE; dlz=0.5; }
        else {
            INNV=IN0(s,t-1); INNE=IN0(s+1,t-1); dlz=2*(zpos-z[t]+0.5*dz[t-1])/(dz[t]+dz[t-1]); } }
    else {
        INNV=IN0(s,t); INNE=IN0(s+1,t);
        if (t==Tbase-2) { INSV=INNV; INSE=INNE; dlz=0.5; }
        else {
            INSV=IN0(s,t+1); INSE=IN0(s+1,t+1); dlz=2*(zpos-z[t]-0.5*dz[t])/(dz[t]+dz[t+1]); }}
    return vx[INNV]*(1-dlx)*(1-dlz)+vx[INNE]*dlx*(1-dlz)+vx[INSV]*(1-dlx)*dlz+vx[INSE]*dlx*dlz;
}

double vztopos(double xpos, double zpos) {
    int INNV,INNE,INSV,INSE;
    double dlx,dlz;
    int iswest;
    int s,t;
    s=get_s(xpos); t=get_t(zpos);
    iswest=xpos<x[s]+dx[s]/2;
    dlz=(zpos-z[t])/dz[t];
    if (iswest) {
        INNE=IN0(s,t); INSE=IN0(s,t+1);
        if (s==0) { INNV=INNE; INSV=INSE; dlx=0.5; }
        else {
            INNV=IN0(s-1,t); INSV=IN0(s-1,t+1); dlx=2*(xpos-x[s]+0.5*dx[s-1])/(dx[s-1]+dx[s]); } }
    else {
        INNV=IN0(s,t); INSV=IN0(s,t+1);
        if (s==Sbase-2) {INNE=INNV; INSE=INSV; dlx=0.5; }
        else {
            INNE=IN0(s+1,t); INSE=IN0(s+1,t+1); dlx=2*(xpos-x[s]-0.5*dx[s])/(dx[s]+dx[s+1]); }}
        return vz[INNV]*(1-dlx)*(1-dlz)+vz[INNE]*dlx*(1-dlz)+vz[INSV]*(1-dlx)*dlz+vz[INSE]*dlx*dlz;
        //return -z[t]*rightboundvel/x[Sbase-1];
}

double vxtomarker(int e) {
    int INNV,INNE,INSV,INSE;
    double dlx,dlz;
    int isnorth;
    int s,t;
    s=mcs[e]; t=mct[e];
    isnorth=mz[e]<z[t]+dz[t]/2;
    dlx=(mx[e]-x[s])/dx[s];
    if (isnorth) {
        INSV=IN0(s,t); INSE=IN0(s+1,t);
        if (t==0) { INNV=INSV; INNE=INSE; dlz=0.5; }
        else {
            INNV=IN0(s,t-1); INNE=IN0(s+1,t-1); dlz=2*(mz[e]-z[t]+0.5*dz[t-1])/(dz[t]+dz[t-1]); } }
    else {
        INNV=IN0(s,t); INNE=IN0(s+1,t);
        if (t==Tbase-2) { INSV=INNV; INSE=INNE; dlz=0.5; }
        else {
            INSV=IN0(s,t+1); INSE=IN0(s+1,t+1); dlz=2*(mz[e]-z[t]-0.5*dz[t])/(dz[t]+dz[t+1]); }}
    return vx[INNV]*(1-dlx)*(1-dlz)+vx[INNE]*dlx*(1-dlz)+vx[INSV]*(1-dlx)*dlz+vx[INSE]*dlx*dlz;
}

double vztomarker(int e) {
    int INNV,INNE,INSV,INSE;
    double dlx,dlz;
    int iswest;
    int s,t;
    s=mcs[e]; t=mct[e];
    iswest=mx[e]<x[s]+dx[s]/2;
    dlz=(mz[e]-z[t])/dz[t];
    if (iswest) {
        INNE=IN0(s,t); INSE=IN0(s,t+1);
        if (s==0) { INNV=INNE; INSV=INSE; dlx=0.5; }
        else {
            INNV=IN0(s-1,t); INSV=IN0(s-1,t+1); dlx=2*(mx[e]-x[s]+0.5*dx[s-1])/(dx[s-1]+dx[s]); } }
    else {
        INNV=IN0(s,t); INSV=IN0(s,t+1);
        if (s==Sbase-2) {INNE=INNV; INSE=INSV; dlx=0.5; }
        else {
            INNE=IN0(s+1,t); INSE=IN0(s+1,t+1); dlx=2*(mx[e]-x[s]-0.5*dx[s])/(dx[s]+dx[s+1]); }}
    return vz[INNV]*(1-dlx)*(1-dlz)+vz[INNE]*dlx*(1-dlz)+vz[INSV]*(1-dlx)*dlz+vz[INSE]*dlx*dlz;
}

void markertorhogrid(int e,double value,double *grid,double *weight) {
    //gives a weighted contribution from a value on a marker to the four sorrounding rho grid points. Division with weights must be performed later
    int s,t;
    int INNV,INNE,INSV,INSE;
    double dlx,dlz,dw;
    s=mcs[e]; t=mct[e];
    INNV=IN0(s,t); INNE=IN0(s+1,t); INSV=IN0(s,t+1); INSE=IN0(s+1,t+1);
    dlx=(mx[e]-x[s])/dx[s];
    dlz=(mz[e]-z[t])/dz[t];
    // NV
    dw=(1-dlx)*(1-dlz); grid[INNV]+=value*dw; weight[INNV]+=dw;
    // NE
    dw=dlx*(1-dlz); grid[INNE]+=value*dw; weight[INNE]+=dw;
    // SV
    dw=(1-dlx)*dlz; grid[INSV]+=value*dw; weight[INSV]+=dw;
    // SE
    dw=dlx*dlz; grid[INSE]+=value*dw; weight[INSE]+=dw;
}

double rhogridtomarker(int e,double *grid) {
    int s,t;
    int INNV,INNE,INSV,INSE;
    double dlx,dlz;
    s=mcs[e]; t=mct[e];
    INNV=IN0(s,t); INNE=IN0(s+1,t); INSV=IN0(s,t+1); INSE=IN0(s+1,t+1);
    dlx=(mx[e]-x[s])/dx[s];
    dlz=(mz[e]-z[t])/dz[t];
    return grid[INNV]*(1-dlx)*(1-dlz)+grid[INNE]*dlx*(1-dlz)+grid[INSV]*(1-dlx)*dlz+grid[INSE]*dlx*dlz;
}

double Pgridtomarker(int e,double *grid) {
    int iswest,isnorth;
    int s,t;
    int INNV,INNE,INSV,INSE;
    double dlx,dlz;
    s=mcs[e]; t=mct[e];
    iswest=mx[e]<x[s]+dx[s]/2;
    isnorth=mz[e]<z[t]+dz[t]/2;
    INNV=0; INNE=0; INSV=0; INSE=0;
    if (iswest) {
        if (isnorth) { //NV
            INSE=IN0(s,t);
            if (s>0) {
                dlx=2*(mx[e]-x[s]+0.5*dx[s-1])/(dx[s-1]+dx[s]);  INSV=IN0(s-1,t);
                if (t>0) { dlz=2*(mz[e]-z[t]+0.5*dz[t-1])/(dz[t-1]+dz[t]); INNE=IN0(s,t-1); INNV=IN0(s-1,t-1); }
                else { dlz=1; } } //weights for northern nodes become 0
            else {
                dlx=1; //weights fo western nodes become 0
                if (t>0) { dlz=2*(mz[e]-z[t]+0.5*dz[t-1])/(dz[t-1]+dz[t]); INNE=IN0(s,t-1); }
                else { dlz=1; } } } //weights for northern nodes become 0
        else { //SV
            INNE=IN0(s,t);
            if (s>0) {
                dlx=2*(mx[e]-x[s]+0.5*dx[s-1])/(dx[s-1]+dx[s]); INNV=IN0(s-1,t);
                if (t<Tbase-2) { dlz=2*(mz[e]-z[t]-0.5*dz[t])/(dz[t]+dz[t+1]); INSE=IN0(s,t+1); INSV=IN0(s-1,t+1); }
                else { dlz=0; } }  //weights for southern nodes become 0
            else {
                dlx=1; //weights fo western nodes become 0
                if (t<Tbase-2) { dlz=2*(mz[e]-z[t]-0.5*dz[t])/(dz[t]+dz[t+1]); INSE=IN0(s,t+1); }
                else { dlz=0; } } } } //weights for southern nodes become 0
    else {
        if (isnorth) { //NE
            INSV=IN0(s,t);
            if (s<Sbase-2) {
                dlx=2*(mx[e]-x[s]-0.5*dx[s])/(dx[s]+dx[s+1]); INSE=IN0(s+1,t);
                if (t>0) { dlz=2*(mz[e]-z[t]+0.5*dz[t-1])/(dz[t-1]+dz[t]); INNV=IN0(s,t-1); INNE=IN0(s+1,t-1); }
                else { dlz=1; } } //weights for northern nodes become 0
            else { 
                dlx=0; //weights fo eastern nodes become 0
                if (t>0) { dlz=2*(mz[e]-z[t]+0.5*dz[t-1])/(dz[t-1]+dz[t]); INNV=IN0(s,t-1); }
                else { dlz=1; } } } //weights for northern nodes become 0
        else { //SE
            INNV=IN0(s,t);
            if (s<Sbase-2) {
                dlx=2*(mx[e]-x[s]-0.5*dx[s])/(dx[s]+dx[s+1]); INNE=IN0(s+1,t);
                if (t<Tbase-2) { dlz=2*(mz[e]-z[t]-0.5*dz[t])/(dz[t]+dz[t+1]); INSV=IN0(s,t+1); INSE=IN0(s+1,t+1); }
                else { dlz=0; } } //weights for southern nodes become 0
            else {
                dlx=0; //weights fo eastern nodes become 0
                if (t<Tbase-2) { dlz=2*(mz[e]-z[t]-0.5*dz[t])/(dz[t]+dz[t+1]); INSV=IN0(s,t+1); }
                else { dlz=0; } }  } } //weights for southern nodes become 0
    return grid[INNV]*(1-dlx)*(1-dlz)+grid[INNE]*dlx*(1-dlz)+grid[INSV]*(1-dlx)*dlz+grid[INSE]*dlx*dlz;
}

void markertoPgrid1(int e,double value,double *grid,double *weight) {
    int iswest,isnorth;
    int s,t;
    int INNV,INNE,INSV,INSE;
    double dlx,dlz,dw;
    s=mcs[e]; t=mct[e];
    iswest=mx[e]<x[s]+dx[s]/2;
    isnorth=mz[e]<z[t]+dz[t]/2;
    INNV=0; INNE=0; INSV=0; INSE=0;
    if (iswest) {
        if (isnorth) { //NV
            INSE=IN0(s,t);
            if (s>0) {
                dlx=2*(mx[e]-x[s]+0.5*dx[s-1])/(dx[s-1]+dx[s]);  INSV=IN0(s-1,t);
                if (t>0) { dlz=2*(mz[e]-z[t]+0.5*dz[t-1])/(dz[t-1]+dz[t]); INNE=IN0(s,t-1); INNV=IN0(s-1,t-1); }
                else { dlz=1; } } //weights for northern nodes become 0
            else {
                dlx=1; //weights fo western nodes become 0
                if (t>0) { dlz=2*(mz[e]-z[t]+0.5*dz[t-1])/(dz[t-1]+dz[t]); INNE=IN0(s,t-1); }
                else { dlz=1; } } } //weights for northern nodes become 0
        else { //SV
            INNE=IN0(s,t);
            if (s>0) {
                dlx=2*(mx[e]-x[s]+0.5*dx[s-1])/(dx[s-1]+dx[s]); INNV=IN0(s-1,t);
                if (t<Tbase-2) { dlz=2*(mz[e]-z[t]-0.5*dz[t])/(dz[t]+dz[t+1]); INSE=IN0(s,t+1); INSV=IN0(s-1,t+1); }
                else { dlz=0; } }  //weights for southern nodes become 0
            else {
                dlx=1; //weights fo western nodes become 0
                if (t<Tbase-2) { dlz=2*(mz[e]-z[t]-0.5*dz[t])/(dz[t]+dz[t+1]); INSE=IN0(s,t+1); }
                else { dlz=0; } } } } //weights for southern nodes become 0
    else {
        if (isnorth) { //NE
            INSV=IN0(s,t);
            if (s<Sbase-2) {
                dlx=2*(mx[e]-x[s]-0.5*dx[s])/(dx[s]+dx[s+1]); INSE=IN0(s+1,t);
                if (t>0) { dlz=2*(mz[e]-z[t]+0.5*dz[t-1])/(dz[t-1]+dz[t]); INNV=IN0(s,t-1); INNE=IN0(s+1,t-1); }
                else { dlz=1; } } //weights for northern nodes become 0
            else { 
                dlx=0; //weights fo eastern nodes become 0
                if (t>0) { dlz=2*(mz[e]-z[t]+0.5*dz[t-1])/(dz[t-1]+dz[t]); INNV=IN0(s,t-1); }
                else { dlz=1; } } } //weights for northern nodes become 0
        else { //SE
            INNV=IN0(s,t);
            if (s<Sbase-2) {
                dlx=2*(mx[e]-x[s]-0.5*dx[s])/(dx[s]+dx[s+1]); INNE=IN0(s+1,t);
                if (t<Tbase-2) { dlz=2*(mz[e]-z[t]-0.5*dz[t])/(dz[t]+dz[t+1]); INSV=IN0(s,t+1); INSE=IN0(s+1,t+1); }
                else { dlz=0; } } //weights for southern nodes become 0
            else {
                dlx=0; //weights fo eastern nodes become 0
                if (t<Tbase-2) { dlz=2*(mz[e]-z[t]-0.5*dz[t])/(dz[t]+dz[t+1]); INSV=IN0(s,t+1); }
                else { dlz=0; } }  } } //weights for southern nodes become 0
    // NV
    dw=(1-dlx)*(1-dlz); grid[INNV]+=value*dw; weight[INNV]+=dw;
    // NE
    dw=dlx*(1-dlz); grid[INNE]+=value*dw; weight[INNE]+=dw;
    // SV
    dw=(1-dlx)*dlz; grid[INSV]+=value*dw; weight[INSV]+=dw;
    // SE
    dw=dlx*dlz; grid[INSE]+=value*dw; weight[INSE]+=dw;

}

void markervisctogrid(int e,double value,double *gridrho,double *gridP,double *weightrho,double *weightP) {
    //gives a weighted contribution from a value on a marker to the nearest visc_s point and to the visc_n in the cell. Division with weights must be performed later
    //for visc_s the weight is as in markertogrid while for visc_n it is 1 at the P node 0.25 at rhode nodes and 0.5 at velocity nodes
    int iswest,isnorth;
    int s,t;
    int INn,INs;
    double dlx,dlz,dw;
    s=mcs[e]; t=mct[e];
    iswest=mx[e]<x[s]+dx[s]/2;
    isnorth=mz[e]<z[t]+dz[t]/2;
    dlx=(mx[e]-x[s])/dx[s];
    dlz=(mz[e]-z[t])/dz[t];
    if (iswest) {
        if (isnorth) { INs=IN0(s,t); dw=(1-dlx)*(1-dlz); }
        else INs=IN0(s,t+1); dw=(1-dlx)*dlz; }
    else {
        if (isnorth) { INs=IN0(s+1,t); dw=dlx*(1-dlz); }
        else INs=IN0(s+1,t+1); dw=dlx*dlz; }
    gridrho[INs]+=value*dw; weightrho[INs]+=dw;
    INn=IN0(s,t); dw=(1-fabs(dlx-0.5))*(1-fabs(dlz-0.5)); // this gives the proper weight for for the visc_n point
    gridP[INn]+=value*dw; weightP[INn]+=dw;
}

void markertoPgrid(int e,double value,double *gridP,double *weightP) {
    //is like markervisctogrid, but only contributions to Pnodes
    int s,t;
    int INn;
    double dlx,dlz,dw;
    s=mcs[e]; t=mct[e];
    dlx=(mx[e]-x[s])/dx[s];
    dlz=(mz[e]-z[t])/dz[t];
    INn=IN0(s,t); dw=(1-fabs(dlx-0.5))*(1-fabs(dlz-0.5)); // this gives the proper weight for for the visc_n point
    gridP[INn]+=value*dw; weightP[INn]+=dw;
}

void markervisctorhogrid(int e,double value,double *gridrho,double *weightrho) {
    // is like markervisctogrid, but only contributions to rhonodes. Not like marker to rhogrid since the cell is divided into 4 areas of contribution
    int iswest,isnorth;
    int s,t;
    int INs;
    double dlx,dlz,dw;
    s=mcs[e]; t=mct[e];
    iswest=mx[e]<x[s]+dx[s]/2;
    isnorth=mz[e]<z[t]+dz[t]/2;
    dlx=(mx[e]-x[s])/dx[s];
    dlz=(mz[e]-z[t])/dz[t];
    if (iswest) {
        if (isnorth) { INs=IN0(s,t); dw=(1-dlx)*(1-dlz); }
        else INs=IN0(s,t+1); dw=(1-dlx)*dlz; }
    else {
        if (isnorth) { INs=IN0(s+1,t); dw=dlx*(1-dlz); }
        else INs=IN0(s+1,t+1); dw=dlx*dlz; }
    gridrho[INs]+=value*dw; weightrho[INs]+=dw;
}

double secdevinv2(double tensorxx,double tensorxz) {
    return tensorxx*tensorxx+tensorxz*tensorxz;
}

double maxvelocity() {
    double maxv=0;
    for (t=0;t<Tbase;t++) { for (s=0;s<Sbase;s++){ st=IN0(s,t); if (maxv<fabs(vx[st])) {maxv=fabs(vx[st]); } if (maxv<fabs(vz[st])) {maxv=fabs(vz[st]); } } }
    return maxv;
}

int get_s(double xpos) {
    int low,high,mid;
    low=0; high=Sbase-1;
    while (high-low>1) { mid=(high+low)/2; if (x[mid]<=xpos) low=mid; else high=mid; }
    return low;
}
int get_t(double zpos) {
    int low,high,mid;
    low=0; high=Tbase-1;
    while (high-low>1) { mid=(high+low)/2; if (z[mid]<=zpos) low=mid; else high=mid; }
    return low;
}

void movemarkers() {
    double optimaldt,dtmove,mxold,mzold;
    int s,t,movecount,numberofmoves,warning;
    warning=0;
    
    optimaldt=mincellsize/maxvelocity();
    numberofmoves=ceil(dt/optimaldt);
    dtmove=dt/numberofmoves;
    for (movecount=0;movecount<numberofmoves;movecount++) {
        //#pragma omp parallel for default(shared) private(e, mxold, mzold)
        for (e=0;e<N;e++) {
            //flytning og opdatering af celleinfo
            mxold=mx[e]; mzold=mz[e];
            s=mcs[e]; t=mct[e];
            mx[e]=mx[e]+dtmove*vxtomarker(e); mz[e]=mz[e]+dtmove*vztomarker(e);
            if (mx[e]>xmax) { mx[e]=2*xmax-mx[e]; warning++; } //den nye position er spejlingen af overskridelsen af xmax: xmax-(mx[e]-xmax)
            if (mx[e]<xmin) { mx[e]=2*xmin-mx[e]; warning++; }
            if (mz[e]>zmax) { mz[e]=2*zmax-mz[e]; warning++; }
            if (mz[e]<zmin) { mz[e]=2*zmin-mz[e]; warning++; }
       /*     
            if (mx[e]>=x[s] && mx[e]<=x[s+1]) { } //markeren forbliver i sin nuv�rende celle
            else { //markeren skal flyttes
                if (mx[e]>mxold) { //celleindex for�ges
                    s++;
                    while (mx[e]>x[s+1]) { //celleindex for�ges indtil man er i den rigtige celle
                        s++; } }
                if (mx[e]<mxold) { //celleindex formindskes
                    s--;
                    while (mx[e]<x[s]) { //celleindex formindskes indtil man er i den rigtige celle
                        s--; } }
            }
            if (mz[e]>=z[t] && mz[e]<=z[t+1]) { } //markeren forbliver i sin nuv�rende celle
            else { //markeren skal flyttes
                if (mz[e]>mzold) { //celleindex for�ges
                    t++;
                    while (mz[e]>z[t+1]) { //celleindex for�ges indtil man er i den rigtige celle
                        t++; } }
                if (mz[e]<mzold) { //celleindex formindskes
                    t--;
                    while (mz[e]<z[t]) { //celleindex formindskes indtil man er i den rigtige celle
                        t--; } }
            }
            mcs[e]=s;
            mct[e]=t; */
            mcs[e]=get_s(mx[e]);
            mct[e]=get_t(mz[e]);
        }
    }
    if (warning>0) {
        fprintf(statusfile,"Markers er 'sparket' ud af griddet %d gange - reducer evt. dt\n",warning);
    }
    fprintf(statusfile,"Number of timesteps used to move markers: %d\n",numberofmoves);
}

void movemarkers1() { // using open boundaries
    double optimaldt,dtmove,mxold,mzold;
    double xstrain,zstrain;
    int movecount,numberofmoves,warning,e;
    int localsleft,localsright;
    double *movedsurfacex,*movedsurfacez;
    movedsurfacex=calloc(Sbase*Nsurfaces,sizeof(double));
    movedsurfacez=calloc(Sbase*Nsurfaces,sizeof(double));
    for (sur=0;sur<Nsurfaces;sur++) {
        for (s=0;s<Sbase;s++){
            movedsurfacex[INsur(s,sur)]=x[s]+vxtopos(x[s],surfaces[INsur(s,sur)])*dt;
            movedsurfacez[INsur(s,sur)]=surfaces[INsur(s,sur)]+vztopos(x[s],surfaces[INsur(s,sur)])*dt; 
        } }
    
    warning=0;
    optimaldt=mincellsize/maxvelocity();
    numberofmoves=ceil(dt/optimaldt);
    dtmove=dt/numberofmoves;
    xmin+=leftboundvel*dt; xmax+=rightboundvel*dt; zmin+=upboundvel*dt; zmax+=downboundvel*dt;
    for (movecount=0;movecount<numberofmoves;movecount++) {
        #pragma omp parallel for schedule(static) private(e) reduction(+:warning)
        for (e=0;e<N;e++) {
            //flytning og opdatering af celleinfo
            mx[e]=mx[e]+dtmove*vxtomarker(e); mz[e]=mz[e]+dtmove*vztomarker(e);
            if (mx[e]>xmax) { mx[e]=2*xmax-mx[e]; warning++; } //den nye position er spejlingen af overskridelsen af xmax: xmax-(mx[e]-xmax)
            if (mx[e]<xmin) { mx[e]=2*xmin-mx[e]; warning++; }
            if (mz[e]>zmax) { mz[e]=2*zmax-mz[e]; warning++; }
            if (mz[e]<zmin) { mz[e]=2*zmin-mz[e]; warning++; }
            //mcs[e]=get_s(mx[e]); mct[e]=get_t(mz[e]);
        }
    }
    //regridding
    xstrain=(xmax-xmin)/(x[Sbase-1]-x[0]); zstrain=(zmax-zmin)/(z[Tbase-1]-z[0]); fprintf(statusfile,"xstrain %g zstrain %g\n",xstrain,zstrain);
    for (s=0;s<Sbase-1;s++) dx[s]=xstrain*dx[s];
    for (t=0;t<Tbase-1;t++) dz[t]=zstrain*dz[t];
    
    //test - for getting a better value of xmax - can also be put 5 lines below
    xmax=xmin; zmax=zmin;
    for (s=0;s<Sbase-1;s++) xmax+=dx[s];
    for (t=0;t<Tbase-1;t++) zmax+=dz[t];
    
    x[0]=xmin; z[0]=zmin; x[Sbase-1]=xmax; z[Tbase-1]=zmax;
    for (s=0;s<Sbase-2;s++) x[s+1]=x[s]+dx[s];
    for (t=0;t<Tbase-2;t++) z[t+1]=z[t]+dz[t];
    for (G=1;G<=multigridlevels;G++){
        S=(Sbase-1)/pow2[G]+1;
        T=(Tbase-1)/pow2[G]+1;
        for (s=0;s<S-1;s++) {
            dx[INdx(s,G)]=dx[INdx(2*s,G-1)]+dx[INdx(2*s+1,G-1)]; }
        for (t=0;t<T-1;t++) {
            dz[INdz(t,G)]=dz[INdz(2*t,G-1)]+dz[INdz(2*t+1,G-1)]; } }
    #pragma omp parallel for schedule(static) private(e)
    for (e=0;e<N;e++) {
        mcs[e]=get_s(mx[e]);
        mct[e]=get_t(mz[e]);
    }
    
    if (warning>0) {
        fprintf(statusfile,"Markers er 'sparket' ud af griddet %d gange - reducer evt. dt\n",warning);
    }
    fprintf(statusfile,"Number of timesteps used to move markers: %d\n",numberofmoves);
    
    for (sur=0;sur<Nsurfaces;sur++) {
        for (s=1;s<Sbase-1;s++){
            localsleft=s; 
            while (x[s]<movedsurfacex[INsur(localsleft,sur)]) {
                if (localsleft>0) localsleft--; else break;
            }
            localsright=s;
            while (x[s]>movedsurfacex[INsur(localsright,sur)]) {
                if (localsright<Sbase) localsright++; else break;
            }
            if (localsleft==localsright) {
                if (localsleft==0) localsright++; else localsleft--;
            }
            surfaces[INsur(s,sur)]=(movedsurfacez[INsur(localsleft,sur)]-movedsurfacez[INsur(localsright,sur)])/(movedsurfacex[INsur(localsleft,sur)]-movedsurfacex[INsur(localsright,sur)])*(x[s]-movedsurfacex[INsur(localsleft,sur)])+movedsurfacez[INsur(localsleft,sur)];
        }
        surfaces[INsur(0,sur)]=movedsurfacez[INsur(0,sur)];
        surfaces[INsur(Sbase-1,sur)]=movedsurfacez[INsur(Sbase-1,sur)];
    }
    
    for (sur=0;sur<Nsurfaces;sur++) {
        for (s=0;s<Sbase;s++){
            if (surfaces[INsur(s,sur)]>zmax) surfaces[INsur(s,sur)]=2*zmax-surfaces[INsur(s,sur)];
            if (surfaces[INsur(s,sur)]<zmin) surfaces[INsur(s,sur)]=2*zmin-surfaces[INsur(s,sur)];
        }}
    
    free(movedsurfacex);
    free(movedsurfacez);
}


void resetbasegrid(double *grid) {
    for (s=0;s<Sbase;s++) { for (t=0;t<Tbase;t++) { grid[IN0(s,t)]=0; } }
}

double minmaxwelltime(double *gridviscn,double *gridviscs, double *gridmyn, double *gridmys) {
    double output=1e100;
    for (s=0;s<Sbase;s++) {
        for (t=0;t<Tbase;t++) {
            st=IN0(s,t);
            if (output>gridviscn[st]/gridmyn[st]) output=gridviscn[st]/gridmyn[st];
            if (output>gridviscs[st]/gridmys[st]) output=gridviscs[st]/gridmys[st];
        } }
    return output;
}

int getnearestmarker(double xpos,double zpos) {
    int nearestmarkerglob,e;
    double mindist2glob=1e100;
    #pragma omp parallel private(e)
 {
        int nearestmarker;
        double dist2;
        double mindist2=1e100;
        #pragma omp for schedule(static)
        for (e=0;e<N;e++) {
            dist2=(xpos-mx[e])*(xpos-mx[e])+(zpos-mz[e])*(zpos-mz[e]);
            if (dist2<mindist2) {mindist2=dist2; nearestmarker=e; }
        }
//         printf("thread: %d nearest marker: %d dist: %f\n", omp_get_thread_num(), nearestmarker, sqrt(mindist2));
        #pragma omp critical 
        {
            if (mindist2<mindist2glob) {mindist2glob=mindist2; nearestmarkerglob=nearestmarker; }
        }
    }
    return nearestmarkerglob;
}

// int getnearestmarker(double xpos,double zpos) {
//     int nearestmarker,e;
//     double dist2;
//     double mindist2=1e100;
//     for (e=0;e<N;e++) {
//         dist2=(xpos-mx[e])*(xpos-mx[e])+(zpos-mz[e])*(zpos-mz[e]);
//         if (dist2<mindist2) {mindist2=dist2; nearestmarker=e; }
//     }
//     return nearestmarker;
// }

void markermomentumparamstogrid() {
    double *Zn,*Zs,*weightviscn,*weightviscs,*weightmyn,*weightmys,*weightoldstressxx,*weightoldstressxz,*weightrho; //these values will not get allocated memory for themselves since they use memory that is not beeing used during this step
    double dtmove;
    int warning,e,s,t;
    dt=dtnom;
    dtmove=4*mincellsize/maxvelocity(); //The time step needed to move the markers through a significant part of the smallest cell with the maximum velocity from previous time step
    if (dtmove<dt) dt=dtmove;
    //if (dtmaxwell<dt) dt=dtmaxwell; //husk at s�tte f�rste dtmaxwell=dtnom
    if (dtnom<dt) dt=dtnom;
    //firstmarkerupdate=1; // skal denne bruges?
    
    calculatemarkerparams(1);
    weightviscn=cx2; weightviscs=cx3; weightmyn=cx4; weightmys=cx5; weightoldstressxx=cx6; weightoldstressxz=cx7; weightrho=cx8;
    //resetting grids and weights
    resetbasegrid(gridviscn); resetbasegrid(gridviscs); resetbasegrid(weightviscn); resetbasegrid(weightviscs); resetbasegrid(gridmys); resetbasegrid(gridmyn); resetbasegrid(weightmys); resetbasegrid(weightmyn); resetbasegrid(gridrho); resetbasegrid(weightrho); resetbasegrid(gridoldstressxx); resetbasegrid(weightoldstressxx); resetbasegrid(gridoldstressxz); resetbasegrid(weightoldstressxz);
    
//#pragma omp parallel for default(shared) private(e)
    #pragma omp parallel private(e)
    {
        int Ncells, Nt, Nblocks, threadmustworkatthisblock, localblock, odd;
        Ncells=Sbase-1; Nt=omp_get_num_threads(); Nblocks=2*Nt;
        
        for (odd=0;odd<2;odd++) {
            #pragma omp barrier
                    threadmustworkatthisblock=2*omp_get_thread_num()+odd;
//             printf("thread: %d Nblocks: %d This block: %d\n", omp_get_thread_num(), Nblocks, threadmustworkatthisblock);
            for (e=0;e<N;e++) { //getting parameters from each marker
                localblock=(mcs[e]*Nblocks)/Ncells;
                if (localblock==threadmustworkatthisblock) {
                    markervisctogrid(e, mvisc[e], gridviscs, gridviscn, weightviscs, weightviscn);
                    markervisctogrid(e, 1/mmy[e], gridmys, gridmyn, weightmys, weightmyn); //for my, harmonic mean is used
                    markertorhogrid(e, mrho[e], gridrho, weightrho);
                    markertorhogrid(e, mstressxz[e], gridoldstressxz, weightoldstressxz);
                    markertoPgrid1(e, mstressxx[e], gridoldstressxx, weightoldstressxx);
                } } }
    }
    warning=0;
    #pragma omp parallel for schedule(static) private(s,t,st,e) reduction(+:warning)
    for (t=0;t<Tbase;t++) {
        for (s=0;s<Sbase;s++) { //reciprocal value of my to obtain harmonic mean
            st=IN0(s,t);
            if (weightviscs[st]>0) {
                gridviscs[st]/=weightviscs[st]; gridmys[st]=weightmys[st]/gridmys[st]; gridrho[st]/=weightrho[st]; gridoldstressxz[st]/=weightoldstressxz[st]; }
            else {
                e=getnearestmarker(x[s],z[t]); gridviscs[st]=mvisc[e]; gridmys[st]=mmy[e]; gridrho[st]=mrho[e]; gridoldstressxz[st]=mstressxz[e]; warning++; }
        }}
    if (warning>0) fprintf(statusfile,"Warning: %d rho-points lack markers. Nearest neighbor interp. performed\n",warning);
    warning=0;
    #pragma omp parallel for schedule(static) private(s,t,st,e) reduction(+:warning)
    for (t=0;t<Tbase-1;t++) {
        for (s=0;s<Sbase-1;s++) { //reciprocal value of my to obtain harmonic mean
            st=IN0(s,t);
            if (weightviscn[st]>0) {
                gridviscn[st]/=weightviscn[st]; gridmyn[st]=weightmyn[st]/gridmyn[st]; gridoldstressxx[st]/=weightoldstressxx[st]; }
            else {
                e=getnearestmarker(x[s]+0.5*dx[s],z[t]+0.5*dz[t]); gridviscn[st]=mvisc[e]; gridmyn[st]=mmy[e]; gridoldstressxx[st]=mstressxx[e]; warning++; }
        }}
    if (warning>0) fprintf(statusfile,"Warning: %d P-points lack markers. Nearest neighbor interp. performed\n",warning);
    
    dtmaxwell=minmaxwelltime(gridviscn,gridviscs,gridmyn,gridmys); //this is the characteristic maxwell time for the next timestep
    
    //calculation of coeffecients
    Zn=cz0; Zs=cz1; //the memory from cz0 and cz1 is used
    #pragma omp parallel for schedule(static) private(s,t,st)
    for (t=0;t<Tbase;t++)	{
        for (s=0;s<Sbase;s++){
            st=IN0(s,t);
            if (s<Sbase-1 && t<Tbase-1) Zn[st]=dt*gridmyn[st]/(dt*gridmyn[st]+gridviscn[st]);
            Zs[st]=dt*gridmys[st]/(dt*gridmys[st]+gridviscs[st]);
            //Zn[st]=1; Zs[st]=1; //skal slettes igen
        }}
    
    #pragma omp parallel for schedule(static) private(s,t,st)
    for (t=0;t<Tbase;t++)	{
        for (s=0;s<Sbase;s++){
            st=IN0(s,t);
            if (s<Sbase-1 && t<Tbase-1) inZn[st]=gridviscn[st]*Zn[st];  //"input" viscosities
            inZs[st]=gridviscs[st]*Zs[st]; //"input" viscosities
            if (s==0 || s==Sbase-1 || t==0 || t==Tbase-2 || t==Tbase-1) { }
            else {
                rightx[st]=-2*(gridoldstressxx[st]*(1-Zn[st])-gridoldstressxx[IN0(s-1,t)]*(1-Zn[IN0(s-1,t)]))/(dx[s]+dx[s-1])-(gridoldstressxz[IN0(s,t+1)]*(1-Zs[IN0(s,t+1)])-gridoldstressxz[st]*(1-Zs[st]))/dz[t]-0.5*(gridrho[st]+gridrho[IN0(s,t+1)])*gx;
            }
            if (t==0 || t==Tbase-1 || s==0 || s==Sbase-2 || s==Sbase-1) { }
            else {
                rightz[st]=-2*(-gridoldstressxx[st]*(1-Zn[st])+gridoldstressxx[IN0(s,t-1)]*(1-Zn[IN0(s,t-1)]))/(dz[t]+dz[t-1])-(gridoldstressxz[IN0(s+1,t)]*(1-Zs[IN0(s+1,t)])-gridoldstressxz[st]*(1-Zs[st]))/dx[s]-0.5*(gridrho[st]+gridrho[IN0(s+1,t)])*gz;
            } } }
}

void stressrotation() {
    int s,t,e;
    double *omega,localomega;
    double tmpstressxx,tmpstressxz;
    double vzV,vzE,vxN,vxS;
    double cos1,sin1,cos2,sin2; // normal cos and sin and cos and sin to double angle
    omega=cx0; //use the memory of cx0 witch is unused
    resetbasegrid(omega); //the rotation of the boundaries is 0

    #pragma omp parallel for schedule(static) private(s, t, st, vzV, vzE, vxN, vxS)
    for (t=1;t<Tbase-1;t++)	{
        for (s=1;s<Sbase-1;s++){
            st=IN0(s,t); vzE=vz[st]; vzV=vz[IN0(s-1,t)]; vxS=vx[st]; vxN=vx[IN0(s,t-1)];
            omega[st]=(vxS-vxN)/(dz[t]+dz[t-1])-(vzE-vzV)/(dx[s]+dx[s-1]); } }
    for (s=1;s<Sbase-1;s++){ omega[IN0(s,0)]=omega[IN0(s,1)]; omega[IN0(s,Tbase-1)]=omega[IN0(s,Tbase-2)]; }
    for (t=1;t<Tbase-1;t++){ omega[IN0(0,t)]=omega[IN0(1,t)]; omega[IN0(Sbase-1,t)]=omega[IN0(Sbase-2,t)]; }
    omega[IN0(0,0)]=0.5*(omega[IN0(1,0)]+omega[IN0(0,1)]);
    omega[IN0(Sbase-1,0)]=0.5*(omega[IN0(Sbase-2,0)]+omega[IN0(Sbase-1,1)]);
    omega[IN0(0,Tbase-1)]=0.5*(omega[IN0(1,Tbase-1)]+omega[IN0(0,Tbase-2)]);
    omega[IN0(Sbase-1,Tbase-1)]=0.5*(omega[IN0(Sbase-2,Tbase-1)]+omega[IN0(Sbase-1,Tbase-2)]);
    
    #pragma omp parallel for default(shared) private(e, localomega, cos1, sin1, cos2, sin2, tmpstressxx, tmpstressxz)
    for (e=0;e<N;e++) { //source: http://www.efunda.com/formulae/solid_mechanics/mat_mechanics/plane_stress.cfm
        localomega=rhogridtomarker(e,omega);
        cos1=cos(dt*localomega); sin1=sin(dt*localomega);
        sin2=2*sin1*cos1; //double angle formula sin(2a)=2sin(a)cos(a)
        cos2=2*cos1*cos1-1; //double angle formula cos(2a)=2cos^2(a)-1
        tmpstressxx=mstressxx[e]*cos2+mstressxz[e]*sin2;
        tmpstressxz=-mstressxx[e]*sin2+mstressxz[e]*cos2;   
        mstressxx[e]=tmpstressxx;
        mstressxz[e]=tmpstressxz;
    }
}

void stressstrainupdate() {
    int s,t,e;
    double *gridstrainxx,*gridstrainxz,*gridstressxx,*gridstressxz,*stressdiffxx,*stressdiffxz,*stressdiffxxA,*stressdiffxzA,*weightstressxx,*weightstressxz;
    double vzV,vzE,vxN,vxS;
    double tstress,stress1; //sltess
    double mtstress,m1stress;// check article for nomenclature
    double Z,expfactor,d;
    //double localstressxx,localstrainxx,localstressdiffxx; //help variables for help with calculating Hr
    d=1;
    gridstrainxx=cx0; gridstrainxz=cx1; gridstressxx=cx2; gridstressxz=cx3; stressdiffxx=cx4; stressdiffxz=cx5;
    gridHs=resx;
    resetbasegrid(gridstrainxx); resetbasegrid(gridstrainxz);
    resetbasegrid(gridstressxx); resetbasegrid(gridstressxz);
    resetbasegrid(stressdiffxx); resetbasegrid(stressdiffxz);
    resetbasegrid(gridHs);
    
    #pragma omp parallel for schedule(static) private(s,t,st)
    for (t=0;t<Tbase-1;t++)	{
        for (s=0;s<Sbase-1;s++){
            st=IN0(s,t);
            gridstrainxx[st]=0.5*((vx[IN0(s+1,t)]-vx[st])/dx[s]-(vz[IN0(s,t+1)]-vz[st])/dz[t]);
            Z=dt*gridmyn[st]/(dt*gridmyn[st]+gridviscn[st]); //Z=1;
            gridstressxx[st]=2*gridviscn[st]*gridstrainxx[st]*Z+gridoldstressxx[st]*(1-Z);
            stressdiffxx[st]=gridstressxx[st]-gridoldstressxx[st];
        }
    }
    
    #pragma omp parallel for schedule(static) private(s,t,st,vzE,vzV,vxS,vxN)
    for (t=1;t<Tbase-1;t++)	{
        for (s=1;s<Sbase-1;s++){
            st=IN0(s,t); vzE=vz[st]; vzV=vz[IN0(s-1,t)]; vxS=vx[st]; vxN=vx[IN0(s,t-1)];
            gridstrainxz[st]=(vxS-vxN)/(dz[t]+dz[t-1])+(vzE-vzV)/(dx[s]+dx[s-1]); } }
    for (s=1;s<Sbase-1;s++){ gridstrainxz[IN0(s,0)]=gridstrainxz[IN0(s,1)]; gridstrainxz[IN0(s,Tbase-1)]=gridstrainxz[IN0(s,Tbase-2)]; }
    for (t=1;t<Tbase-1;t++){ gridstrainxz[IN0(0,t)]=gridstrainxz[IN0(1,t)]; gridstrainxz[IN0(Sbase-1,t)]=gridstrainxz[IN0(Sbase-2,t)]; }
    gridstrainxz[IN0(0,0)]=0.5*(gridstrainxz[IN0(1,0)]+gridstrainxz[IN0(0,1)]);
    gridstrainxz[IN0(Sbase-1,0)]=0.5*(gridstrainxz[IN0(Sbase-2,0)]+gridstrainxz[IN0(Sbase-1,1)]);
    gridstrainxz[IN0(0,Tbase-1)]=0.5*(gridstrainxz[IN0(1,Tbase-1)]+gridstrainxz[IN0(0,Tbase-2)]);
    gridstrainxz[IN0(Sbase-1,Tbase-1)]=0.5*(gridstrainxz[IN0(Sbase-2,Tbase-1)]+gridstrainxz[IN0(Sbase-1,Tbase-2)]);
    
    #pragma omp parallel for schedule(static) private(s, t, st)
    for (t=0;t<Tbase;t++)	{
        for (s=0;s<Sbase;s++){
            st=IN0(s,t);
            Z=dt*gridmys[st]/(dt*gridmys[st]+gridviscs[st]); //Z=1;
            gridstressxz[st]=2*gridviscs[st]*gridstrainxz[st]*Z+gridoldstressxz[st]*(1-Z);
            stressdiffxz[st]=gridstressxz[st]-gridoldstressxz[st]; } }
    //#pragma omp parallel for schedule(static, 65) default(shared) private(s, t, st)
    
    #pragma omp parallel for schedule(static) private(s, t, st)
    for (t=1;t<Tbase-1;t++)	{
        for (s=1;s<Sbase-1;s++){
            st=IN0(s,t); //Hs=sigmaxx^2/visc+sigmaxz^2/visc
            gridHs[st]=0.25*(gridstressxx[st]*gridstressxx[st]/gridviscn[st]+gridstressxx[IN0(s-1,t)]*gridstressxx[IN0(s-1,t)]/gridviscn[IN0(s-1,t)]+gridstressxx[IN0(s,t-1)]*gridstressxx[IN0(s,t-1)]/gridviscn[IN0(s,t-1)]+gridstressxx[IN0(s-1,t-1)]*gridstressxx[IN0(s-1,t-1)]/gridviscn[IN0(s-1,t-1)])+gridstressxz[st]*gridstressxz[st]/gridviscs[st];
        }}
    weightstressxx=cx6; weightstressxz=cx7; stressdiffxxA=cx8; stressdiffxzA=cx9;
    resetbasegrid(stressdiffxxA); resetbasegrid(stressdiffxzA);
    resetbasegrid(weightstressxx); resetbasegrid(weightstressxz);
    //#pragma omp parallel for default(shared) private(e,mtstress,m1stress)
    
    #pragma omp parallel private(e, mtstress, m1stress)
    {
        int Ncells, Nt, Nblocks, threadmustworkatthisblock, localblock, odd;
        Ncells=Sbase-1; Nt=omp_get_num_threads(); Nblocks=2*Nt;
        
        for (odd=0;odd<2;odd++) {
            #pragma omp barrier
                    threadmustworkatthisblock=2*omp_get_thread_num()+odd;
//             printf("thread: %d Nblocks: %d This block: %d\n", omp_get_thread_num(), Nblocks, threadmustworkatthisblock);
            for (e=0;e<N;e++) {
                localblock=(mcs[e]*Nblocks)/Ncells;
                if (localblock==threadmustworkatthisblock) {
                    mtstress=mstressxx[e]+Pgridtomarker(e, stressdiffxx);
                    m1stress=Pgridtomarker(e, gridstressxx);
                    mstressxx[e]=m1stress-(m1stress-mtstress)*exp(-d*dt*mmy[e]/mvisc[e]); //tstressD
                    markertoPgrid1(e, mstressxx[e]-mtstress, stressdiffxxA, weightstressxx); //deltastressA. NOTICE function number 1
                    mtstress=mstressxz[e]+rhogridtomarker(e, stressdiffxz);
                    m1stress=rhogridtomarker(e, gridstressxz);
                    mstressxz[e]=m1stress-(m1stress-mtstress)*exp(-d*dt*mmy[e]/mvisc[e]); //tstressD
                    markertorhogrid(e, mstressxz[e]-mtstress, stressdiffxzA, weightstressxz); //deltastressA
                }}
        } }

    #pragma omp parallel for schedule(static) private(s, t, st)
    for (t=0;t<Tbase;t++)	{
        for (s=0;s<Sbase;s++){
            st=IN0(s, t);
            if (weightstressxx[st]>0) stressdiffxxA[st]/=weightstressxx[st];
            if (weightstressxz[st]>0) stressdiffxzA[st]/=weightstressxz[st]; } } //it is not a big problem if a node is not assigned a weight. the correction will then just be zero
   
    #pragma omp parallel for schedule(static) private(e)
    for (e=0;e<N;e++) {
        mstressxx[e]=mstressxx[e]-Pgridtomarker(e,stressdiffxxA);
        mstressxz[e]=mstressxz[e]-rhogridtomarker(e,stressdiffxzA); //final correction
        //new that the below is in the same loop as the above
        mstrainxx[e]=Pgridtomarker(e, gridstrainxx);
        mstrainxz[e]=rhogridtomarker(e, gridstrainxz);
        if (mfinitestrain[e]>=0) mfinitestrain[e]+=dt*sqrt(secdevinv2(mstrainxx[e]-Pgridtomarker(e, stressdiffxx)/(2*dt*mmy[e]), mstrainxz[e]-rhogridtomarker(e, stressdiffxz)/(2*dt*mmy[e]))); //the first time strain will be remembered, but this should not matter much since the first time step is very small
        mP[e]=Pgridtomarker(e, P);
    }

    stressintegral=0;
    for (t=0;t<Tbase-1;t++)	{
        st=IN0(Sbase-2,t);
        stressintegral+=gridstressxx[st]*dz[t]; }
    fprintf(statusfile,"stress integral %g\n",stressintegral);
}

double stressstrainupdate_readonly() {
    int s, t, st;
    double *gridstrainxx, *gridstrainxz, *gridstressxx, *gridstressxz;
    double vzV, vzE, vxN, vxS;
    double Z;
    double stressintegral;
    double d = 1;
    
    // Allocate temporary grids
    gridstrainxx = cx0;
    gridstrainxz = cx1;
    gridstressxx = cx2;
    gridstressxz = cx3;
    
    resetbasegrid(gridstrainxx);
    resetbasegrid(gridstrainxz);
    resetbasegrid(gridstressxx);
    resetbasegrid(gridstressxz);
    
    // Calculate gridstrainxx and gridstressxx
    #pragma omp parallel for schedule(static) private(s, t, st, Z)
    for (t = 0; t < Tbase - 1; t++) {
        for (s = 0; s < Sbase - 1; s++) {
            st = IN0(s, t);
            gridstrainxx[st] = 0.5 * ((vx[IN0(s + 1, t)] - vx[st]) / dx[s] - 
                                      (vz[IN0(s, t + 1)] - vz[st]) / dz[t]);
            Z = dt * gridmyn[st] / (dt * gridmyn[st] + gridviscn[st]);
            gridstressxx[st] = 2 * gridviscn[st] * gridstrainxx[st] * Z + 
                               gridoldstressxx[st] * (1 - Z);
        }
    }
    
    // Calculate gridstrainxz
    #pragma omp parallel for schedule(static) private(s, t, st, vzE, vzV, vxS, vxN)
    for (t = 1; t < Tbase - 1; t++) {
        for (s = 1; s < Sbase - 1; s++) {
            st = IN0(s, t);
            vzE = vz[st];
            vzV = vz[IN0(s - 1, t)];
            vxS = vx[st];
            vxN = vx[IN0(s, t - 1)];
            gridstrainxz[st] = (vxS - vxN) / (dz[t] + dz[t - 1]) + 
                               (vzE - vzV) / (dx[s] + dx[s - 1]);
        }
    }
    // Fill boundaries
    for (s = 1; s < Sbase - 1; s++) {
        gridstrainxz[IN0(s, 0)] = gridstrainxz[IN0(s, 1)];
        gridstrainxz[IN0(s, Tbase - 1)] = gridstrainxz[IN0(s, Tbase - 2)];
    }
    for (t = 1; t < Tbase - 1; t++) {
        gridstrainxz[IN0(0, t)] = gridstrainxz[IN0(1, t)];
        gridstrainxz[IN0(Sbase - 1, t)] = gridstrainxz[IN0(Sbase - 2, t)];
    }
    gridstrainxz[IN0(0, 0)] = 0.5 * (gridstrainxz[IN0(1, 0)] + gridstrainxz[IN0(0, 1)]);
    gridstrainxz[IN0(Sbase - 1, 0)] = 0.5 * (gridstrainxz[IN0(Sbase - 2, 0)] + gridstrainxz[IN0(Sbase - 1, 1)]);
    gridstrainxz[IN0(0, Tbase - 1)] = 0.5 * (gridstrainxz[IN0(1, Tbase - 1)] + gridstrainxz[IN0(0, Tbase - 2)]);
    gridstrainxz[IN0(Sbase - 1, Tbase - 1)] = 0.5 * (gridstrainxz[IN0(Sbase - 2, Tbase - 1)] + gridstrainxz[IN0(Sbase - 1, Tbase - 2)]);
    
    // Calculate gridstressxz
    #pragma omp parallel for schedule(static) private(s, t, st, Z)
    for (t = 0; t < Tbase; t++) {
        for (s = 0; s < Sbase; s++) {
            st = IN0(s, t);
            Z = dt * gridmys[st] / (dt * gridmys[st] + gridviscs[st]);
            gridstressxz[st] = 2 * gridviscs[st] * gridstrainxz[st] * Z + 
                               gridoldstressxz[st] * (1 - Z);
        }
    }
    
    // Calculate and return stress integral
    stressintegral = 0;
    for (t = 0; t < Tbase - 1; t++) {
        st = IN0(Sbase - 2, t);
        stressintegral += gridstressxx[st] * dz[t];
    }
    
    return stressintegral;
}

void printages(int name) { 
    FILE *file;
    char a[100];
    sprintf(a,"outfiles/age%d.txt",name);
    file = fopen(a,"wb");
    fwrite(mage,sizeof(double),N,file);
    fclose(file);
}

void printX(int name) { 
    FILE *file;
    char a[100];
    sprintf(a,"outfiles/X%d.txt",name);
    file = fopen(a,"wb");
    fwrite(mX,sizeof(double),N,file);
    fwrite(mdm,sizeof(double),N,file);
    fwrite(mXComp[1],sizeof(double),N,file);
    fwrite(mRComp[1],sizeof(double),N,file);
    fclose(file);
}

void printstate(int name) { 
    FILE *file;
    char a[100];
    sprintf(a,"outfiles/%d.txt",name);
    file = fopen(a,"wb");
    fwrite(&N,sizeof(int),1,file);
    fwrite(&Sbase,sizeof(int),1,file);
    fwrite(&Tbase,sizeof(int),1,file);
    fwrite(&multigridlevels,sizeof(int),1,file);
    fwrite(&leftbound,sizeof(int),1,file);
    fwrite(&rightbound,sizeof(int),1,file);
    fwrite(&upbound,sizeof(int),1,file);
    fwrite(&downbound,sizeof(int),1,file);
    fwrite(&time,sizeof(double),1,file);
    fwrite(&dtnom,sizeof(double),1,file); dtmaxwell=dtnom;
    fwrite(dx,sizeof(double),Sbase-1,file);
    fwrite(dz,sizeof(double),Tbase-1,file);
    fwrite(mx,sizeof(double),N,file);
    fwrite(mz,sizeof(double),N,file);
    fwrite(mcs,sizeof(int),N,file);
    fwrite(mct,sizeof(int),N,file);
    fwrite(mtype,sizeof(int),N,file);
    fwrite(mTemp,sizeof(double),N,file);
    fwrite(mstressxx,sizeof(double),N,file);
    fwrite(mstressxz,sizeof(double),N,file);
    fwrite(mstrainxx,sizeof(double),N,file);
    fwrite(mstrainxz,sizeof(double),N,file);
    fwrite(mfinitestrain,sizeof(double),N,file);
    fwrite(vx,sizeof(double),Sbase*Tbase,file);
    fwrite(vz,sizeof(double),Sbase*Tbase,file);
    fwrite(P,sizeof(double),Sbase*Tbase,file);
    fwrite(&Nsurfaces,sizeof(int),1,file);
    fwrite(surfaces,sizeof(double),Nsurfaces*Sbase,file);

    fclose(file);
    
    printages(name);
    printX(name);
}

/*void stressstrainupdate() {
    double *gridstrainxx,*gridstrainxz,*gridstressxx,*gridstressxz,*stressdiffxx,*stressdiffxz,*stressdiffxxA,*stressdiffxzA,*weightstressxx,*weightstressxz;
    double vzV,vzE,vxN,vxS;
    double tstress,stress1;
    double mtstress,m1stress;// check article for nomenclature
    double Z,expfactor,d;
    d=0;
    gridstrainxx=cx0; gridstrainxz=cx1; gridstressxx=cx2; gridstressxz=cx3; stressdiffxx=cx4; stressdiffxz=cx5;
    resetbasegrid(gridstrainxx); resetbasegrid(gridstrainxz);
    resetbasegrid(gridstressxx); resetbasegrid(gridstressxz);
    resetbasegrid(stressdiffxx); resetbasegrid(stressdiffxz);
    for (s=0;s<Sbase-1;s++){
        for (t=0;t<Tbase-1;t++)	{
            st=IN0(s,t);
            gridstrainxx[st]=0.5*((vx[IN0(s+1,t)]-vx[st])/dx[s]-(vz[IN0(s,t+1)]-vz[st])/dz[t]);
            Z=dt*gridmyn[st]/(dt*gridmyn[st]+gridviscn[st]); //Z=1;
            gridstressxx[st]=2*gridviscn[st]*gridstrainxx[st]*Z+gridoldstressxx[st]*(1-Z);
            stressdiffxx[st]=gridstressxx[st]-gridoldstressxx[st];
        }
    }
    for (s=1;s<Sbase-1;s++){
        for (t=1;t<Tbase-1;t++)	{
            st=IN0(s,t); vzE=vz[st]; vzV=vz[IN0(s-1,t)]; vxS=vx[st]; vxN=vx[IN0(s,t-1)];
            gridstrainxz[st]=(vxS-vxN)/(dz[t]+dz[t-1])+(vzE-vzV)/(dx[s]+dx[s-1]); } }
    for (s=1;s<Sbase-1;s++){ gridstrainxz[IN0(s,0)]=gridstrainxz[IN0(s,1)]; gridstrainxz[IN0(s,Tbase-1)]=gridstrainxz[IN0(s,Tbase-2)]; }
    for (t=1;t<Tbase-1;t++){ gridstrainxz[IN0(0,t)]=gridstrainxz[IN0(1,t)]; gridstrainxz[IN0(Sbase-1,t)]=gridstrainxz[IN0(Sbase-2,t)]; }
    gridstrainxz[IN0(0,0)]=0.5*(gridstrainxz[IN0(1,0)]+gridstrainxz[IN0(0,1)]);
    gridstrainxz[IN0(Sbase-1,0)]=0.5*(gridstrainxz[IN0(Sbase-2,0)]+gridstrainxz[IN0(Sbase-1,1)]);
    gridstrainxz[IN0(0,Tbase-1)]=0.5*(gridstrainxz[IN0(1,Tbase-1)]+gridstrainxz[IN0(0,Tbase-2)]);
    gridstrainxz[IN0(Sbase-1,Tbase-1)]=0.5*(gridstrainxz[IN0(Sbase-2,Tbase-1)]+gridstrainxz[IN0(Sbase-1,Tbase-2)]);
    for (s=0;s<Sbase;s++){
        for (t=0;t<Tbase;t++)	{
            st=IN0(s,t);
            Z=dt*gridmys[st]/(dt*gridmys[st]+gridviscs[st]); //Z=1;
            gridstressxz[st]=2*gridviscs[st]*gridstrainxz[st]*Z+gridoldstressxz[st]*(1-Z);
            stressdiffxz[st]=gridstressxz[st]-gridoldstressxz[st]; } }
    
    
    weightstressxx=cx6; weightstressxz=cx7; stressdiffxxA=cx8; stressdiffxzA=cx9;
    resetbasegrid(stressdiffxxA); resetbasegrid(stressdiffxzA);
    resetbasegrid(weightstressxx); resetbasegrid(weightstressxz);
    for (e=0;e<N;e++) {
        mtstress=mstressxx[e]+Pgridtomarker(e,stressdiffxx);
        m1stress=Pgridtomarker(e,gridstressxx);
        mstressxx[e]=m1stress-(m1stress-mtstress)*exp(-d*dt*mmy[e]/mvisc[e]); //tstressD
        markertoPgrid1(e,mstressxx[e]-mtstress,stressdiffxxA,weightstressxx); //deltastressA. NOTICE function number 1
        mtstress=mstressxz[e]+rhogridtomarker(e,stressdiffxz);
        m1stress=rhogridtomarker(e,gridstressxz);
        mstressxz[e]=m1stress-(m1stress-mtstress)*exp(-d*dt*mmy[e]/mvisc[e]); //tstressD
        markertorhogrid(e,mstressxz[e]-mtstress,stressdiffxzA,weightstressxz); //deltastressA
    }
    for (s=0;s<Sbase;s++){ for (t=0;t<Tbase;t++)	{
        st=IN0(s,t);
        if (weightstressxx[st]>0) stressdiffxxA[st]/=weightstressxx[st];
        if (weightstressxz[st]>0) stressdiffxzA[st]/=weightstressxz[st]; } } //it is not a big problem if a node is not assigned a weight. the correction will then just be zero
    for (e=0;e<N;e++) {
        mstressxx[e]=mstressxx[e]-Pgridtomarker(e,stressdiffxxA);
        mstressxz[e]=mstressxz[e]-rhogridtomarker(e,stressdiffxzA); //final correction
    }

    for (e=0;e<N;e++) {
        mstrainxx[e]=Pgridtomarker(e,gridstrainxx);
        mstrainxz[e]=rhogridtomarker(e,gridstrainxz);
        //mstressxx[e]=Pgridtomarker(e,gridstressxx);
        //mstressxz[e]=rhogridtomarker(e,gridstressxz);        
        //mstressxz[e]=1;//mstrainxz[e];
        mP[e]=Pgridtomarker(e,P); }
    
}*/


