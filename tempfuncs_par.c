#include <stdio.h>
#include <math.h>
#include <omp.h>
#include <stdlib.h>
#include "main.h"
#include "tempfuncs_par.h"
#include "momentummarkerfuncs_par.h"
#include "tempfuncs_par.h"

void tridiag(double *x, double *a, double *b, double *c, double *d, double *bdot, double *ddot, int length) {
    int i;
    double  m;
    
    bdot[0]=b[0];
    ddot[0]=d[0];
    
    for (i=1; i<length; i++) {
        m=a[i]/bdot[i-1];
        bdot[i]=b[i]-m*c[i-1];
        ddot[i]=d[i]-m*ddot[i-1]; }
    x[length-1]=ddot[length-1]/bdot[length-1];
    for (i=length-2; i>=0; i--) {
        x[i]=(ddot[i]-c[i]*x[i+1])/bdot[i];
    }
    
    return;
}

int max(int a, int b) {
    if (a<b)	{
        return(b);	}
    else {
        return(a); }
}

void gettempmarkerparams(double *mC,double *mk,double *mHr) {
    for (e=0;e<N;e++) {
        switch(mtype[e]) {
            case 0:
                mk[e]=3;
                mC[e]=800*2800;
                mHr[e]=0;
                mTemp[e]=Ttop;
                break;
	    	case 16: /* ice */
                mk[e]=3;
                mC[e]=800*2800;
                mHr[e]=0;
                mTemp[e]=Ttop;
                break;
                
            case 1: case 11:
                mk[e]=2.5;
                mC[e]=800;
                mHr[e]=1.25e-6;
                break;
                
            case 2: case 12:
                mk[e]=2.5;
                mC[e]=800;
                mHr[e]=0.75e-6;
                break;

            case 3: case 4: case 13: case 14:
                mk[e]=4.08*pow((298/(mTemp[e]+273)),0.406);
//                 mC[e]=1000;
                mC[e]=(238.64-20.013e2*pow(mTemp[e]+273,-0.5)-11.642e7*pow(mTemp[e]+273,-3))*7.1077;
                mHr[e]=0;
                break;
                
            case 5: case 15:
                mk[e]=2;
                mC[e]=800;
                mHr[e]=0;
                break;
                      
                default:
                mk[e]=2;
                mC[e]=800;
                mHr[e]=0;
        }
    }
}

double Ts_Katz(double inputP,double F,double X, double M) {
    double A1,A2,A3,B1,B2,B3,Bt1,Bt2,C1,C2,C3,D,R1,R2,dt1,dte,xs1,xs2,xse;
    double Ts,R,Fcpxout,Tsol,Tmid,Tliq,Tcpx,taperwidth,diffstep,F0,F1,F2,F3,Ts0,Ts1,Ts2,Ts3,k1,k2,a,b,t,xm,xsat,DT,P;
    A1=1085.661000; A2=132.899000; A3=-5.104000; B1=1475.000000; B2=80.000000; B3=-3.200000; Bt1=1.500000; Bt2=1.500000; C1=1780.000000; C2=45.000000; C3=-2.000000; D=0.010000; R1=0.500000; R2=0.080000; dt1=43.000000; dte=0.750000; xs1=0.120000; xs2=0.010000; xse=0.600000;
    P=inputP;
    if (P>10.0) P=10;
    if (P<0.0) P=0.0;
// Reaction coef for cpxout
    R = R1 + R2*P;
    Fcpxout=M/R;
//Dry solidus calculations based on depletion
    Tsol=A1+A2*P+A3*P*P;
    Tmid=B1+B2*P+B3*P*P;
    Tliq=C1+C2*P+C3*P*P;
// Different Ts depending on melt fraction
    if (F<0) Ts=Tsol;
    if (F>=0 && F<=Fcpxout) Ts=pow(F, 1.0/Bt1)*(Tmid-Tsol)+Tsol;
    if (F>Fcpxout && F<=1) {
        Tcpx=pow(Fcpxout, 1.0/Bt1)*(Tmid-Tsol) + Tsol;
        Ts=pow((F-Fcpxout)/(1.0-Fcpxout), 1.0/Bt2)*(Tliq-Tcpx)+Tcpx; }
    if (F>1) Ts=Tliq;
// Cubic spline interpolation around Fcpxout
    taperwidth=0.01; diffstep=1e-8;
    if (fabs(F-Fcpxout)<taperwidth) {
        F0=Fcpxout-taperwidth-diffstep; F1=Fcpxout-taperwidth; F2=Fcpxout+taperwidth; F3=Fcpxout+taperwidth+diffstep;
        Ts0=pow(F0, 1.0/Bt1)*(Tmid-Tsol)+Tsol;
        Ts1=pow(F1, 1.0/Bt1)*(Tmid-Tsol)+Tsol;
        Tcpx=pow(Fcpxout, 1.0/Bt1)*(Tmid-Tsol) + Tsol;
        Ts2=pow((F2-Fcpxout)/(1.0-Fcpxout), 1.0/Bt2)*(Tliq-Tcpx)+Tcpx;
        Ts3=pow((F3-Fcpxout)/(1.0-Fcpxout), 1.0/Bt2)*(Tliq-Tcpx)+Tcpx;
        k1=(Ts1-Ts0)/diffstep; k2=(Ts3-Ts2)/diffstep;
//method from http://en.wikipedia.org/wiki/Spline_interpolation
        a=k1*(F2-F1)-(Ts2-Ts1); b=-k2*(F2-F1)+(Ts2-Ts1);
        t=(F-F1)/(F2-F1);
        Ts=(1.0-t)*Ts1+t*Ts2+t*(1.0-t)*(a*(1.0-t)+b*t);
    }
//temperature decrease of solidus due to water content
    if (X>0) {
        xm=X/(D+F*(1-D));
        xsat = xs1*pow(P, xse)+xs2*P;
        if (xm>xsat) xm=xsat;
        DT=dt1*pow(xm*100, dte);
        Ts=Ts-DT;
    }
    return Ts;
}

double T_G2_kdp(double inputP,double F) {
    double Aa,Bb,Cc,Dd,Ee,Ff,P,Tsol,Tliq,Ts;
    P=inputP;
    if (P>10.0) P=10;
    if (P<0.0) P=0.0;
    Aa = 73.68; Bb = 26.32; Cc = 1175; Dd = 114; Ee = 920; Ff = 130;
    Tsol = Ee+Ff*P;
    Tliq = Cc+Dd*P;
    return (2*Tsol*Aa - Tliq*Bb + Tsol*Bb + Tliq*sqrt(Bb*Bb + 400*F*Aa) - Tsol*sqrt(Bb*Bb + 400*F*Aa))/(2*Aa);
}

double Ts_simple(double inputP,double F) {
	double Ts,dTdF_P,dTdP_F,Ts0;
	Ts0 = 1081;
	dTdP_F = 132e-9; /*solidus pressure gradient*/
	dTdF_P = 300; /*solidus depletion gradient spinel lherzolite*/
	Ts = Ts0 + dTdF_P*F + dTdP_F*inputP;
	return(Ts);
}

double get_Ts(double X, double press, double Xs, int type) {
  double Ts,aX,bX;
    aX = -43;
    bX = 0.75;
    if (type==3 || type==4 || type==7 || type==8) {
       	Ts=Ts_Katz(press*1e-9,(1.0-1.0/X),0,0.17);
		// Ts=1100+880*(X-1)+1.0503e-7*press;
		//Ts = Ts_simple(press,(1.0-1.0/X));
       	// effect of water on the solidus - I do it this way rather than line 107 (John)
       	Ts += aX*pow(0.01*Xs,bX);
    }

    if (type==6 || type==10)
        Ts=T_G2_kdp(press*1e-9,(1.0-1.0/X));
//         Ts=767+200e-9*press+150*(1-1/X);
//         Ts=3000;
    return Ts;
    //if (X>1.02) return Ts;
    //else return Ts-200*(1.02-X)/0.02;
}

double get_dTs_dphi(double X, double press, double Xs, double phi, int type) {
    double Ts,Tsplus,dX;
    dX=1e-8;
    Ts=get_Ts(X,press,Xs,type); Tsplus=get_Ts(X+dX,press,Xs,type);
    return (Tsplus-Ts)/dX*X/(1-phi);
    //if (X>1.02) return 440*X/(1-phi);
    //else return 1440*X/(1-phi);
}

double get_dTs_dphi2(double F, double Xs) {
    double dXsdF,dTdXs,dTdF_P,dTdF_X;
	double aX,bX,Dbulk;
	aX = -43;
	bX = 0.75;
	Dbulk = 0.01;
	dTdF_P = 300;
	if (Xs > 1e-6) {
		dXsdF = -Xs*(1/Dbulk-1)*pow((1-F),(1/Dbulk-2)); /*slope (derivative) of X(F)*/
		dTdXs = bX*aX*pow(0.01*Xs,(bX-1)); /*slope of Ts(X) */
		dTdF_X = dTdXs*dXsdF; /*combined to get slope of Ts( F(X) )*/
		dTdF_P = dTdF_P + dTdF_X;
	}
	return(dTdF_P);
}

void meltupdate() {
    int e,t; /* to stop warnings - but does it matter? */
    double *griddm, *weightdm,maxdm,integratedmelt;
    maxdm=0;
    #pragma omp parallel private(e)
    {
      double deltaS, L, dTs_dphi, Ts, C, D, phi, dm,localmaxdm,localmaxmeltpos;
        localmaxdm=0;localmaxmeltpos=0;
        #pragma omp for schedule(static)
        for (e=0;e<N;e++) {
            dm=0;
            if (mtype[e]==3 || mtype[e]==4 || mtype[e]==7 || mtype[e]==8 || mtype[e]==6 || mtype[e]==10) {
                C=(238.64-20.013e2*pow(mTemp[e]+273,-0.5)-11.642e7*pow(mTemp[e]+273,-3))*7.1077;
                phi=1.0-1.0/mX[e];
                Ts=get_Ts(mX[e], fabs(mP[e]), mWater[e], mtype[e]);
                dTs_dphi=get_dTs_dphi(mX[e], mP[e], mWater[e], phi, mtype[e]);
                //dTs_dphi=get_dTs_dphi2(1.0-1.0/mX[e], mWater[e]);
				deltaS=300;
               
                L=(mTemp[e]+273)*deltaS;
                if (mTemp[e]>Ts && phi<0.99) {
                    dm=(mTemp[e]-Ts)/(L/C+dTs_dphi);
                    mX[e]+=mX[e]/(1-phi)*dm;
                    mTemp[e]-=L/C*dm;
                    if (mtype[e]==3 || mtype[e]==4 || mtype[e]==6 || mtype[e]==9) mtype[e]+=4;
                }
                else {
                    if (mtype[e]==7 || mtype[e]==8 || mtype[e]==10 || mtype[e]==13) mtype[e]-=4;
                }
            }
	    D = 0.01; // bulk water partition coefficient
	    mWater[e] = mWater[e]/(1 + dm*(1/D - 1)); // update water
            mdm[e]=dm;
            // This step involves the calculation of of the marker with highest dm. Its lateral position is remembered.
            // I did this because I wanted to add the melt as sediments above this points (not seen here)
            // "local" refers to that each thread has finds a maximum. The "critical" pragma below identifies the global maximum.
            if (dm>localmaxdm) {localmaxdm=dm; localmaxmeltpos=mx[e]; }
        }
        #pragma omp critical
        if(localmaxdm > maxdm) { maxdm = localmaxdm; maxmeltpos=localmaxmeltpos;  }
    }
    // Here marker values of dm (mdm) are interpolated to a grid (without dividing correcting for weights yet).
    // Each marker point updates the value of the 4 grid points that corner the cell of the marker (the indices for the cell are mcs[e],mct[e])
    // Were it not parallelized, it could simply look like this:
    // for (e=0;e<N;e++) markertorhogrid(e, mdm[e], griddm, weightdm);   
    // Extra stuff (e.g. the loop with "odd") is to make sure that several threads don't write to the same grid point at the same time (possible race condition)
    // The idea with the "red/black coloring" (i.e. odd/even) is that two or more threads never work in neighboring cells
   
    griddm=cx0; weightdm=cz0; resetbasegrid(griddm); resetbasegrid(weightdm);
    #pragma omp parallel private(e)
    {
        int Ncells, Nt, Nblocks, threadmustworkatthisblock, localblock, odd;
        Ncells=Sbase-1; Nt=omp_get_num_threads(); Nblocks=2*Nt;
        for (odd=0;odd<2;odd++) {
            #pragma omp barrier
	    threadmustworkatthisblock=2*omp_get_thread_num()+odd;
            for (e=0;e<N;e++)    {
	        meltproduction = 0;
		if (mdm[e] > 0) meltproduction = mdm[e]; /* small degree melting is not pooled */
                localblock=(mcs[e]*Nblocks)/Ncells;
                if (localblock==threadmustworkatthisblock) {
                    markertorhogrid(e, meltproduction, griddm, weightdm);
                }
	    }
	}
    }
    // Here weights are divided by, finishing the interpolation to the grid
    #pragma omp parallel for schedule(static) private(s, t, st)
    for (t=0;t<Tbase;t++) {
        for (s=0;s<Sbase;s++) {
            st=IN0(s, t);
            if (weightdm[st]>0) griddm[st]/=weightdm[st]; else { griddm[st]=0; }
        }
    }
    //And here the grid is integraded
    integratedmelt=0;
    #pragma omp parallel for schedule(static) private(s, t, st) reduction(+:integratedmelt)
    for (t=0;t<Tbase-1;t++) {
        for (s=0;s<Sbase-1;s++) {
            integratedmelt+=0.25*(griddm[IN0(s, t)]+griddm[IN0(s+1, t)]+griddm[IN0(s, t+1)]+griddm[IN0(s+1, t+1)])*dx[s]*dz[t];
        }
    }
    //Conversion from mass to volume, making a gross assumption about mantle and crustal densities, respectively.
    meltproductionVOL=3300/3000*integratedmelt;
}


void updatetemp() {
    int s,t,e;
    double *a,*b,*c,*d,*x,*bdot,*ddot;
    double kplus,kcenter,kminus,alpha;
    double *gridC,*gridk,*gridoldTemp,*gridTemp,*gridTempdiff,*gridTempdiffrem,*gridTempdiffsub,*gridHr,*gridH,*gridHa,*weightC,*weightk,*weightoldTemp,*weightTemp,*weightHr;
    double standardC,standardk,standardTemp,standardHr,dt0,dimless;
    double mTempgm,mTempdiffsub;// taras chapt 11 for nomenclature
    dimless=1; // the name d is already taken
    
    gettempmarkerparams(mC,mk,mHr);
    
    gridC=cx0; gridk=cx1; gridoldTemp=cx2; gridTemp=cx3; gridTempdiff=cx4; gridHr=cx5; gridH=cx6;
    weightC=cz0; weightk=cz1; weightoldTemp=cz2; weightTemp=cz3; weightHr=cz4;
    gridHa=resz;
    resetbasegrid(gridC); resetbasegrid(gridk); resetbasegrid(gridoldTemp); resetbasegrid(gridTemp); resetbasegrid(gridTempdiff); resetbasegrid(gridHr); resetbasegrid(gridHa);
    resetbasegrid(weightC); resetbasegrid(weightk); resetbasegrid(weightoldTemp); resetbasegrid(weightTemp); resetbasegrid(weightHr);
    
    standardC=0; standardk=0; standardTemp=0; standardHr=0;
    #pragma omp parallel private(e) reduction(+:standardC,standardk,standardTemp,standardHr)
 {
        int Ncells, Nt, Nblocks, threadmustworkatthisblock, localblock, odd;
        Ncells=Sbase-1; Nt=omp_get_num_threads(); Nblocks=2*Nt;
        for (odd=0;odd<2;odd++) {
            #pragma omp barrier
                    threadmustworkatthisblock=2*omp_get_thread_num()+odd;
            for (e=0;e<N;e++) {
                localblock=(mcs[e]*Nblocks)/Ncells;
                if (localblock==threadmustworkatthisblock) {
                    standardC+=mC[e]; standardk+=mk[e]; standardTemp+=mTemp[e]; standardHr+=mHr[e];
                    markertorhogrid(e, mC[e], gridC, weightC);
                    markertorhogrid(e, mk[e], gridk, weightk);
                    markertorhogrid(e, mHr[e], gridHr, weightHr);
                    markertorhogrid(e, mTemp[e], gridoldTemp, weightoldTemp);
                } }
        }   }
    standardC/=N; standardk/=N; standardTemp/=N; standardHr/=N;
    
    #pragma omp parallel for schedule(static) private(s, t, st)
    for (t=0;t<Tbase;t++) {
        for (s=0;s<Sbase;s++) {
            st=IN0(s,t);
            if (weightC[st]>0) gridC[st]/=weightC[st]; else { gridC[st]=standardC; }
            if (weightk[st]>0) gridk[st]/=weightk[st]; else { gridk[st]=standardk; }
            if (weightHr[st]>0) gridHr[st]/=weightHr[st]; else { gridHr[st]=standardHr; }
            if (weightoldTemp[st]>0) gridoldTemp[st]/=weightoldTemp[st]; else { gridoldTemp[st]=standardTemp; }      
        }}
    
    #pragma omp parallel for schedule(static) private(s, t, st)
    for (t=1;t<Tbase-1;t++) {
        for (s=1;s<Sbase-1;s++) {
            st=IN0(s,t); //LAPPEL�SNING for term. ekspansion. Bedst at interpolere fra markers
            gridHa[st]=0.5*gridoldTemp[st]*(2.77e-5+0.97e-8*(gridoldTemp[st]+273)-0.32*pow(gridoldTemp[st]+273,-2))*gridrho[st]*((vx[st]+vx[IN0(s,t-1)])*gx+(vz[st]+vz[IN0(s-1,t)])*gz);
    }}
    #pragma omp parallel for schedule(static) private(s, t, st)
    for (t=0;t<Tbase;t++) {
        for (s=0;s<Sbase;s++) {
            st=IN0(s,t);
            gridH[st]=gridHr[st]+gridHs[st]+gridHa[st]; //add more eventually
        }}

    #pragma omp parallel private(s,t,st,kplus,kcenter,kminus,alpha,a,b,c,d,x,bdot,ddot)
    {
        a = calloc(max(Sbase, Tbase), sizeof(double));
        b = calloc(max(Sbase, Tbase), sizeof(double));
        c = calloc(max(Sbase, Tbase), sizeof(double));
        d = calloc(max(Sbase, Tbase), sizeof(double));
        x = calloc(max(Sbase, Tbase), sizeof(double));
        bdot = calloc(max(Sbase, Tbase), sizeof(double));
        ddot = calloc(max(Sbase, Tbase), sizeof(double));
        
        #pragma omp for schedule(static)
        for (s=0; s<Sbase; s++) {
            for (t=1; t<Tbase-1; t++) {
                st=IN0(s, t);
                kplus=gridk[IN0(s, t+1)]; kcenter=gridk[st]; kminus=gridk[IN0(s, t-1)];
                alpha=1/(gridrho[st]*gridC[st]*(dz[t-1]+dz[t]));
                a[t]=alpha*(kminus+kcenter)/dz[t-1];
                b[t]=-(1/dt+alpha*(kminus+kcenter)/dz[t-1]+alpha*(kplus+kcenter)/dz[t]);
                c[t]=alpha*(kplus+kcenter)/dz[t];
                d[t]=-1/dt*gridoldTemp[st]-0.5*gridH[st]/(gridrho[st]*gridC[st]); //Kildeleddet bidrager halvt for hver dimension
            }
            b[0]=1;
            c[0]=0;
            d[0]=Ttop;
            a[Tbase-1]=0;
            b[Tbase-1]=1;
            d[Tbase-1]=Tbottom;
            
            tridiag(x, a, b, c, d, bdot, ddot, Tbase);
            for (t=0; t<Tbase; t++) {
                gridTemp[IN0(s, t)]=x[t];
            }
        }

        #pragma omp for schedule(static)
        for (t=1; t<Tbase-1; t++) {
            for (s=1; s<Sbase-1; s++) {
                st=IN0(s, t);
                kplus=gridk[IN0(s+1, t)]; kcenter=gridk[st]; kminus=gridk[IN0(s-1, t)];
                alpha=1/(gridrho[st]*gridC[st]*(dx[s-1]+dx[s]));
                a[s]=alpha*(kminus+kcenter)/dx[s-1];
                b[s]=-(1/dt+alpha*(kminus+kcenter)/dx[s-1]+alpha*(kplus+kcenter)/dx[s]);
                c[s]=alpha*(kplus+kcenter)/dx[s];
                d[s]=-1/dt*gridTemp[st]-0.5*gridH[st]/(gridrho[st]*gridC[st]);
            }
            b[0]=-1;
            c[0]=1;
            d[0]=0;
            a[Sbase-1]=-1;
            b[Sbase-1]=1;
            d[Sbase-1]=0;
            
            tridiag(x, a, b, c, d, bdot, ddot, Sbase);
            for (s=0; s<Sbase; s++) {
                gridTemp[IN0(s, t)]=x[s];
            }
        }
        
        free(a);
        free(b);
        free(c);
        free(d);
        free(x);
        free(bdot);
        free(ddot);
    }
    
    gridTempdiffrem=gridHr; //gridHr is not used anymore
    gridTempdiffsub=gridHa; //gridHa is not used anymore
    resetbasegrid(weightTemp); resetbasegrid(gridTempdiffrem); resetbasegrid(gridTempdiffsub);
    
    
    #pragma omp parallel private(e,dt0,mTempgm,mTempdiffsub)
 {
        int Ncells, Nt, Nblocks, threadmustworkatthisblock, localblock, odd;
        Ncells=Sbase-1; Nt=omp_get_num_threads(); Nblocks=2*Nt;
        for (odd=0;odd<2;odd++) {
            #pragma omp barrier
                    threadmustworkatthisblock=2*omp_get_thread_num()+odd;
            for (e=0;e<N;e++) {
                localblock=(mcs[e]*Nblocks)/Ncells;
                if (localblock==threadmustworkatthisblock) {
                    dt0=mC[e]*mrho[e]/(2*mk[e]*(1/(dx[mcs[e]]*dx[mcs[e]])+1/(dz[mct[e]]*dz[mct[e]])));
                    mTempgm=rhogridtomarker(e, gridoldTemp); //mTempgm is the temperature interpolated from the markers to the grid and back to the markers. It is thus a measure of the numerical diffusion
                    mTempdiffsub=(mTempgm-mTemp[e])*(1-exp(-dimless*dt/dt0));
                    mTemp[e]+=mTempdiffsub; //a portion (depending on 'dimless') of this numerical diffusion is added to the solution
                    markertorhogrid(e, mTempdiffsub, gridTempdiffsub, weightTemp); //The numerical diffusion is interpolated back to the grid
                } }}
    }
            
    #pragma omp parallel for schedule(static) private(s, t, st)
    for (s=0;s<Sbase;s++){ for (t=0;t<Tbase;t++)	{
        st=IN0(s,t);
        if (weightTemp[st]>0) gridTempdiffsub[st]/=weightTemp[st];  //it is not a big problem if a node is not assigned a weight. the correction will then just be zero
        gridTempdiffrem[st]=gridTemp[st]-gridoldTemp[st]-gridTempdiffsub[st]; //Numerical is subtracted from the temperature correction
    } }
    #pragma omp parallel for schedule(static) private(e)
    for (e=0;e<N;e++) {
            mTemp[e]+=rhogridtomarker(e,gridTempdiffrem); //final correction
    }
}

void outputtempprofiles() {
    FILE *file0, *file1, *file2, *file3, *file4;
    char *writemode;
    double *gridoldTemp, horav, centertemp;
    gridoldTemp=cx2; //assumes that updatetemp has just been called so that cx2 still holds the old temperature field
    if (iterationnumber==0) writemode="wb"; else writemode="ab";
    
    file0 = fopen("outfiles/times_temp.txt", writemode); //redundant!
    file1 = fopen("outfiles/zhist.txt", writemode);
    file2 = fopen("outfiles/avtemphist.txt", writemode);
    file3 = fopen("outfiles/centertemphist.txt", writemode);
    file4 = fopen("outfiles/meltvol.txt",writemode);
    
    fwrite(&time, sizeof(double), 1, file0);
    fwrite(z, sizeof(double), Tbase, file1);
    fwrite(&meltproductionVOL, sizeof(double), 1, file4);

    for (t=0;t<Tbase;t++) {
        horav=0;
        for (s=0;s<Sbase;s++) {
            st=IN0(s, t);
            horav+=gridoldTemp[st];
        }
        horav/=Sbase; centertemp=gridoldTemp[IN0((Sbase-1)/2, t)];
        fwrite(&horav, sizeof(double), 1, file2);
        fwrite(&centertemp, sizeof(double), 1, file3);
    }
    
    
    fclose(file0); fclose(file1); fclose(file2); fclose(file3); fclose(file4);
}

