function [mTemp,mtype,surface,midcrust,moho,topast]=steadystatelithos1(mx,mz,xgrid);
S=201; T=201; dim=[S T];
xmax=max(xgrid); zmax=max(mz); N=length(mx);
dx=xmax/(S-1); dz=zmax/(T-1); x=0:dx:xmax; z=0:dz:zmax;
transitionwidth=20e3;
rhouc=2800; rholc=3000; rhom=3300; alpha=0.5; %the fraction of upper to lower crust
thinliththickness=180e3; thickliththickness=0; meanthickness=0.5*(thinliththickness+thickliththickness);
liththickness=thinliththickness-thickliththickness*exp(-(x-xmax/2).^2./(2*transitionwidth^2));

thincrustthickness=40e3; thickcrustthickness=5e3; meanthickness=0.5*(thincrustthickness+thickcrustthickness);
crustthickness=thincrustthickness+thickcrustthickness*exp(-(x-xmax/2).^2./(2*transitionwidth^2));
% crustthickness=meanthickness+(meanthickness-thincrustthickness)*erf(1/transitionwidth*(x-xmax/2));
% crustthickness=thincrustthickness+thickcrustthickness*exp(-(x-xmax/2).^2./(2*transitionwidth^2))-2*thickcrustthickness*exp(-(x-xmax/2).^2./(2/10*transitionwidth^2));

Tsurface=0;

C1=(alpha*rhouc+(1-alpha)*rholc)*thincrustthickness;
C2=(alpha*rhouc+(1-alpha)*rholc)*crustthickness;
uplift=(C1-C2+(crustthickness-thincrustthickness)*rhom)/rhom;
weaklayerthickness=40000-uplift;

%midcrustdepth=0.5*alpha*crustthickness.*(1-erf(1/transitionwidth*(x-xmax/2)))+weaklayerthickness;
midcrustdepth=alpha*crustthickness+weaklayerthickness;
mohodepth=crustthickness+weaklayerthickness;
lithdepth=liththickness+weaklayerthickness;


[Z,X] = meshgrid(z,x);

k=3*ones(S,T); H=zeros(S,T);
for s=1:S
    for t=1:T
        if z(t)<midcrustdepth(s) & z(t)>=weaklayerthickness(s)
            H(s,t)=1.25e-6;
            k(s,t)=2.5;
        end
        if z(t)<mohodepth(s) & z(t)>=midcrustdepth(s)
            H(s,t)=1.25e-6;
            k(s,t)=2.5;
        end
    end
end
        
boundary=nan(S,T);
for s=1:S
    Ta=adiabat_empirical(z,1300,40e3);
    for t=1:T
        if s==1 | s==S
            boundary(s,t)=inf;
        end
        if t==1 & ~(s==1 || s==S)
            boundary(s,t)=Tsurface;
        end
        if z(t)>lithdepth(s) & ~(s==1 || s==S)
            %boundary(s,t)=(z(t)-lithdepth(s)-weaklayerthickness(s))*0.25e-3+1300;
            %boundary(s,t)=(z(t)-weaklayerthickness(s))*0.63332e-3+1315; %husk!!
            boundary(s,t)=Ta(t);
        end        
        if z(t)<weaklayerthickness(s) & ~(s==1 || s==S)
            boundary(s,t)=Tsurface;
        end
    end
end

R=zeros(S*T,1);

rowind=zeros(5*S*T,1);
colind=zeros(5*S*T,1);
vals=zeros(5*S*T,1);
fillcount=1;

for s=1:S
    for t=1:T
        if isfinite(boundary(s,t))
            rowind(fillcount)=sub2ind(dim,s,t);
            colind(fillcount)=sub2ind(dim,s,t);
            vals(fillcount)=1;
            R(sub2ind(dim,s,t))=boundary(s,t);
            fillcount=fillcount+1;
        end
        if isinf(boundary(s,t))
            rowind(fillcount)=sub2ind(dim,s,t);
            colind(fillcount)=sub2ind(dim,s,t);
            vals(fillcount)=1;
            fillcount=fillcount+1;

            rowind(fillcount)=sub2ind(dim,s,t);
            if (s>1)colind(fillcount)=sub2ind(dim,s-1,t); else colind(fillcount)=sub2ind(dim,s+1,t); end
            vals(fillcount)=-1;
            fillcount=fillcount+1;

            R(sub2ind(dim,s,t))=0;
        end

        if isnan(boundary(s,t))
            a1=-0.5*(k(s+1,t)+k(s,t))/dx^2;
            a2=-0.5*(k(s-1,t)+k(s,t))/dx^2;
            a3=-0.5*(k(s,t+1)+k(s,t))/dz^2;
            a4=-0.5*(k(s,t-1)+k(s,t))/dz^2;

            rowind(fillcount)=sub2ind(dim,s,t);
            colind(fillcount)=sub2ind(dim,s,t);
            vals(fillcount)=-(a1+a2+a3+a4);
            fillcount=fillcount+1;

            rowind(fillcount)=sub2ind(dim,s,t);
            colind(fillcount)=sub2ind(dim,s+1,t);
            vals(fillcount)=a1;
            fillcount=fillcount+1;

            rowind(fillcount)=sub2ind(dim,s,t);
            colind(fillcount)=sub2ind(dim,s-1,t);
            vals(fillcount)=a2;
            fillcount=fillcount+1;

            rowind(fillcount)=sub2ind(dim,s,t);
            colind(fillcount)=sub2ind(dim,s,t+1);
            vals(fillcount)=a3;
            fillcount=fillcount+1;

            rowind(fillcount)=sub2ind(dim,s,t);
            colind(fillcount)=sub2ind(dim,s,t-1);
            vals(fillcount)=a4;
            fillcount=fillcount+1;

            R(sub2ind(dim,s,t))=H(s,t);
        end
    end
end
ind=find(rowind>0 & colind>0);
L=sparse(rowind(ind),colind(ind),vals(ind),S*T,S*T);

solution=L\R;
temp=zeros(S,T);
for s=1:S for t=1:T
        temp(s,t)=solution(sub2ind(dim,s,t));
    end
end

% for s=2:floor(S/2)
%     temp(s,:)=temp(1,:);
% end
% for s=floor(S/2)+1:S-1
%     temp(s,:)=temp(S,:);
% end
    

T1=interp2(X',Z',temp',x,weaklayerthickness+dz); T2=interp2(X',Z',temp',x,weaklayerthickness+2*dz); ki=interp2(X',Z',k',x,weaklayerthickness+1.5*dz);
surfaceheatflow=ki.*(T2-T1)/dz;

mTemp=zeros(1,N);
mtype=zeros(1,N)+3; %3 is the mantle
Ndivisions=16
for i=0:Ndivisions-1
    index=floor(i*N/Ndivisions)+1:floor((i+1)*N/Ndivisions);
    [min(index) max(index)]
    
    
    mTemp(index)=interp2(X',Z',temp',mx(index),mz(index));




    Iweaklayer=find(interp1(x,weaklayerthickness,mx(index))>mz(index)); mtype(index(Iweaklayer))=0; %0 is the weak layer
    Iucrust=find(interp1(x,weaklayerthickness,mx(index))<=mz(index) & interp1(x,midcrustdepth,mx(index))>mz(index)); mtype(index(Iucrust))=1; %1 is the upper crust
    Ilcrust=find(interp1(x,midcrustdepth,mx(index))<=mz(index) & interp1(x,mohodepth,mx(index))>mz(index)); mtype(index(Ilcrust))=2; %2 is the lower crust
    Iast=find(interp1(x,lithdepth,mx(index))<mz(index)); mtype(index(Iast))=4; %4 is the asthenosphere

    salthickness=1e3; preriftthickness=1e3;
    %%Iprerift=find(mz(index)>interp1(x,weaklayerthickness,mx(index)) & mz(index)<=interp1(x,weaklayerthickness+preriftthickness,mx(index))); mtype(index(Iprerift))=5;
    %%Isalt=find(mz(index)>interp1(x,weaklayerthickness+preriftthickness,mx(index)) & mz(index)<=interp1(x,weaklayerthickness+preriftthickness+salthickness,mx(index))); mtype(index(Isalt))=6;
end

surface=interp1(x,weaklayerthickness,xgrid);
midcrust=interp1(x,midcrustdepth,xgrid);
moho=interp1(x,mohodepth,xgrid);
topast=interp1(x,lithdepth,xgrid);

figure(1); contour(X,Z,temp,10);  set(gca,'yDir','rev'); shading flat; hold on; plot(x,lithdepth,'k','linewidth',3); plot(x,mohodepth,'k','linewidth',3); plot(x,midcrustdepth,'k','linewidth',3); plot(x,weaklayerthickness,'k','linewidth',3); colorbar;
figure(2); plot(x,surfaceheatflow);
end
