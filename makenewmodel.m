clear all; close all;
S=257; T=257;
%S=2049; T=609;
dtnom=0.5*10^5*60*60*24*365.25;
xmax=2000e3; zmax=600e3; %med ekstra bufferzone som ikke er en del af modellen
[dz z]=customgrid(0,zmax,[1 1 1],[20e3 40e3],[10e3 10e3],T);
[dx x]=customgrid(0,xmax,[1 1 1],[50e3 250e3],[50e3 50e3],S);

maxlength0=sum(dx); dx=xmax/maxlength0*dx;
maxlength0=sum(dz); dz=zmax/maxlength0*dz;

m_pr_cell=12*ones(S-1,T-1); minmarkersincell=min(m_pr_cell(:)); if minmarkersincell<4 display('Number of markers is too low - add more'); end
N=sum(m_pr_cell(:));
mx=zeros(1,N); mz=zeros(1,N); mcs=zeros(1,N); mct=zeros(1,N);

mtype=zeros(1,N);
time=0;
multigridlevels=1;

% multigridlevels=4;
bound=[1 1 1 1];

x=zeros(1,S); z=zeros(1,T); xvz=zeros(1,S)+0.5*dx(1); zvx=zeros(1,T)+0.5*dz(1);
for s=1:S-1 x(s+1)=x(s)+dx(s); end
xvx=x;
for s=1:S-2 xvz(s+1)=xvz(s)+0.5*(dx(s)+dx(s+1)); end
xvz(S)=xvz(S-1)+dx(S-1);
for t=1:T-1 z(t+1)=z(t)+dz(t); end 
zvz=z; 
for t=1:T-2 zvx(t+1)=zvx(t)+0.5*(dz(t)+dz(t+1)); end
zvx(T)=zvx(T-1)+dz(T-1);
xP=xvz; zP=zvx;

e=1; %rand=0.5;
for s=1:S-1
    for t=1:T-1
        mx(e)=x(s)+0.5*dx(s)*rand;
        mz(e)=z(t)+0.5*dz(t)*rand;
        mcs(e)=s-1;
        mct(e)=t-1;
        mtype(e)=(-1)^s*(-1)^t;
        e=e+1;

        mx(e)=x(s)+0.5*dx(s)+0.5*dx(s)*rand;
        mz(e)=z(t)+0.5*dz(t)*rand;
        mcs(e)=s-1;
        mct(e)=t-1;
        mtype(e)=(-1)^s*(-1)^t;
        e=e+1;

        mx(e)=x(s)+0.5*dx(s)*rand;
        mz(e)=z(t)+0.5*dz(t)+0.5*dz(t)*rand;
        mcs(e)=s-1;
        mct(e)=t-1;
        mtype(e)=(-1)^s*(-1)^t;
        e=e+1;
        
        mx(e)=x(s)+0.5*dx(s)+0.5*dx(s)*rand;
        mz(e)=z(t)+0.5*dz(t)+0.5*dz(t)*rand;
        mcs(e)=s-1;
        mct(e)=t-1;
        mtype(e)=(-1)^s*(-1)^t;
        e=e+1;
        for m=5:m_pr_cell(s,t)
            mx(e)=x(s)+dx(s)*rand;
            mz(e)=z(t)+dz(t)*rand;
            mcs(e)=s-1;
            mct(e)=t-1;
            mtype(e)=(-1)^s*(-1)^t;
            e=e+1;
        end
    end
end

clear mtype;
mfinitestrain=-1e-6*ones(1,N);


[mTemp,mtype,surface,midcrust,moho,topast]=steadystatelithos1(mx,mz,x);
% slopeangle=30;
% slope=tan(slopeangle*pi/180);
% ind=find(mtype==3 & abs(mz-moho((S-1)/2)-slope*(mx-xmax/2))<10e3/cos(slopeangle*pi/180));
% mtype(ind)=6;

Nsurfaces=4; surfaces=zeros(S,Nsurfaces);
surfaces(:,4)=surface;
surfaces(:,3)=midcrust;
surfaces(:,2)=moho;
surfaces(:,1)=topast;

fid = fopen('0.txt','wb');
fwrite(fid,N,'int32');
fwrite(fid,S,'int32');
fwrite(fid,T,'int32');
fwrite(fid,multigridlevels,'int32');
fwrite(fid,bound,'int32');
fwrite(fid,time,'real*8');
fwrite(fid,dtnom,'real*8');
fwrite(fid,dx,'real*8');
fwrite(fid,dz,'real*8');
fwrite(fid,mx,'real*8');
fwrite(fid,mz,'real*8');
fwrite(fid,mcs,'int32');
fwrite(fid,mct,'int32');
fwrite(fid,mtype,'int32');
fwrite(fid,mTemp,'real*8');

%clear('mx','mz','mtype','mTemp');
mstressxx=zeros(1,N); mstressxz=zeros(1,N);
mstrainxx=zeros(1,N); mstrainxz=zeros(1,N);
% mfinitestrain=-1e-6*ones(1,N);

fwrite(fid,mstressxx,'real*8');
fwrite(fid,mstressxz,'real*8');
fwrite(fid,mstrainxx,'real*8');
fwrite(fid,mstrainxz,'real*8');
fwrite(fid,mfinitestrain,'real*8');

clear('mstressxx','mstressxz','mstrainxx','mstrainxz','mfinitestrain');
vx=1e-11*randn(S,T); vz=1e-11*randn(S,T); P=1e-11*rand(S,T); %tilf�ldige startv�rdier

fwrite(fid,vx(:),'real*8');
fwrite(fid,vz(:),'real*8');
fwrite(fid,P(:),'real*8');
fwrite(fid,Nsurfaces,'int32');
for sur=1:Nsurfaces
    fwrite(fid,surfaces(:,sur),'real*8');
end
fclose(fid);


% I1=find(mtype==1);
% I2=find(mtype==-1);
% plot(mx(I1),mz(I1),'r.','MarkerSize',4); hold on;
% plot(mx(I2),mz(I2),'b.','MarkerSize',4); hold off;
% figure;
% plot(mx(:),mz(:),'r.'); hold on;
% [X,Z] = meshgrid(x,z);
% plot(X(:),Z(:),'g.'); hold off; axis ij;
figure;
scalardotplot(mx,mz,mTemp,[0 1850]);
%clear ('mTemp','mcs','mct','vx','vz','P');
figure;
scalardotplot(mx,mz,mtype,[0 6]); hold on;
% for s=1:S
%     plot([x(s) x(s)],[z(1) z(T)],'k','linewidth',0.1);
% end
% for t=1:T
%     plot([x(1) x(S)],[z(t) z(t)],'k','linewidth',0.1);
% end
% axis equal;

% figure(1); plot(x,surfaces(:,1),'k','linewidth',3); plot(x,surfaces(:,2),'k','linewidth',3); plot(x,surfaces(:,3),'k','linewidth',3);

% out=readtest; outint=griddata(mx,mz,out,X,Z);
% figure(2); pcolor(X,Z,vx'); figure(3); pcolor(X,Z,outint);

% [Xvx,Zvx] = meshgrid(xvx,zvx);
% out=readtest; outint=griddata(mx,mz,out,Xvx,Zvx);
% figure(2); pcolor(Xvx,Zvx,vx'); figure(3); pcolor(Xvx,Zvx,outint);
% figure(3); plot(Xvx(4,:),vx(4,:),'b.'); hold on; plot(Xvx(4,:),outint(:,4),'r.');
% figure(4); plot(Zvx(:,4),vx(:,4),'b.'); hold on; plot(Zvx(:,4),outint(4,:),'r.');

% [Xvz,Zvz] = meshgrid(xvz,zvz);
% out=readtest; outint=griddata(mx,mz,out,Xvz,Zvz);
% figure(2); pcolor(Xvz,Zvz,vz'); figure(3); pcolor(Xvz,Zvz,outint);
% figure(3); plot(Xvz(4,:),vz(4,:),'b.'); hold on; plot(Xvz(4,:),outint(:,4),'r.');
% figure(4); plot(Zvz(:,4),vz(:,4),'b.'); hold on; plot(Zvz(:,4),outint(4,:),'r.');

% out=readtest; range=[0 1]; mtressgrid=griddata(mx,mz,mstressxx,X,Z);
% figure(2); pcolor(X,Z,mtressgrid); figure(3); pcolor(X,Z,out');

% [XP,ZP] = meshgrid(xP,zP);
% out=readtest; mtressgrid=griddata(mx,mz,mstressxx,XP,ZP);
% figure(2); pcolor(XP,ZP,mtressgrid); figure(3); pcolor(XP,ZP,out');
% figure(4); plot(XP(4,:),mtressgrid(4,:),'b.'); hold on; plot(XP(4,:),out(:,4),'r.');
% figure(5); plot(ZP(:,4),mtressgrid(:,4),'b.'); hold on; plot(ZP(:,4),out(4,:),'r.');

% [XP,ZP] = meshgrid(xP,zP);
% out=readtest; mtressgrid=griddata(mx,mz,mstressxx,XP,ZP); outgrid=griddata(mx,mz,out,XP,ZP);
% figure(2); pcolor(XP,ZP,mtressgrid); figure(3); pcolor(XP,ZP,outgrid);
% figure(4); plot(XP(4,:),mtressgrid(4,:),'b.'); hold on; plot(XP(4,:),outgrid(4,:),'r.');
% figure(5); plot(ZP(:,4),mtressgrid(:,4),'b.'); hold on; plot(ZP(:,4),outgrid(:,4),'r.');

% [X,Z] = meshgrid(x,z);
% out=readtest; mtressgrid=griddata(mx,mz,mstressxx,X,Z); outgrid=griddata(mx,mz,out,X,Z);
% figure(2); pcolor(X,Z,mtressgrid); figure(3); pcolor(X,Z,outgrid);
% figure(3); plot(X(4,:),mtressgrid(4,:),'b.'); hold on; plot(X(4,:),outgrid(4,:),'r.');
% figure(4); plot(Z(:,4),mtressgrid(:,4),'b.'); hold on; plot(Z(:,4),outgrid(:,4),'r.');

% [X,Z] = meshgrid(x,z);
% outgrid=readtest; mtressgrid=griddata(mx,mz,mstressxz,X,Z);
% figure(2); pcolor(X,Z,mtressgrid); figure(3); pcolor(X,Z,outgrid');
% figure(4); plot(X(4,:),mtressgrid(4,:),'b.'); hold on; plot(X(4,:),outgrid(:,4),'r.');
% figure(5); plot(Z(:,4),mtressgrid(:,4),'b.'); hold on; plot(Z(:,4),outgrid(4,:),'r.');

% figure(3); plot(X(4,:),mtressgrid(4,:),'b.'); hold on; plot(X(4,:),out(:,4),'r.');
% figure(4); plot(Z(:,4),mtressgrid(:,4),'b.'); hold on;
% plot(Z(:,4),out(4,:),'r.');

% exit;
