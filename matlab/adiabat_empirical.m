function Tp=adiabat_empirical(z,Tp,top)
Tpemp=[1200 1250 1300 1350];
slopeemp=[0.4418 0.4648 0.4878 0.5108]*1e-3;
P=polyfit(Tpemp,slopeemp,1);

Tp=Tp+(z-top)*polyval(P,Tp);

% Tpemp=[1200 1250 1300 1350]; Tbottomemp=[1408 1469 1530 1590]; P=polyfit(Tpemp,Tbottomemp,1); polyval(P,1400)
    
    
% g=9.82;
% T0=1400;
% N=length(z);
% T=zeros(1,N);
% T(1)=T0;
% 
% T=T+273;
% 
% for t=1:N-1
%     Cu=(238.64-20.013e2*T(t)^-0.5-11.642e7*T(t)^-3)*7.1077;
%     alphau=2.77e-5+0.97e-8*T(t)-0.32*T(t)^-2;
%     %Cu=1e3; alphau=3e-5;
%     T1=T(t)+T(t)*alphau*g/Cu*(z(t+1)-z(t));
%     Cd=(238.64-20.013e2*T1^-0.5-11.642e7*T1^-3)*7.1077;
%     alphad=2.77e-5+0.97e-8*T1-0.32*T1^-2;
%     %Cd=1e3; alphad=3e-5;
%     T2=T(t)+T(t)*alphad*g/Cd*(z(t+1)-z(t));
%     if (z(t)>top)
%         T(t+1)=0.5*(T1+T2);
%     else
%         T(t+1)=T(t);
%     end
% end
% T=T-273;