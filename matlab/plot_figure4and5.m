addpath('~/code/FD/Code_Kenni/tools/');

% plot crustal thickness against time 
clear all;
close all;

%% indian nMORB, N. Atlantic nMORB, N. Atlantic eMORB (Klein, 2004) nMORB
%% (Kelemen et al., 2004)

LaN = [4.46     1.88    11.5    4.13];
CeN = [12.21    5.99    26      11.46];
PrN = [NaN      NaN     NaN     1.84];
NdN = [8.47     6.07    17.1    9.30];
SmN = [2.72     2.22    4.38    2.96];
EuN = [1.07     0.9     1.54    1.13];
GdN = [3.49     3.5     5.26    4.04];
TbN = [NaN      NaN     NaN     0.68];
DyN = [4.04     4.46    5.24    4.70];
HoN = [NaN      NaN     NaN     0.93];
ErN = [2.56     2.57    2.83    3.00];
TmN = [NaN      NaN     NaN     0.47];
YbN = [2.25     2.72    2.68    2.72];
LuN = [0.0354   NaN     NaN     0.41];

%% erta ale

LaE = [11.0262	8.00125	27.9279	7.7616	16.5669];
CeE = [25.7184	18.4353	57.6414	18.396	38.4468];
PrE = [];	
NdE = [14.4478	11.396	27.621	11.004	21.2382];
SmE = [3.4498	2.98775	5.4405	2.898	4.8369] ;
EuE = [1.2032	1.073	1.7577	1.0584	1.6353] ;
GdE = [3.6754	3.33925	4.8453	3.2676	4.7817] ;
TbE = [0.6486	0.5365	0.7254	0.546	0.7383] ;
DyE = [3.6754	3.626	4.1385	3.5616	4.5816] ;
HoE = [];
ErE = [2.0492	2.109	2.1204	2.0496	2.4771] ;
TmE = [];
YbE = [1.7672	1.8685	1.767	1.7892	2.1045] ;
LuE = [0.2726	0.2775	0.2697	0.2772	0.3243] ;
%crustE = [6 10];
crustE = [0 4]; % assuming stretch factor of 2

%% dabbahu

LaD = [11.3832353453	10.1363045046	10.6415699247	16.9314850358	16.516032461	11.4256152999	11.2509987712] ;
CeD = [26.4053503118	23.7483405674	24.895015317	37.585149517	36.6639156637	26.2583881585	26.233384358]  ;
PrD = [3.6443911019     3.2795129015	3.4160282497	4.8878506686	4.7424755469	3.5721582037	3.5767246542]  ;
NdD = [16.5383918973	14.9455316704	15.5886257472	20.8642454282	20.4319658819	16.1225978188	16.1070087624] ;
SmD = [4.2003646711     3.8739664245	4.0235962599	5.1459870433	4.9836737973	4.0454183185	4.0437977777]  ;
EuD = [1.5309280203     1.4407697352	1.5024131725	1.6167649432	1.5809151391	1.5216139265	1.5134557848]  ;
GdD = [4.5901371331     4.2117003811	4.3942564874	5.4761114442	5.3579764498	4.4553490297	4.4734200335]  ;
TbD = [0.7465774856     0.6920744181	0.7287166169	0.9262588888	0.8982326991	0.7314379349	0.7370750832]  ;
DyD = [4.6208889163     4.2683811145	4.4262238352	5.6441110175	5.5017773898	4.461709021     4.5009237165]  ;
HoD = [0.9268266042     0.8605892872	0.8846430702	1.1712304906	1.1287961413	0.8961500828	0.9147174234]  ;
ErD = [2.4363301576     2.230588513     2.3563699235	3.08515992      2.985761103     2.331080681     2.3704707107]  ;
TmD = [0.3338097082     0.314133063     0.3271393325	0.4410557064	0.4328597172	0.3276856623	0.3325569482]  ;
YbD = [2.0378558848     1.9311376458	1.9674758886	2.6835901534	2.6071560057	1.9994853172	2.0045363089]  ;
LuD = [0.3198956367     0.2847261397	0.2977928069	0.4214295886	0.4015829314	0.3010571469	0.300830097]   ;
%crustD = [11 15];
crustD = [8.5 11]; % assuming stretch factor of 2

%% henities bay

LaH = [9.7 5.5];
SmH = [3   2.1];
YbH = [2.0 1.5];
crustH = [20 25];

%% false bay

LaF = [16  11];
SmF = [4.9 3.7];
YbF = [3.3 2.6];
crustF = [6 8];

%% First get the variation in temperature

dirname = {'./Water100ppm/2D_Tp1350_riftweak/outfiles/'... 
    './Water100ppm/2D_Tp1400_riftweak/outfiles/'...
    './Water100ppm/2D_Tp1450_riftweak/outfiles/'...
    './Water100ppm/2D_Tp1350_riftstrong/outfiles/'...
    './Water100ppm/2D_Tp1400_riftstrong/outfiles/'...
    './Water100ppm/2D_Tp1450_riftstrong/outfiles/'};

name = {'1350{^o}C weak'... 
    '1400{^o}C weak'...
    '1450{^o}C weak'...
    '1350{^o}C strong'...
    '1400{^o}C strong'...
    '1450{^o}C strong'};


S=readdata_name(dirname{1},0,2);

yr = 60*60*24*365;
dt = 0.001*1e6*yr;
hlfsp = 7.5e-3;
my = 1e-6;
km = 1e-3;

windowSize = 10;

fig1 = figure(101);
set(fig1,'position',[10 10 800 800],'PaperPositionMode','auto');
fig3 = figure(103);
set(fig3,'position',[20 10 800 800],'PaperPositionMode','auto');

cm = colormap('jet');

for n = 1:length(dirname)
    
    fid=fopen([dirname{n} 'meltvol.txt'],'rb');
    melt1=fread(fid,'real*8');
    fclose(fid);
    
    fid=fopen([dirname{n} 'bulk_comp.txt'],'rb');
    bulk1=fread(fid,'real*8');
    fclose(fid);
    
    fid=fopen([dirname{n} 'times.txt'],'rb');
    times1=fread(fid,'real*8');
    fclose(fid);
    
    times=times1(1:end-1);
    Ntimes=length(times);
    
    fid=fopen([dirname{n} 'surface.txt'],'rb');
    surface1=fread(fid,'real*8');
    surface1=surface1(1:S*Ntimes);
    fclose(fid);
    surface1=reshape(surface1,S,Ntimes);
    
    fid=fopen([dirname{n} 'basement.txt'],'rb');
    basement1=fread(fid,'real*8');
    basement1=basement1(1:S*Ntimes);
    fclose(fid);
    basement1=reshape(basement1,S,Ntimes);
    
    fid=fopen([dirname{n} 'moho.txt'],'rb');
    moho1=fread(fid,'real*8');
    moho1=moho1(1:S*Ntimes);
    fclose(fid);
    moho1=reshape(moho1,S,Ntimes);
    
    cthick1 = zeros(size(times));
    for i = 1:length(times)
        cthick1(i) = min(abs(moho1(:,i)-basement1(:,i)));
    end
    stretch1 = cthick1(1)./cthick1;
        
    figure(fig1)
    subplot(2,2,1)
    plot(my*times1./yr,.5*km*yr*filter(ones(1,windowSize)/windowSize,1,melt1)./(dt*hlfsp),'color',cm(round(length(cm)/length(dirname)*n),:));
    set(gca,'xlim',[0 40])
    hold on
        
    subplot(2,2,2)
    plot(my*times1(1:end-1)./yr,1e-3*cthick1,'color',cm(round(length(cm)/length(dirname)*n),:));
    set(gca,'xlim',[0 40])
    hold on

    subplot(2,2,3)
    semilogx(stretch1,.5*km*yr*filter(ones(1,windowSize)/windowSize,1,melt1(1:end-1))./(dt*hlfsp),'color',cm(round(length(cm)/length(dirname)*n),:));
    set(gca,'xlim',[1 1e4])
    hold on
    
    figure(fig3)
    subplot(2,2,1)
    semilogx([1 1e4],[LaN(1)/YbN(1) LaN(1)/YbN(1)],'color',[.5 .5 .5])
    semilogx(stretch1,filter(ones(1,windowSize)/windowSize,1,bulk1(1:14:end-14)./bulk1(13:14:end-14)),'color',cm(round(length(cm)/length(dirname)*n),:));
    set(gca,'xlim',[1 1e4])
    set(gca,'ylim',[0 20])
    hold on
    
    subplot(2,2,2)
    semilogx([1 1e4],[SmN(1)/YbN(1) SmN(1)/YbN(1)],'color',[.5 .5 .5])
    semilogx(stretch1,filter(ones(1,windowSize)/windowSize,1,bulk1(5:14:end-14)./bulk1(13:14:end-14)),'color',cm(round(length(cm)/length(dirname)*n),:));
    set(gca,'xlim',[1 1e4])
    set(gca,'ylim',[0 6])
    hold on
    
    subplot(2,2,3)
    plot(.5*km*yr*filter(ones(1,windowSize)/windowSize,1,melt1(1:101))./(dt*hlfsp),...
        filter(ones(1,windowSize)/windowSize,1,bulk1(1:14:101*14)./bulk1(13:14:101*14+12)),'color',cm(round(length(cm)/length(dirname)*n),:))
    set(gca,'xlim',[0 25])
    set(gca,'ylim',[0 20])
    hold on

    subplot(2,2,4)
    plot(.5*km*yr*filter(ones(1,windowSize)/windowSize,1,melt1(1:101))./(dt*hlfsp),...
        filter(ones(1,windowSize)/windowSize,1,bulk1(5:14:101*14+4)./bulk1(13:14:101*14+12)),'color',cm(round(length(cm)/length(dirname)*n),:))
    set(gca,'xlim',[0 25])
    set(gca,'ylim',[0 6])
    hold on
           
    clearvars times1 melt1 cthick1 bulk1 stretch1;
    
end

dirname = {'./Extension_rate/Weak_010_1350/outfiles/'... 
    './Extension_rate/Weak_015_1350/outfiles/'...
    './Extension_rate/Weak_020_1350/outfiles/'...
    './Extension_rate/Weak_040_1350/outfiles/'... 
    './Extension_rate/Weak_060_1350/outfiles/'...
    './Extension_rate/Weak_100_1350/outfiles/'...
    './Extension_rate/Strong_010_1350/outfiles/'...
    './Extension_rate/Strong_015_1350/outfiles/'...
    './Extension_rate/Strong_020_1350/outfiles/'...
    './Extension_rate/Strong_040_1350/outfiles/'...
    './Extension_rate/Strong_060_1350/outfiles/'...
    './Extension_rate/Strong_100_1350/outfiles/'};

fullsp        = [10 15 20 40 60 100 10 15 20 40 60 100];
steadymelt    = zeros(length(dirname));
rift_duration = zeros(length(dirname));

for n = 1:length(dirname)
    
    fid=fopen([dirname{n} 'meltvol.txt'],'rb');
    melt1=fread(fid,'real*8');
    fclose(fid);
    
    fid=fopen([dirname{n} 'bulk_comp.txt'],'rb');
    bulk1=fread(fid,'real*8');
    fclose(fid);
    
    fid=fopen([dirname{n} 'times.txt'],'rb');
    times1=fread(fid,'real*8');
    fclose(fid);
    
    hlfsp = .5e-3*fullsp(n);
    melt2 = .5*km*yr*melt1./(dt*hlfsp);
    steadymelt(n) = mean(melt2(end-40:end));
    test = 0;
    m = 0;
    while (test == 0)
        m = m + 1;
        if (melt2(m) < 8)
            test = 0;
        else
            test = 1;
        end
    end
    rift_duration(n) = my*times1(m)./yr;
    
end
    

figure(fig1)
subplot(2,2,1)
xlabel('Time (Myr)');
ylabel('Melt thickness (km)');
subplot(2,2,2)
xlabel('Time (Myr)');
ylabel('Crustal thickness (km)');
legend(name','location','northeast')

subplot(2,2,3)
xlabel('Stretch factor');
ylabel('Melt thickness (km)'); 

subplot(2,2,4)
plot(fullsp(1:6),rift_duration(1:6),'ok')
hold on
plot(fullsp(7:12),rift_duration(7:12),'or')
legend('Weak','Strong','location','northeast')
plot(fullsp(1:6),rift_duration(1:6),'k')
plot(fullsp(7:12),rift_duration(7:12),'r')
xlabel('Extension rate (mmyr^{-1})');
ylabel('Rift duration (My)');
legend('Weak','Strong','location','northeast')

figure(fig3)
subplot(2,2,1)
semilogx([1 1e4],[mean(LaE./YbE) mean(LaE./YbE)],'k--');
hold on
semilogx([1 1e4],[mean(LaD./YbD) mean(LaD./YbD)],'k');
xlabel('Stretch factor');
ylabel('La/Yb');

subplot(2,2,2)
semilogx([1 1e4],[mean(SmE./YbE), mean(SmE./YbE)],'k--');
hold on
semilogx([1 1e4],[mean(SmD./YbD), mean(SmD./YbD)],'k');
xlabel('Stretch factor');
ylabel('Sm/Yb');

subplot(2,2,3)
errorbarxy(.5,mean(LaE./YbE),.5,.5,std(LaE./YbE),std(LaE./YbE))
hold on
errorbarxy(4.5,mean(LaD./YbD),4.5,4.5,std(LaD./YbD),std(LaD./YbD))
xlabel('Melt thickness (km)');
ylabel('La/Yb');
legend(name','location','northeast')

subplot(2,2,4)
errorbarxy(.5,mean(SmE./YbE),.5,.5,std(SmE./YbE),std(SmE./YbE))
hold on
errorbarxy(4.5,mean(SmD./YbD),4.5,4.5,std(SmD./YbD),std(SmD./YbD))
xlabel('Melt thickness (km)');
ylabel('Sm/Yb');





nameeps = './plots/figure04.eps';
print(fig1,'-depsc',nameeps);
nameeps = './plots/figure05.eps';
print(fig3,'-depsc',nameeps);
