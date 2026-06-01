addpath('~/code/Code_Kenni/tools/');

% plot crustal thickness against time 
clear all;
close all;

home = pwd;

dirname = {'./Ice/2D_Tp1400_riftweak_noice/'...
    './Ice/2D_Tp1400_riftweak_ice/'...
    './Ice/2D_Tp1400_riftweak_step_dt1e2/'...
    './Ice/more_output/'};

name = {'1400 no ice'...
    '1400 ice'...
    '1400 step ice'...
    'more'};

smallprint = [200 200 200 100];

S=readdata_name(dirname{1},0,2);

yr = 60*60*24*365;
my = 1e-6;
km = 1e-3;

fullsp = 10*ones(size(dirname));

windowSize = 10;

fig1 = figure(101);
set(fig1,'position',[10 10 800 800],'PaperPositionMode','auto');
fig2 = figure(102);
set(fig2,'position',[10 10 800 800],'PaperPositionMode','auto');

cm = colormap('jet');

for n = 1:length(dirname)
    
    hlfsp = .5e-3*fullsp(n);
    
    fid=fopen([dirname{n} 'outfiles/meltvol.txt'],'rb');
    melt1=fread(fid,'real*8');
    fclose(fid);
    
    fid=fopen([dirname{n} 'outfiles/bulk_comp.txt'],'rb');
    bulk1=fread(fid,'real*8');
    fclose(fid);
    
    fid=fopen([dirname{n} 'outfiles/times.txt'],'rb');
    times1=fread(fid,'real*8');
    fclose(fid);
    
    times=times1(1:end-1);
    Ntimes=length(times);
    
    fid=fopen([dirname{n} 'outfiles/surface.txt'],'rb');
    surface1=fread(fid,'real*8');
    surface1=surface1(1:S*Ntimes);
    fclose(fid);
    surface1=reshape(surface1,S,Ntimes);
    
    fid=fopen([dirname{n} 'outfiles/basement.txt'],'rb');
    basement1=fread(fid,'real*8');
    basement1=basement1(1:S*Ntimes);
    fclose(fid);
    basement1=reshape(basement1,S,Ntimes);
    
    fid=fopen([dirname{n} 'outfiles/moho.txt'],'rb');
    moho1=fread(fid,'real*8');
    moho1=moho1(1:S*Ntimes);
    fclose(fid);
    moho1=reshape(moho1,S,Ntimes);
    
    cthick1 = zeros(size(times));
    for i = 1:length(times)
        cthick1(i) = min(abs(moho1(:,i)-basement1(:,i)));
    end
    stretch1 = cthick1(1)./cthick1;
    
    cd(dirname{n});
    !grep dtnom status.txt | awk '{print $4}' > dt.txt
    dtall = load('dt.txt','-ascii');
    !rm dt.txt
    cd(home);
    dt = zeros(size(times1));
    dt(1) = dtall(1)*1e6*yr;
    dt(2:length(dtall(smallprint(n)-1:smallprint(n):end))+1) = dtall(smallprint(n)-1:smallprint(n):end)*1e6*yr;
    
    figure(fig1)
    subplot(2,2,1)
    plot(my*times1./yr,.5*km*yr*melt1./(dt.*hlfsp),'color',cm(round(length(cm)/length(dirname)*n),:));
    set(gca,'xlim',[0 40])
    hold on
        
    subplot(2,2,2)
    % plot(my*times1(2:end)./yr,1e-3*cthick1,'color',cm(round(length(cm)/length(dirname)*n),:));
    plot(my*times1./yr,.5*km*yr*melt1,'color',cm(round(length(cm)/length(dirname)*n),:));
    set(gca,'xlim',[0 40])
    hold on
    
    subplot(2,2,3)
    semilogx(stretch1,.5*km*yr*filter(ones(1,windowSize)/windowSize,1,melt1(2:end))./(dt(2:end).*hlfsp),'color',cm(round(length(cm)/length(dirname)*n),:));
    set(gca,'xlim',[1 1e4])
    hold on
    
    figure(fig2)
    subplot(2,1,1)
    plot(my*times1./yr,.5*km*yr*melt1./(dt.*hlfsp),'color',cm(round(length(cm)/length(dirname)*n),:));
    set(gca,'xlim',[5 22])
    hold on
    
    if n == length(dirname)
        cd(dirname{n});
        !grep deglaciate status.txt | awk '{print $12}' > deglaciate.txt
        deglaciate = load('deglaciate.txt','-ascii');
        !rm deglaciate.txt
        cd(home);
        stepice = zeros(size(times1));
        stepice(1) = 2e3*deglaciate(1);
        stepice(2:end) = 2e3*deglaciate(smallprint(n):smallprint(n):end);
        figure(fig2)
        subplot(2,1,2)
        plot(my*times1./yr,stepice,'k');%,'color',cm(round(length(cm)/length(dirname)*n),:));
        set(gca,'xlim',[5 22])
        set(gca,'ylim',[0 5e3])
    end    
            
    % clearvars times1 melt1 cthick1 bulk1 stretch1;
    
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
xlabel('Extension rate (mmyr^{-1})');
ylabel('Rift duration (km)');

figure(fig2)
subplot(2,1,1)
xlabel('Time (Myr)');
ylabel('Melt thickness (km)');
legend(name','location','northeast')
grid on
subplot(2,1,2)
xlabel('Time (Myr)');
ylabel('Ice thickness (m)');

nameeps = './plots/melt-thickness.eps';
print(fig2,'-depsc',nameeps);

