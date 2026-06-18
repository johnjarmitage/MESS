function scalardotplot(mx,mz,mscalar,range)
N=length(mx);

cmap=colormap;
M=size(cmap);
M=M(1);

clf;
hold on;

Ncolors=100;
step=(range(2)-range(1))/Ncolors;
scalarintervals=range(1):step:range(2);

for i=1:Ncolors
    smin=scalarintervals(i);
    smax=scalarintervals(i+1);
    thiscolor=find(mscalar>=smin & mscalar<smax);
    if i==Ncolors thiscolor=find(mscalar>=smin & mscalar<=smax); end
    ndimscalar=(0.5*(smin+smax)-range(1))/(range(2)-range(1));
    if ndimscalar>1
        ndimscalar=1;
        display(sprintf(['Grænser overskredet - max temp: ' num2str(smax)]));
    end
    if ndimscalar<0;
        ndimscalar=0;
        display(sprintf(['Grænser overskredet - min temp: ' num2str(smin)]));
    end
    color(1)=interp1(0:1/(M-1):1,cmap(:,1),ndimscalar);
    color(2)=interp1(0:1/(M-1):1,cmap(:,2),ndimscalar);
    color(3)=interp1(0:1/(M-1):1,cmap(:,3),ndimscalar);
    plot(mx(thiscolor),mz(thiscolor),'.','color',color,'markersize',20);
end

caxis(range); colorbar;
axis ij;

end
