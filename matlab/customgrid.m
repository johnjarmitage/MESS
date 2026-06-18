function [dx x]=customgrid(xmin,xmax,amp,my,sigma,I)

xsample=xmin:(xmax-xmin)/1e5:xmax;
f=amp(1)+0.5*(erf((xsample-my(1))/(sigma(1)*2^0.5))+1)*(amp(2)-amp(1))+0.5*(erf((xsample-my(2))/(sigma(2)*2^0.5))+1)*(amp(3)-amp(2));

F=zeros(size(f));
for h=2:length(f)
    F(h)=F(h-1)+0.5*(f(h-1)+f(h))*(xsample(h)-xsample(h-1));
end
F=F/F(end);

x=zeros(1,I);
Fstep=1/(I-1);
x(1)=xmin; currentF=0;
for i=2:I-1
    currentF=(i-1)*Fstep;
    left=1; right=length(F);
    while right-left>1
        middle=floor(0.5*(right+left));
        if F(middle)<currentF
            left=middle;
        else
            right=middle;
        end
    end
    x(i)=xsample(left)+(currentF-F(left))*(xsample(right)-xsample(left))/(F(right)-F(left));
end
x(I)=xmax;

dx=diff(x);