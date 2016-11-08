close all
x1=0:50;
y1=besselj(0,x1);
x2=2.35:.001:2.45;
y2=besselj(0,x2);
fig1=figure(1);
plot(x1,y1,'b','linewidth',2)
hold on
plot(x1,0*x1,':k')
set(gca,'fontsize',15)
title ('bessel function')
xlabel('X')
ylabel('Y')
fig2=figure(2);
plot(x2,y2,'r')
hold on
plot(x2,0*x2,':k')
title ('close-up')
[h_m h_i]=inset(fig1,fig2);

set(h_i,'xtick',2.35:.025:2.45,'xlim',[2.35,2.45])