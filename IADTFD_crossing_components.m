
close all;
clear all;
addpath 'D:\tfsa_5-5\windows\win64_bin\'
t=0:1/128:1-1/128;
f=0:0.5:64-0.5;
M=128;
n=0:M-1;

Sig1 =1*exp(1i*(2*pi*(5*t +15*t.^3)));
Sig2 =1*exp(1i*(2*pi*(12*t +15*t.^3)));
Sig3 =1*exp(1i*(2*pi*(50*t -1*13*t.^3)));
Sig4=exp(1i*(2*pi*(60*t )));
IF_O(:,1)=45*t.^2+12;
IF_O(:,1)=45*t.^2+5;
IF_O(:,2)=-39*t.^2+50;
%IF_O(:,3)=60;
SigA=1*Sig1+0*Sig2+1*Sig3+0*Sig4;
%x=cos(0.0005*2*pi*n.^3/M+2*2*pi*n/M)+cos(0.0005*2*pi*n.^3/M+30*2*pi*n/M);%+exp(1i*0.001*2*pi*n.^3/M+1i*15*2*pi*n/M);
%x=cos(0.001*2*pi*n.^3/M+2*2*pi*n/M)+cos(-0.001*2*pi*n.^3/M+45*2*pi*n/M);%+exp(1i*0.001*2*pi*n.^3/M+1i*15*2*pi*n/M);
%x=exp(1i*0.0005*2*pi*n.^3/M+1i*2*2*pi*n/M)+exp(1i*0.0005*2*pi*n.^3/M+1i*10*2*pi*n/M)+exp(1i*52*2*pi*n/M);
%x=exp(1i*0.001*2*pi*n.^3/M+1i*2*2*pi*n/M)+exp(-1i*0.001*2*pi*n.^3/M+1i*45*2*pi*n/M);

x=exp(1i*0.0004*2*pi*n.^3/M+1i*2*2*pi*n/M)+exp(1i*0.0004*2*pi*n.^3/M+1i*20*2*pi*n/M)+1*exp(-1i*0.0012*2*pi*n.^3/M+1i*60*2*pi*n/M) +0*exp(1i*25*2*pi*n/M);
x=exp(1i*0.0012*2*pi*n.^3/M+1i*2*2*pi*n/M)+0*exp(1i*0.0008*2*pi*n.^3/M+1i*20*2*pi*n/M)+1*exp(-1i*0.0012*2*pi*n.^3/M+1i*60*2*pi*n/M) +0*exp(1i*25*2*pi*n/M);


x=SigA;

N=length(x);
iii= randsample(128,32) ;
iii=[5:15 25:30 50:60  75:80  110:120];
Ispec=quadtfd(x,length(x)-1,1,'specx',45,'hamm',128);
Ispec=Ispec/sum(abs(Ispec(:)));
Imb=quadtfd(x,length(x)/4-1,1,'mb',0.1,128);

Imb(Imb<0)=0;
Imb=Imb/sum(abs(Imb(:)));Imb(Imb<0)=0;
Imb=Imb/sum(abs(Imb(:)));
%x(iii)=0;
Fre_Grid=1*length(x);
Fre_Line  = 2*pi*(1:Fre_Grid)/Fre_Grid;
%I1=min(HTFD_new1(x,2,30,64/1),HTFD_new1(x,2,50,84/1));
I1=HTFD_new1(x,2,30,64/1);
I1=I1/sum(abs(I1(:)));
I=HTFD_new1(x,2,20,48/1);
for i=1:4
    I= post_processing_directional(I,2,20,48);
    if i==1
        II=I;
        II(II<0)=0;
        II=II/sum(abs(II(:)));
    end
end
I(I<0)=0;
II(II<0)=0;
I=I/sum(abs(I(:)));
t=1:1:128;
f=0:1/256:0.5-1/256;
figure; imagesc(t,f,I1)

set(gcf,'Position',[20 100 640 500]);
xlabel('Time / Sec','FontSize',20,'FontName','Times New Roman');
ylabel('Frequency / Hz','FontSize',20,'FontName','Times New Roman');
set(gca,'YDir','normal');
title('(a)','FontSize',24,'FontName','Times New Roman');
set(gca,'FontSize',20);

figure; imagesc(t,f,II)

set(gcf,'Position',[20 100 640 500]);
xlabel('Time / Sec','FontSize',20,'FontName','Times New Roman');
ylabel('Frequency / Hz','FontSize',20,'FontName','Times New Roman');
set(gca,'YDir','normal');
title('(b)','FontSize',24,'FontName','Times New Roman');
set(gca,'FontSize',20);

figure; imagesc(t,f,I)

set(gcf,'Position',[20 100 640 500]);
xlabel('Time / Sec','FontSize',20,'FontName','Times New Roman');
ylabel('Frequency / Hz','FontSize',20,'FontName','Times New Roman');
set(gca,'YDir','normal');
title('(c)','FontSize',24,'FontName','Times New Roman');
set(gca,'FontSize',20);

figure; imagesc(t,f,Ispec)

set(gcf,'Position',[20 100 640 500]);
xlabel('Time / Sec','FontSize',20,'FontName','Times New Roman');
ylabel('Frequency / Hz','FontSize',20,'FontName','Times New Roman');
set(gca,'YDir','normal');
title('(d)','FontSize',24,'FontName','Times New Roman');
set(gca,'FontSize',20);

figure; imagesc(t,f,Imb)

set(gcf,'Position',[20 100 640 500]);
xlabel('Time / Sec','FontSize',20,'FontName','Times New Roman');
ylabel('Frequency / Hz','FontSize',20,'FontName','Times New Roman');
set(gca,'YDir','normal');
title('(e)','FontSize',24,'FontName','Times New Roman');
set(gca,'FontSize',20);

sum(sum(abs(Imb.^0.5))).^2
sum(sum(abs(Ispec.^0.5))).^2
sum(sum(abs(I.^0.5))).^2
sum(sum(abs(II.^0.5))).^2
sum(sum(abs(I1.^0.5))).^2

