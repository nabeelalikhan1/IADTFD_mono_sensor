
close all;
clear all;
addpath 'D:\tfsa_5-5\windows\win64_bin\'
SampFreq=128;

M=128;
n=0:M-1;
t = 0:1/SampFreq:1-1/SampFreq;
%x=cos(0.0005*2*pi*n.^3/M+2*2*pi*n/M)+cos(0.0005*2*pi*n.^3/M+30*2*pi*n/M);%+exp(1i*0.001*2*pi*n.^3/M+1i*15*2*pi*n/M);
%x=cos(0.001*2*pi*n.^3/M+2*2*pi*n/M)+cos(-0.001*2*pi*n.^3/M+45*2*pi*n/M);%+exp(1i*0.001*2*pi*n.^3/M+1i*15*2*pi*n/M);
%x=exp(1i*0.0005*2*pi*n.^3/M+1i*2*2*pi*n/M)+exp(1i*0.0005*2*pi*n.^3/M+1i*10*2*pi*n/M)+exp(1i*52*2*pi*n/M);
%x=exp(1i*0.001*2*pi*n.^3/M+1i*2*2*pi*n/M)+exp(-1i*0.001*2*pi*n.^3/M+1i*45*2*pi*n/M);
%n=-M/2:M/2-1;
%x=exp(1i*0.001*2*pi*n.^3/M+1i*2*2*pi*n/M)+exp(-1i*0.001*2*pi*n.^3/M+1i*45*2*pi*n/M);% +exp(1i*55*2*pi*n/M);
x=1*exp(2*pi*1i*(0.35*n.^3/(3*M*M)+0.006*n.^2/(2*M)+0.08*n))+1*exp(2*pi*1i*(0.35*n.^3/(3*M*M)+0.002*n.^2/(2*M)+0.05*n));
%x=exp(1i*0.001*2*pi*n.^3/M+1i*4*pi*n/M)+0*exp(1i*0.0007*2*pi*n.^3/M+1i*8*2*pi*n/M)+1*exp(-1i*0.001*2*pi*n.^3/M+1i*0.95*pi*n);
x=1*exp(2*pi*1i*(0.35*n.^3/(3*M*M)+0.0825*n))+1*exp(2*pi*1i*(0.35*n.^3/(3*M*M)+0.05*n));

%x=1*exp(2*pi*1i*(0.5*n.^3/(3*M*M)+0.11*n))+1*exp(2*pi*1i*(0.5*n.^3/(3*M*M)+0.05*n));
%x=exp(1i*0.001*2*pi*n.^3/M+1i*2*2*pi*n/M)+exp(-1i*0.001*2*pi*n.^3/M+1i*45*2*pi*n/M);
%  Sig1 =1*exp(1i*(2*pi*(5*t +15*t.^3)));
%     Sig2 =1*exp(1i*(2*pi*(10*t +15*t.^3)));
%     Sig3 =1*exp(1i*(2*pi*(50*t -1*13*t.^3)));
%     Sig4=exp(1i*(2*pi*(60*t )));
%     IF_O(:,1)=45*t.^2+9;
%     IF_O(:,2)=45*t.^2+5;
%     %IF_O(:,3)=60;
%     SigA=1*Sig1+1*Sig2+0*Sig3+0*Sig4;
%
%
% x=SigA;


N=length(x);
iii= randsample(128,32) ;
Ispec=quadtfd(x,length(x)-1,1,'specx',25,'hamm',128);
Ispec=Ispec/sum(abs(Ispec(:)));
Imb=quadtfd(x,length(x)/4-1,1,'mb',0.1,128);

Imb(Imb<0)=0;
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
%for iii=1:5
%   Inew= post_processing_directional(adtfd,2,20,64);
%end
%I=adtfd;
t=1:128;
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
