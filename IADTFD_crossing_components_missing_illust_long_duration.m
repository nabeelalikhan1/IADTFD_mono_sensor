
close all;
clear all;
addpath 'D:\tfsa_5-5\windows\win64_bin\'
t=0:1/512:1-1/512;
f=0:0.5:256-0.5;
M=128;
n=0:M-1;

Sig1 =1*exp(1i*(2*pi*(5*t +65*t.^3)));
Sig2 =1*exp(1i*(2*pi*(12*t +15*t.^3)));
Sig3 =1*exp(1i*(2*pi*(250*t -1*65*t.^3)));
Sig4=exp(1i*(2*pi*(60*t )));
IF_O(:,1)=45*t.^2+12;
IF_O(:,1)=45*t.^2+5;
IF_O(:,2)=-39*t.^2+50;
%IF_O(:,3)=60;
SigA=1*Sig1+0*Sig2+1*Sig3+0*Sig4;


x=SigA;
N=length(x);
iii= randsample(512,400) ;
x(iii)=0;
Ispec=quadtfd(x,length(x)-1,1,'specx',145,'hamm',512);
Ispec=Ispec/sum(abs(Ispec(:)));
Imb=quadtfd(x,length(x)/4-1,1,'mb',0.051,512);

Imb(Imb<0)=0;
Imb=Imb/sum(abs(Imb(:)));Imb(Imb<0)=0;
Imb=Imb/sum(abs(Imb(:)));
%x(iii)=0;
Fre_Grid=1*length(x);
Fre_Line  = 2*pi*(1:Fre_Grid)/Fre_Grid;
%I1=min(HTFD_new1(x,2,30,64/1),HTFD_new1(x,2,50,84/1));
I1=HTFD_new1(x,2,25,64*4);
I1=I1/sum(abs(I1(:)));
I=I1;
for i=1:5
    I= post_processing_directional(I,2,25,64*4);
    if i==1
        II=I;
        II(II<0)=0;
        II=II/sum(abs(II(:)));
    end
end
I(I<0)=0;
II(II<0)=0;
I=I/sum(abs(I(:)));
t=1:1:512;
f=0:1/1024:0.5-1/1024;
figure; imagesc(t,f,I1)

set(gcf,'Position',[20 100 640 500]);
xlabel('Time / Sec','FontSize',20,'FontName','Times New Roman');
ylabel('Frequency / Hz','FontSize',20,'FontName','Times New Roman');
set(gca,'YDir','normal');
title('(a)','FontSize',24,'FontName','Times New Roman');
set(gca,'FontSize',20);

figure; imagesc(t,f,I)

set(gcf,'Position',[20 100 640 500]);
xlabel('Time / Sec','FontSize',20,'FontName','Times New Roman');
ylabel('Frequency / Hz','FontSize',20,'FontName','Times New Roman');
set(gca,'YDir','normal');
title('(b)','FontSize',24,'FontName','Times New Roman');
set(gca,'FontSize',20);

sum(sum(abs(I.^0.5))).^2
sum(sum(abs(I1.^0.5))).^2

