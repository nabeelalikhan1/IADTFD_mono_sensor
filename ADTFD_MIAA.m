
close all;
clear all;
addpath 'D:\tfsa_5-5\windows\win64_bin\'

M=128;
n=0:M-1;
%x=cos(0.0005*2*pi*n.^3/M+2*2*pi*n/M)+cos(0.0005*2*pi*n.^3/M+30*2*pi*n/M);%+exp(1i*0.001*2*pi*n.^3/M+1i*15*2*pi*n/M);
%x=cos(0.001*2*pi*n.^3/M+2*2*pi*n/M)+cos(-0.001*2*pi*n.^3/M+45*2*pi*n/M);%+exp(1i*0.001*2*pi*n.^3/M+1i*15*2*pi*n/M);
%x=exp(1i*0.0005*2*pi*n.^3/M+1i*2*2*pi*n/M)+exp(1i*0.0005*2*pi*n.^3/M+1i*10*2*pi*n/M)+exp(1i*52*2*pi*n/M);
%x=exp(1i*0.001*2*pi*n.^3/M+1i*2*2*pi*n/M)+exp(-1i*0.001*2*pi*n.^3/M+1i*45*2*pi*n/M);

%x=exp(1i*0.001*2*pi*n.^3/M+1i*2*2*pi*n/M)+exp(-1i*0.001*2*pi*n.^3/M+1i*45*2*pi*n/M);% +exp(1i*55*2*pi*n/M);
x=1*exp(2*pi*1i*(0.3*n.^3/(3*M*M)+0.006*n.^2/(2*M)+0.08*n))+1*exp(2*pi*1i*(0.3*n.^3/(3*M*M)+0.002*n.^2/(2*M)+0.05*n));
%x=exp(1i*0.001*2*pi*n.^3/M+1i*4*pi*n/M)+0*exp(1i*0.0007*2*pi*n.^3/M+1i*8*2*pi*n/M)+1*exp(-1i*0.001*2*pi*n.^3/M+1i*0.95*pi*n);

N=length(x);
iii= randsample(128,32) ;

%x(iii)=0;

                  I1=HTFD_new1(x,2,20,64/1);
I=I1;
for i=1:10
   I= post_processing_directional(I,2,12,32);
 if i==1
   II=I;
end
end
  %for iii=1:5
 %   Inew= post_processing_directional(adtfd,2,20,64);
%end
%I=adtfd;
figure;tfsapl(x,I1);
figure;tfsapl(x,II);

figure;tfsapl(x,I);
