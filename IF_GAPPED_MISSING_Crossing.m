clc
clear
close all
SampFreq = 128;
addpath('D:\tfsa_5-5\windows\win64_bin');

t = 0:1/SampFreq:1-1/SampFreq;



% COde for 3 component


Sig1 =1*exp(1i*(2*pi*(5*t +15*t.^3)));
Sig2 =1*exp(1i*(2*pi*(12*t +15*t.^3)));
Sig3 =1*exp(1i*(2*pi*(50*t -1*13*t.^3)));
Sig4=exp(1i*(2*pi*(60*t )));
IF_O(:,1)=45*t.^2+12;
IF_O(:,1)=45*t.^2+5;
IF_O(:,2)=-39*t.^2+50;
%IF_O(:,3)=60;
SigA=1*Sig1+0*Sig2+1*Sig3+0*Sig4;
num=2;
IF_O=IF_O/(SampFreq/2);


iiii=0;
t=16:128-15;
jjjj=0;
NS=200/1;
%IF_O(:,3)=90*t.^2/2+15;
for snr=0:10:30
    jjjj=jjjj+1;
    iiii=0;
    for N_S=4:4:16
        iiii=iiii+1;
        
        for k1=1:NS
            Sig=SigA;
            Sig=awgn(Sig,snr,'measured');
            p=[];
            for i=1:4
                pp = 32*(i-1)+ randperm(32-N_S-1,1);
                p1=pp:1:pp+N_S;
                p=[ p p1];
            end
            Sig(p)=0;
            [NA]=find(Sig~=0);
            for kkkkk=0:1
                
                % ORIGINAL
                delta=5;
                alpha = 5;
                if kkkkk==0   %ADTFD+VITERBI
                    [tfd,orient]=HTFD_new1(Sig,2,20,48);
                    
                    for i=1:4
                        [tfd,orient]= post_processing_directional(tfd,2,20,48/1);
                    end
                    
                elseif kkkkk==1 %the new algorithm
                    [tfd,orient]=HTFD_new1(Sig,2,30,84);
                end
                [fmult,out, peaks] = component_linking_new(tfd,orient,0.1,length(tfd)/4,10);
                [fmult]= merge_IFs(fmult,orient,20,20,length(Sig)/2);
                findex=fill_zeros(fmult);
                findex1=zeros(num,length(Sig));
                [aa,~]=size(findex);
                
                for ii=1:aa
                    findex1(ii,1:length(findex))=findex(ii,:);
                    
                end
                findex=findex1;
                msee=0.1*ones(1,num);
                dis=0;
                [aa,~]=size(findex);
                %findex=findex/(SampFreq/2);
                
                for ii=1:num
                    if ii<=aa
                        
                        IF=findex(ii,:)/(length(Sig));
                        %t=t(5:end-5);
                        for i=1:num
                            c(i)=sum(abs(IF(t).'-IF_O(t,i)).^2);
                        end
                        [a1 b1]=min(c);
                        if msee(b1)>=a1(1)/length(t)
                            msee(b1)=a1(1)/length(t);
                        end
                        
                    end
                end
               % msee1(kk+1,k1)=mean(msee);
            
            if kkkkk==0
                mse_IADTFD1(k1)=mean(msee);
                
            elseif kkkkk==1
                mse_ADTFD1(k1)=mean(msee);
                
            end
            
            
            
        end
        
        
    end
    mse_IADTFD(jjjj,iiii)=mean(mse_IADTFD1);
    mse_ADTFD(jjjj,iiii)=mean(mse_ADTFD1);
    
    end
mse_ADTFD(jjjj,:)
mse_IADTFD(jjjj,:)
end

figure;
plot(4:4:16,mse_IADTFD(1,:),'k','linewidth',3);
hold on;
plot(4:4:16,mse_ADTFD(1,:),'r','linewidth',3);
xlabel('Number of missing samples in each gap')
ylabel('Mean absolute errror')
title('a')
legend('IADTFD','ADTFD');
set(gca,'fontsize', 24)

figure;
plot(4:4:16,mse_IADTFD(2,:),'k','linewidth',3);
hold on;
plot(4:4:16,mse_ADTFD(2,:),'r','linewidth',3);
xlabel('Number of missing samples in each gap')
ylabel('Mean absolute errror')
title('b')
legend('IADTFD','ADTFD');
set(gca,'fontsize', 24)

figure;
plot(4:4:16,mse_IADTFD(3,:),'k','linewidth',3);
hold on;
plot(4:4:16,mse_ADTFD(3,:),'r','linewidth',3);
xlabel('Number of missing samples in each gap')
ylabel('Mean absolute errror')
title('c')
legend('IADTFD','ADTFD');
set(gca,'fontsize', 24)

figure;
plot(4:4:16,mse_IADTFD(4,:),'k','linewidth',3);
hold on;
plot(4:4:16,mse_ADTFD(4,:),'r','linewidth',3);
xlabel('Number of missing samples in each gap')
ylabel('Mean absolute errror')
title('d')
legend('IADTFD','ADTFD');
set(gca,'fontsize', 24)

