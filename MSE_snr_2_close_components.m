clc
clear
close all
SampFreq = 128;
addpath('D:\tfsa_5-5\windows\win64_bin');

t = 0:1/SampFreq:1-1/SampFreq;



% COde for 3 component


    Sig1 =1*exp(1i*(2*pi*(5*t +15*t.^3)));
    Sig2 =1*exp(1i*(2*pi*(10*t +15*t.^3)));
    Sig3 =1*exp(1i*(2*pi*(50*t -1*13*t.^3)));
    Sig4=exp(1i*(2*pi*(60*t )));
    IF_O(:,1)=45*t.^2+10;
    IF_O(:,2)=45*t.^2+5;
    %IF_O(:,3)=60;
    SigA=1*Sig1+1*Sig2+0*Sig3+0*Sig4;
    num=2;
IF_O=IF_O/(SampFreq/2);


% HADTFD BASED

iiii=0;
dis=1;
%Sig =1*Sig1+Sig3;


NS=200;
t=16:128-15;

iiii=0;



for snr=-10:2:10
    iiii=iiii+1;
    
    for k1=1:NS
        Sig =SigA;
        
        Sig=awgn(Sig,snr,'measured');
        [tfd,orient]=HTFD_new1(Sig,2,30,84);
        
        %   tfsapl(Sig,tfd);
        for kk=0:1
            if kk==0
                [fmult,out, peaks] = component_linking_new(tfd,orient,0.1,length(tfd)/8,10);
                [fmult]= merge_IFs(fmult,orient,20,30,length(Sig)/2);

            else
                        [tfd,orient]=HTFD_new1(Sig,2,20,48);

                for i=1:4
                [tfd,orient]= post_processing_directional(tfd,2,20,48/1);
                end
                 [fmult,out, peaks] = component_linking_new(tfd,orient,0.1,length(tfd)/8,10);
                [fmult]= merge_IFs(fmult,orient,20,30,length(Sig)/2);
                %[fmult,~, ~] = component_linking(tfd,0.1,64);
            
            end
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
                    if dis==1
                        figure;
                        plot(t,IF(t),'-',t,IF_O(t,b1),'d');
                    end
                end
            end
            msee1(kk+1,k1)=mean(msee);
        end
    end
    mse_adtfd_new(iiii)=mean(msee1(1,:));
    mse_adtfd_old(iiii)=mean(msee1(2,:));
end



snr=-10:2:10;
plot(snr,10*(log10(mse_adtfd_new)),'-.b+','linewidth',4);
hold on;
plot(snr,10*(log10(mse_adtfd_old)),'--md','linewidth',4);
xlabel('Signal to noise ratio');
ylabel('Mean square error in dB');
legend('ADTFD','Iterative ADTFD');