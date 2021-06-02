%--------------------------------------------------------------------------
% Newborn EEG classification
%
% input:    Newborn EEG data
% output:   Display the classification results
%
% By Dr. Nabeel Ali Khan
%--------------------------------------------------------------------------

%%-------------------------------------------------------------------------
%% Newborn EEG Database
%%
%% The N.mat (S.mat) file contains 50 non-seizure (seizure) segments.
%% The segments have been inspected visually and picked.
%% The EEG segments have been band-pass ?ltered in the range 0.5 to 10 Hz
%% and down-sampled from 256 Hz to 20 Hz.
%%-------------------------------------------------------------------------

clear all;
N_C=4;
load seizure_samples
load non_seizure_samples

nS=200;
%addpath 'E:\TFSA7\TFSA7\'

cfeatures_vector_S=[];
cfeatures_vector_N=[];

art=0;
for j=1:1:nS
    for class=1:2
        
        %%--------------------------------------------
        %% Time-frequency signal representation
        %%--------------------------------------------
        if class==1 signal=seizure_samples(j,:); end
        if class==2 signal=non_seizure_samples(j,:); end
       % signal=decimate(signal,2);
        signalf=filter([1 -1],1,signal);
        
        % Code to estimate Spectral Centroids%
        signal_freq=abs(fft(hilbert(signal)));
        signal_freq=signal_freq(1:128);
        mean_freq=sum(signal_freq.*(0:127))/sum(signal_freq); %Spectral Centroid
        mean_freq1=sum(signal_freq(1:round(mean_freq)).*(0:round(mean_freq)-1))/sum(signal_freq(1:round(mean_freq)));
        mean_freq2=sum(signal_freq(1:round(mean_freq1)).*(0:round(mean_freq1)-1))/sum(signal_freq(1:round(mean_freq1)));
        
        % Code to estimate Spectral Flatness%
        signal_freq(signal_freq==0)=eps;
%        SF=prod(abs(signal_freq).^(1/256))/sum(abs(signal_freq));


%        [I,O]=HTFD_new_EEG(signalf,2,20,64);%48
%Inew= post_processing_directional(wvv,2,20,64);
Inew=I;
for i=1:4
[Inew,orient]= post_processing_directional_EEG(Inew,2,20,64);%48
end
O=O*3;
orient=orient*3;
I1=zeros(size(Inew));
I1(or(orient<15,orient>165))=Inew(or(orient<15,orient>165));
I2=zeros(size(Inew));
I2(and(orient>80,orient<100))=Inew(and(orient>80,orient<100));
(sum(I1(:))+sum(I2(:)))/sum(Inew(:))
TFC(1)=(sum(I1(:))+sum(I2(:)))/sum(Inew(:));
Inew(Inew==0)=eps;

SF=prod(abs(Inew).^(1/256))/sum(abs(Inew));


        %TFC(1)=min(mean(std(IF.')),mean(std(GD.'))); % Proposed feature
        TFC(2)=entropy(signal);
        TFC(3)=mean_freq; % Spectral Centroid
        TFC(4)=mean(abs(diff(sign(signalf))));  %Zero crossing rate
        TFC(5)=iqr(signalf) ; % Inter Quartile Range
        TFC(6)=mean_freq1; % Spectral Centroid 1
        TFC(7)=mean_freq2; %Spectral centroid 2
        TFC(8)=SF; % Spectral FLux
        TFC(9)=skewness(signal);
        TFC(10)=sum(I1(:))/sum(Inew(:));
        TFC(11)=sum(I2(:))/sum(Inew(:));
 I1=zeros(size(I));
I1(or(orient<15,orient>165))=I(or(orient<15,orient>165));
I2=zeros(size(Inew));
I2(and(orient>80,orient<100))=I(and(orient>80,orient<100));
(sum(I1(:))+sum(I2(:)))/sum(I(:))
TFC(12)=       (sum(I1(:))+sum(I2(:)))/sum(Inew(:));
        
        %% feature extraction
        
        if class==1
           
            features_vector_S(j,:)=  TFC;
        else
            features_vector_N(j,:)=  TFC;
            
        end
    end
end




F=[features_vector_S;features_vector_N];

mask=[zeros(1,nS) ones(1,nS)];
for k = 1:size(F,2)
    N_thresh = 100000; % Number of the thresholding levels
    [Sen,Spe] = roc_rates_function(F(:,k),mask,N_thresh);
    auc(k) = trapz(1-Spe,Sen);
    
    if auc(k)<0.5
        auc(k) = 1 - auc(k);
    else
    end
end
auc(1)
Cross_validation_leave_one(features_vector_N(:,[1]),features_vector_S(:,[1]),1)

%save('Seizure_samples','Seizure_samples');
%save('Back_samples','Back_samples');
