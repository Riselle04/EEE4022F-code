% All steps mentioned in the comments refer to the steps for the MSA
% mentioned in the Design section of the report 
% Note the Yuan_Autofocus refers to the selected MSA autofocus technique
close all;
clear all; 

% Load and process Spanish dataset of a ship at sea
Dataset = load('Altera_CWLFM_Dataset.mat');
HRR_profiles = Dataset.HRR_profiles(1:636,:);
ProfileRepetitionFreq = Dataset.ProfileRepitionFreq;
Range_axis = Dataset.Range_axis; 

%Plot HRR profiles
NumProfiles_axis = 1:size(HRR_profiles,1); 
figure; imagesc(Range_axis, NumProfiles_axis, 20*log10(abs(HRR_profiles))); 
colorbar
xlabel('Range (m)');
ylabel('Profile number');
title('Unaligned range profiles');
axis xy;
colormap('jet');

%Obtain ISAR image generation parameters
StartProfile = 1;  
StopProfile = 256;
CPTWL = StopProfile-StartProfile+1;
OrderOfFit =3;
scattererSelValue = 0.16;

overlap = 0.5*CPTWL;  
hop = CPTWL - overlap;
noFrames=1+floor((size(HRR_profiles,1)-CPTWL)/hop); 
nextSeg=0;
 
 for i=1:noFrames
      
      %applying translational motion compensation techniques to data subset
      HRR_profilesToAlign = HRR_profiles((1+nextSeg):(CPTWL+nextSeg), :); 
      alignedHRR_profiles= Haywood_RangeAlignment(HRR_profilesToAlign, OrderOfFit);  
      phaseCorrected_HRR_profiles = Yuan_Autofocus(alignedHRR_profiles,scattererSelValue);
      
      %windowing and ISAR image generation process
      WindowFunction = repmat(hamming(size(phaseCorrected_HRR_profiles,1)), 1, size(phaseCorrected_HRR_profiles,2)); 
      ISAR_image_linear = fftshift(fft(phaseCorrected_HRR_profiles.*WindowFunction, [], 1),1); 
      
      %displaying ISAR image
      n=size(ISAR_image_linear,1);   
      f = ((-CPTWL/2:1:(CPTWL/2-1))*ProfileRepetitionFreq/n);  
      figure
      max_ISAR_dB = max(max(20*log10(abs(ISAR_image_linear))));      
      DynamicRange_dB = 50;                                          
      clims = [(max_ISAR_dB - DynamicRange_dB) max_ISAR_dB];        
      imagesc(Range_axis, f, 20*log10(abs(ISAR_image_linear)),clims);     
      colormap('jet');      
      xlabel('Range (m)'); 
      ylabel('Frequency (Hz)');  
      colorbar;
      drawnow; 
      
      IE = ImageEntropy(ISAR_image_linear);
      IC = ImageContrast(ISAR_image_linear);
      fprintf('IE for ISAR image %i: \n',IE);
      fprintf('IC for ISAR image %i: \n', IC);
      
      nextSeg = nextSeg+hop;  %determines the next segment interval
 end


%Plot aligned HRR profiles

  HRRProfiles2Align = circshift(alignedHRR_profiles, [0 -30]);   
  figure;
  axes('fontsize',12);
  imagesc(20*log10(abs(HRRProfiles2Align)));
  xlabel('Range (bins)');
  ylabel('Number of HRR profiles');
  title('Aligned Range Profiles');
  colorbar;
  colormap('jet');
  axis xy;
  

function alignedHRR_profiles=Haywood_RangeAlignment(HRR_profiles, OrderOfFit)

N = size(HRR_profiles,1); 
n = size(HRR_profiles,2); 
ref_profile = HRR_profiles(1,:);  
timeD = zeros(N,1);
alignedHRR_profiles =zeros(N,n);
m_1 = 0:(n-1);

for i=1:N
    % returns cross correlation and lags at which the correlations are computed
    [r,lags] = xcorr(abs(ref_profile),abs(HRR_profiles(i,:)));
    maxV= max(r);
    findM = find(r==maxV);
    timeD(i,:) = lags(findM);
end

%using polyfit and polyval to obtain non-integer values for time delay
x = 1:N;
y = timeD.';
coeff = polyfit(x,y,OrderOfFit);
outp = polyval(coeff,x);

 for j=1:N   % performs range alignment
     rshift = exp(-1i*(2*pi*outp(j).*m_1/n));  % determines range shifts required
     alignedHRR_profiles(j,:)=ifft(rshift.*fft(HRR_profiles(j,:)));% applying shift in the frequency domain

 end

end

function phaseCorrected_HRR_profiles = Yuan_Autofocus(alignedHRR_profiles,scattererSelValue)

phaseShift=0;
phaseShift(1)=0;

%Step 1: determine the mean and variance of all range bins
mean_rb = mean(abs(alignedHRR_profiles),1); 
var_rb  = var(abs(alignedHRR_profiles));
ScattererArr = (var_rb)./((var_rb)+mean_rb.^2); %selection criteria for choosing multiple scatterers        

%calculating number of scatterers that meet selection criterion
sel_rangeBins = find(ScattererArr<scattererSelValue); %condition for selecting dominant scatterers
no_selRangeBins = length(sel_rangeBins); 
disp('Number of selected range bins:'); disp(no_selRangeBins);

%Step 2: creating array only containing range profiles of selected scatterers
for i = 1:size(sel_rangeBins,2)
    sel_rangeBinsArr(:,i) = alignedHRR_profiles(:,sel_rangeBins(i));
end

ref_rb = sel_rangeBinsArr(1,:); 

%Steps 3-4: determining phase difference for each range profile w.r.t. refrence range profile
for rangeProfile = 2:size(sel_rangeBinsArr,1) %looping through for all N range profiles
    product = conj(ref_rb).*sel_rangeBinsArr(rangeProfile,:);
    avgScatterers = sum(product)/size(product,2); 
    phaseShift(rangeProfile) = angle(avgScatterers); %extract only the phase of the complex values and use this for compensation
                                        
end

%Step 5: applying phase correction on range profiles
CorrectionVector = exp(-1i*phaseShift.');                                        
CompensationMatrix = repmat(CorrectionVector, 1, size(alignedHRR_profiles,2));   
phaseCorrected_HRR_profiles = alignedHRR_profiles.*CompensationMatrix;
end

function IE = ImageEntropy(ISARimage)
%calculates image entropy using given equation
ISARimage_Norm = abs(ISARimage).^2/(sum(sum(abs(ISARimage).^2)));
IE = -sum(sum(ISARimage_Norm.*log(ISARimage_Norm)));
end

function IC = ImageContrast(ISARimage_linear)
%calculates image contrast using given equation
B = abs((ISARimage_linear)).^2; 
C_df= mean(B,1);
C = mean(C_df);
E = sqrt(mean(mean( (B-C).^2))); 
IC = abs(E/C);

end