%Note: All steps indicated in code refer to steps for the range alignment,
%Dominant Scatterer Selection Autofocus and ISAR image generation as seen
%in the System Design section of the report.

close all;
clear all; 

StartProfile=1;
StopProfile = 636; 

% Load and process Altera dataset of a ship at sea (i.e.Spanish ship)
Dataset = load('Altera_CWLFM_Dataset.mat');
HRR_profiles = Dataset.HRR_profiles(StartProfile: StopProfile, :);    
ProfileRepetitionFreq = Dataset.ProfileRepitionFreq;
Range_axis = Dataset.Range_axis; 

% Plot all unaligned range profiles
NumProfiles_axis = 1:size(HRR_profiles,1);    
figure; imagesc(Range_axis, NumProfiles_axis, 20*log10(abs(HRR_profiles))); 
colorbar
xlabel('Range (m)');
ylabel('Profile number');
title('Unaligned range profiles');
axis xy;
colormap('jet');

%ISAR image generation process

% Steps 1-4: Calculate image generation parameters
OrderOfFit = 3;
window_size = 256;
overlap = 0.5*window_size;  

hop = window_size - overlap;
noFrames=1+floor((size(HRR_profiles,1)-window_size)/hop); 
nextSeg=0;
 
 for i=1:noFrames  %loops for the number of images generated
    
      %Step 5: Extract subset of range profiles for image generation
      HRR_profilesToAlign = HRR_profiles((1+nextSeg):(window_size+nextSeg), :); 
     
      %Step 6: Perform translational motion compensation techniques
      alignedHRR_profiles= Haywood_RangeAlignment(HRR_profilesToAlign, OrderOfFit);  
      phaseCorrected_HRR_profiles = Haywood_Autofocus(alignedHRR_profiles);
      
      %Step 7: windowing 
      WindowFunction = repmat(hamming(window_size),1,size(phaseCorrected_HRR_profiles,2)); % generates hamming window 
      windowed_HRRP = phaseCorrected_HRR_profiles.*WindowFunction; %windows extracted segment  
      
      % Step 8: Performing FFT for ISAR image generation and displaying ISAR image
      ISAR_image_linear = fftshift(fft(windowed_HRRP,[],1),1); %performs DFT on windowed segment  
      
      %axes for ISAR image
      n=size(ISAR_image_linear,1);  
      f = ((-window_size/2:1:(window_size/2-1))*ProfileRepetitionFreq/n);  
    
      figure
      max_ISAR_dB = max(max(20*log10(abs(ISAR_image_linear))));      
      DynamicRange_dB = 50;                                          
      clims = [ (max_ISAR_dB - DynamicRange_dB) max_ISAR_dB];        
      imagesc(Range_axis, f, 20*log10(abs(ISAR_image_linear)),clims);     
      colormap('jet');                                              
      xlabel('Range(m)')
      ylabel('Doppler Frequency (Hz)');  
      colorbar;
      drawnow;  
      
      % calculate image entropy and contrast for comparing performance
      IE = ImageEntropy(ISAR_image_linear);
      IC = ImageContrast(ISAR_image_linear);
      fprintf('IE for ISAR image %i: \n',IE);
      fprintf('IC for ISAR image %i: \n', IC);
      
      nextSeg = nextSeg+hop;  %Step 9: Determines the next segment interval
 end
 
% Displaying Range-Doppler image without Autofocus for all range profiles
n = size(phaseCorrected_HRR_profiles,1);  
f = ((-n/2:1:(n/2-1))*ProfileRepetitionFreq/n);

WindowFunction = repmat(hamming(size(alignedHRR_profiles,1)), 1, size(alignedHRR_profiles,2)); % multiply by a window function before applying the FFT
ISARimageNoAuto = fftshift(fft(alignedHRR_profiles.*WindowFunction, [], 1),1); 
ISARimageNoAuto_dB = 20*log10(abs(ISARimageNoAuto));

% Plot ISAR image 
figure
max_ISAR_dB = max(max(20*log10(abs(ISARimageNoAuto))));      
DynamicRange_dB = 50;                                          
clims = [ (max_ISAR_dB - DynamicRange_dB) max_ISAR_dB];        
imagesc(Range_axis, f, 20*log10(abs(ISARimageNoAuto)),clims);     
colormap('jet');                                              
xlabel('Range(m)'); 
ylabel('Doppler Frequency (Hz)');  
colorbar;

% Plot aligned HRR profiles
alignedHRR_profiles= Haywood_RangeAlignment(HRR_profilesToAlign, OrderOfFit);
NumProfiles_axis1 = 1:size(alignedHRR_profiles,1); 
figure; imagesc(Range_axis, NumProfiles_axis1, 20*log10(abs(alignedHRR_profiles))); 
colorbar
xlabel('Range (m)');
ylabel('Profile number');
title('Aligned range profiles');
axis xy;
colormap('jet');
% -------------------------------------------------------------------------------------------------------------------------
   
%Haywood Range Alignment technique
function alignedHRR_profiles=Haywood_RangeAlignment(HRR_profiles, OrderOfFit)
N = size(HRR_profiles,1); 
n = size(HRR_profiles,2); 
ref_profile = HRR_profiles(1,:);  

timeD = zeros(N,1);
alignedHRR_profiles =zeros(N,n);
m_1 = 0:(n-1);

 % Steps 1-3:Performing cross-correlation and estimating time delays
for i=1:N
    [r,lags] = xcorr(abs(ref_profile),abs(HRR_profiles(i,:)));
    maxV= max(r);  
    findM = find(r==maxV); %time delay indicated by the peak position after cross-correlation
    timeD(i,:) = lags(findM);
end

% Step 4: Curve fitting process to ontain non-integer values for time delay
x = 1:N;
y = timeD.';
coeff = polyfit(x,y,OrderOfFit);
outp = polyval(coeff,x);

% Step 5-6: Applying relevant range shifts to profiles and storing resulting array
 for j=1:N 
     rshift = exp(-1i*(2*pi*outp(j).*m_1/n));  
     alignedHRR_profiles(j,:)=ifft(rshift.*fft(HRR_profiles(j,:)));
 end

end

% Haywood autofocus technique
function phaseCorrected_HRR_profiles = Haywood_Autofocus(alignedHRR_profiles)

n = size(alignedHRR_profiles,2); %calculates total number of range bins
var_arr= var(abs(alignedHRR_profiles));  %calculates variance array from aligned profiles
lowest=max(var_arr);
 
power_arr = abs(sum(alignedHRR_profiles.^2,1)); %creates array for the power of each range bin
for m=1:n %loops through each range bin
    power_m = power_arr(m);  
    overal_p = (sum(power_arr)-power_m)/(n-1); %calculates overal average power excluding power of selected range bin
    
    % Step 2 - checks minimum variance and adequate power conditions
     if (var_arr(m)<lowest)&&( power_m > overal_p)
          lowest = var_arr(m);
          ref_rb =m; %reference range bin chosen as one that meets the above criteria
           
     end

end
fprintf('ref range bin:');disp(ref_rb);

% Step 3-4: Determines phase compensation vector and performs phase
% correction on range profiles
 conj_ref_profiles = conj(alignedHRR_profiles(:,ref_rb)); 
 CompensationMatrix = repmat(conj_ref_profiles,1,n);   
 phaseCorrected_HRR_profiles = alignedHRR_profiles.*CompensationMatrix;

end

%Image entropy calculations
function IC = ImageContrast(ISARimage_linear)

B = abs((ISARimage_linear)).^2; 
C_df= mean(B,1);
C = mean(C_df);
E = sqrt(mean(mean( (B-C).^2))); 
IC = abs(E/C);

end

%Image contrast calculations
function IE = ImageEntropy(ISARimage)

ISARimage_Norm = abs(ISARimage).^2/(sum(sum(abs(ISARimage).^2)));
IE = -sum(sum(ISARimage_Norm.*log(ISARimage_Norm)));

end



