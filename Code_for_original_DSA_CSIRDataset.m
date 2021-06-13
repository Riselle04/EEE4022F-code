clear all;
close all;

c = 299792458;

% Load CSIR dataset file Esperance ship 
EsperanceDataset = load('Esp_CAP77-58INB_1_P728_B1_STC76_HRR');

% Obtain the relevant data from the Esperance Dataset
RadarDataset.HRR_profiles = EsperanceDataset.Settings.HRR.HRR_calib_velcomppeak_win.';
[NumberofHRRprofiles, NumberofPulsesinISARBurst] = size(RadarDataset.HRR_profiles);
frequence_step_MHz = EsperanceDataset.Settings.Pattern.FStep(2)*1e6;
RangeResolution_m = c/(2*NumberofPulsesinISARBurst*frequence_step_MHz);
RadarDataset.Range_axis = (0:1:(NumberofPulsesinISARBurst-1))*RangeResolution_m;
RadarDataset.BRF_Hz = 154;

WindowMatrix = EsperanceDataset.Settings.HRR.Win_Func;
sb_Calib_Matrix =  EsperanceDataset.Settings.HRR.corr_filter_vector;
HRR_profiles = ifft(EsperanceDataset.Settings.HRR.TgtBin_velcomp_peak.*WindowMatrix.*sb_Calib_Matrix).';
xaxis = RadarDataset.Range_axis;
yaxis = (1:1:size(RadarDataset.HRR_profiles,1));
    
% Plot of Non-aligned High Range Resolution (HRR) profiles
figure;
axes('fontsize',12);
imagesc(20*log10(abs(HRR_profiles)));
xlabel('Range (bins)');
ylabel('Number of HRR profiles');
title('All HRR profiles: non-aligned');
colorbar;
colormap('jet');
axis xy;
    
%Select ISAR image generation parameters
CPTWL = 128;
StartProfile = 200;  
StopProfile = StartProfile + CPTWL -1;

OrderOfFit = 1;
overlap = 0.5*CPTWL;  
hop = CPTWL - overlap;
noFrames=1+floor(((size(HRR_profiles,1)-StartProfile)-CPTWL)/hop);


% Plot subset of non-aligned High Range Resolution (HRR) profiles plotted
HRRProfiles_Subset = HRR_profiles(StartProfile:StopProfile,:);  
HRRProfiles2Process = circshift(HRRProfiles_Subset, [0 -30]);  
figure;
axes('fontsize',12);
imagesc(20*log10(abs(HRRProfiles2Process)));
xlabel('Range (bins)');
ylabel('Number of HRR profiles');
title('Subset of HRR profiles: non-aligned');
colorbar;
colormap('jet');
axis xy;
    
  
  for i=1:noFrames %loops for the number of images generated
      
      % extract subset of range profiles and perform translational motion compensation
      HRRProfiles_Subset = HRR_profiles(StartProfile:StopProfile,:); 
      alignedHRR_profiles = Haywood_RangeAlignment(HRRProfiles_Subset,OrderOfFit);
      phaseCorrected_HRR_profiles = Haywood_Autofocus(alignedHRR_profiles);
      
      %windowing and ISAR image generation operations
      WindowFunction = repmat(hamming(CPTWL),1,size(phaseCorrected_HRR_profiles,2));    
      windowed_HRRP = phaseCorrected_HRR_profiles.*WindowFunction; 
      windowed_HRRP_noPC = alignedHRR_profiles.*WindowFunction;
      ISAR_image_linear = fftshift(fft(windowed_HRRP,[],1),1); % performs FFT to obtain ISAR image 
      ISAR_image_linear_noAuto= fftshift(fft(windowed_HRRP_noPC,[],1),1);
      
      %axes of ISAR image
      n=size(ISAR_image_linear,1);  
      f = ((-CPTWL/2:1:(CPTWL/2-1))*RadarDataset.BRF_Hz/n);  
     
      % displaying ISAR image with autofocus
      figure
      max_ISAR_dB = max(max(20*log10(abs(ISAR_image_linear))));      
      DynamicRange_dB = 30;                                          
      clims = [ (max_ISAR_dB - DynamicRange_dB) max_ISAR_dB];        
      imagesc(xaxis, f, 20*log10(abs(ISAR_image_linear)),clims);     
      colormap('jet');                                              
      xlabel('Range (m)'); 
      ylabel('Doppler Frequency (Hz)');  
      colorbar;
      drawnow; 
      pause; %click enter to proceed                                                        
      
      % displays ISAR image without autofocus
      figure
      max_ISAR_dB = max(max(20*log10(abs(ISAR_image_linear_noAuto))));      
      DynamicRange_dB = 30;                                        
      clims = [ (max_ISAR_dB - DynamicRange_dB) max_ISAR_dB];        
      imagesc(xaxis, f, 20*log10(abs(ISAR_image_linear_noAuto)),clims);     
      colormap('jet');                                            
      title(sprintf('ISAR image without autofocus %i for profiles %i:%i',i,StartProfile,StopProfile));
      xlabel('Range (m)'); 
      ylabel('Doppler Frequency (Hz)');  
      colorbar;
      drawnow; 
      pause;                                                               
      
      % change interval used for extracting range profiles
      StartProfile = StartProfile + hop;
      StopProfile = StopProfile + hop;
 end


% plotting range aligned subset of profiles
 figure; 
 AlignedHRRProfiles2Process = circshift(alignedHRR_profiles, [0 -15]); 
 axes('fontsize',12);
 imagesc(20*log10(abs(AlignedHRRProfiles2Process)));
 xlabel('Range (bins)');
 ylabel('Number of HRR profiles');
 title('Subset of HRR profiles: Aligned');
 colorbar;
 colormap('jet');
 axis xy;
 
function alignedHRR_profiles=Haywood_RangeAlignment(HRR_profiles, OrderOfFit)

N = size(HRR_profiles,1); %total number of range profiles
n = size(HRR_profiles,2); %total number of range bins
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

%obtains non-integer values for time delay using curve fitting
x = 1:N;
y = timeD.';
coeff = polyfit(x,y,OrderOfFit);
outp = polyval(coeff,x);

 for j=1:N   
     rshift = exp(-1i*(2*pi*outp(j).*m_1/n)); %determining range shift needed for profiles
     alignedHRR_profiles(j,:)=ifft(rshift.*fft(HRR_profiles(j,:)));% applying shift in the frequency domain

 end

end

function phaseCorrected_HRR_profiles = Haywood_Autofocus(alignedHRR_profiles)
n = size(alignedHRR_profiles,2); 

var_arr= var(abs(alignedHRR_profiles));  
lowest=max(var_arr);
 
% average power calculation for all HRRP
power_arr = abs(sum(alignedHRR_profiles.^2,1)); %vector containing the power for each range bin

%calculation for reference range bin
for m=1:n %loops through each range bin
    power_m = power_arr(m); 
    overal_p = (sum(power_arr)-power_m)/(n-1); %calculates overal average power excluding power of selected range bin
     
     if (var_arr(m)<lowest)&&( power_m > overal_p)% checks that DSA criteria is met
          lowest = var_arr(m);
          ref_rb =m; %reference range bin chosen as one that meets the above criteria
     end

end

%determing and applying required phase correction 
 conj_ref_profiles = conj(alignedHRR_profiles(:,ref_rb)); 
 CompensationMatrix = repmat(conj_ref_profiles,1,n);    
 phaseCorrected_HRR_profiles = alignedHRR_profiles.*CompensationMatrix; 

end
