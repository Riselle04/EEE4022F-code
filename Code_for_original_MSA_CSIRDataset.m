
clear all;
close all;

c = 299792458;

% Load CSIR dataset for Esperance ship

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
    
% Plot Non-aligned High Range Resolution (HRR) profiles
figure;
axes('fontsize',12);
imagesc(20*log10(abs(HRR_profiles)));
xlabel('Range (bins)');
ylabel('Number of HRR profiles');
title('All HRR profiles: non-aligned');
colorbar;
colormap('jet');
axis xy;

%Obtain ISAR image generation parameters
CPTWL = 128; % in samples
StartProfile = 200;  
StopProfile = StartProfile + CPTWL -1;

OrderOfFit = 1;
scattererSelValue = 0.16;
overlap = 0.5*CPTWL;  
hop = CPTWL - overlap;
noFrames=1+floor(((size(HRR_profiles,1)-StartProfile)-CPTWL)/hop);


% Plot subset of non-aligned High Range Resolution (HRR) profiles
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

  
  for i=1:noFrames
     
      %applying translational motion compensation and windowing techniques to data subset
      HRRProfiles_Subset = HRR_profiles(StartProfile:StopProfile,:); 
      alignedHRR_profiles = Haywood_RangeAlignment(HRRProfiles_Subset,OrderOfFit);
      phaseCorrected_HRR_profiles = Yuan_Autofocus(alignedHRR_profiles,scattererSelValue);
    
      WindowFunction = repmat(hamming(CPTWL),1,size(phaseCorrected_HRR_profiles,2));    
      win_x = phaseCorrected_HRR_profiles.*WindowFunction; %windows extracted segment      
      ISAR_image_linear = fftshift(fft(win_x,[],1),1); %performs DFT on windowed segment  
      
      %displaying ISAR image
      n=size(ISAR_image_linear,1);  
      f = ((-CPTWL/2:1:(CPTWL/2-1))*frequence_step_MHz/n);  
     
      figure
      max_ISAR_dB = max(max(20*log10(abs(ISAR_image_linear))));     
      DynamicRange_dB = 30;                                         
      clims = [ (max_ISAR_dB - DynamicRange_dB) max_ISAR_dB];        
      imagesc(xaxis, f, 20*log10(abs(ISAR_image_linear)),clims);     
      colormap('jet');                                             
      xlabel('Range (m)'); 
      ylabel('Frequency (Hz)');  
      colorbar;
      drawnow; 
      pause;  % Press enter to proceed                                                     
                                                                    
                                                                    
      %determining next data subset interval
      StartProfile = StartProfile + hop;
      StopProfile = StopProfile + hop;
 end


 
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

 for j=1:N   % looping through HRR_profiles
     rshift = exp(-1i*(2*pi*outp(j).*m_1/n)); % determines range shifts required 
     alignedHRR_profiles(j,:)=ifft(rshift.*fft(HRR_profiles(j,:)));% applying shift in the frequency domain

 end

end

function phaseCorrected_HRR_profiles = Yuan_Autofocus(alignedHRR_profiles,scattererSelValue)

phaseShift=0;
phaseShift(1)=0;

%determine the mean and variance of all range bins
mean_rb = mean(abs(alignedHRR_profiles),1); 
var_rb  = var(abs(alignedHRR_profiles));
ScattererArr = (var_rb)./(var_rb+(mean_rb.^2)); %selection criteria for choosing multiple scatterers  

sel_rangeBins = find(ScattererArr<scattererSelValue); %condition for selecting dominant scatterers
no_selRangeBins = length(sel_rangeBins);

%creating array only containing range profiles of selected scatterers
for i = 1:size(sel_rangeBins,2)
    sel_rangeBinsArr(:,i) = alignedHRR_profiles(:,sel_rangeBins(i));
end

ref_rb = sel_rangeBinsArr(1,:);

%determining phase difference for each range profile w.r.t. refrence range profile
for rangeProfile = 2:size(sel_rangeBinsArr,1) %looping through for all N range profiles
    product = conj(ref_rb).*sel_rangeBinsArr(rangeProfile,:);
    avgScatterers = sum(product)/size(product,2); 
    phaseShift(rangeProfile) = angle(avgScatterers); %extract only the phase of the complex values and use this for compensation
end

%applynig phase correction
CorrectionVector = exp(-1i*phaseShift.');  
CompensationMatrix = repmat(CorrectionVector, 1, size(alignedHRR_profiles,2));
phaseCorrected_HRR_profiles = alignedHRR_profiles.*CompensationMatrix;
end

function IC = ImageContrast(ISARimage_linear)

B = abs((ISARimage_linear)).^2; 
C_df= mean(B,1);
C = mean(C_df);
E = sqrt(mean(mean( (B-C).^2))); 
IC = abs(E/C);

end

function IE = ImageEntropy(ISARimage)

ISARimage_Norm = abs(ISARimage).^2/(sum(sum(abs(ISARimage).^2)));
IE = -sum(sum(ISARimage_Norm.*log(ISARimage_Norm)));

end