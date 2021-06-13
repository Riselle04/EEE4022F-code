  clear all;
  close all;

  c = 299792458;

  % Load file
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

    
  CPTWL = 128; 
  OrderOfFit = 1;
  
  %Array containing start and stop profiles for focused CSIR images
  FocusedImages = [584 2760 3272 3336 3400 4872 4936;711 2887 3399 3463 3527 4999 5063];
  CSIR_imageNo = [7 41 49 50 51 74 75];
  ImageToPlot = 1;   %7 selected images - change variable from 1 -7 to view  images
   
  %extracting required profiles for selected image
  StartProfile = FocusedImages(1,ImageToPlot); 
  StopProfile = FocusedImages(2,ImageToPlot);     
  
  %Applying translational motion compensation techniques, windowing and DFT
  HRRProfiles_Subset = HRR_profiles(StartProfile:StopProfile,:);  
  alignedHRR_profiles = Haywood_RangeAlignment(HRRProfiles_Subset,OrderOfFit);
  phaseCorrected_HRR_profiles = Haywood_Autofocus(alignedHRR_profiles);

  WindowFunction = repmat(hamming(CPTWL),1,size(phaseCorrected_HRR_profiles,2));    
  windowed_HRRP = phaseCorrected_HRR_profiles.*WindowFunction; %windows extracted segment      
  ISAR_image_linear = fftshift(fft(windowed_HRRP,[],1),1); %performs DFT on windowed segment 
  
  %range shifts required to move target to the center of the image
     if (ImageToPlot == 1)
        ISAR_image_linear = circshift(ISAR_image_linear, [0 -20]);       
   elseif (ImageToPlot == 2)
        ISAR_image_linear = circshift(ISAR_image_linear, [0 25]);      
   elseif (ImageToPlot == 3)
        ISAR_image_linear = circshift(ISAR_image_linear, [0 20]);       
   elseif (ImageToPlot == 4)
       ISAR_image_linear = circshift(ISAR_image_linear, [0 20]);       
   elseif (ImageToPlot == 5)
        ISAR_image_linear = circshift(ISAR_image_linear, [0 20]);       
   elseif (ImageToPlot == 6)
        ISAR_image_linear = circshift(ISAR_image_linear, [0 -5]);       
   elseif (ImageToPlot == 7)
       ISAR_image_linear = circshift(ISAR_image_linear, [0 -5]);       
   elseif (ImageToPlot == 8)   
       ISAR_image_linear = circshift(ISAR_image_linear, [0 -5]);       
     end
   
  %Displaying ISAR image
  n=size(ISAR_image_linear,1);  
  f = ((-CPTWL/2:1:(CPTWL/2-1))*RadarDataset.BRF_Hz/n);          

  figure
  max_ISAR_dB = max(max(20*log10(abs(ISAR_image_linear)))); 
  ISAR_dB = 20*log10(abs(ISAR_image_linear));                   
  ISAR_dB  = ISAR_dB  - max_ISAR_dB;      
  DynamicRange_dB = 30;                                          
  clims = [ -DynamicRange_dB  0];                              
  imagesc(xaxis, f, ISAR_dB ,clims);     
  colormap('jet');                                              
  xlabel('Range (m)');                                         
  ylabel('Doppler frequency (Hz)');                           
  colorbar;
  

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

 for j=1:N   % looping through range profiles
     rshift = exp(-1i*(2*pi*outp(j).*m_1/n)); %calculating required range shift for range profiles 
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
     if (var_arr(m)<lowest)&&( power_m > overal_p)% checks that DS criteria is met
          lowest = var_arr(m);
          ref_rb =m; %reference range bin chosen as one that meets the above criteria
       
     end

end

 % calculating and applying phase correction
 conj_ref_profiles = conj(alignedHRR_profiles(:,ref_rb)); 
 CompensationMatrix = repmat(conj_ref_profiles,1,n);    
 phaseCorrected_HRR_profiles = alignedHRR_profiles.*CompensationMatrix; 

end
