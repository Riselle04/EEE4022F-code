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
  scattererSelValue = 0.16;
  %Array containing start and stop profiles for focused CSIR images
  FocusedImages = [584 2760 3272 3336 3400 4872 4936;711 2887 3399 3463 3527 4999 5063];
  CSIR_imageNo = [7 41 49 50 51 74 75];
  ImageToPlot = 1;   %7 focused images - change variable from 1 -7 to view  images
                    
  StartProfile = FocusedImages(1,ImageToPlot);
  StopProfile = FocusedImages(2,ImageToPlot);     
  
  %applying translational motion compensation and windowing techniques to data subset
  HRRProfiles_Subset = HRR_profiles(StartProfile:StopProfile,:);   
  alignedHRR_profiles = Haywood_RangeAlignment_Modified(HRRProfiles_Subset,OrderOfFit);
  phaseCorrected_HRR_profiles = Yuan_Autofocus(alignedHRR_profiles,scattererSelValue);
  
  WindowFunction = repmat(hamming(CPTWL),1,size(phaseCorrected_HRR_profiles,2));    
  win_x = phaseCorrected_HRR_profiles.*WindowFunction; %windows extracted segment      
  ISAR_image_linear = fftshift(fft(win_x,[],1),1); %performs DFT on windowed segment 
  
  IE = ImageEntropy(ISAR_image_linear);
  IC = ImageContrast(ISAR_image_linear);
  
  %determining requierd shifts to center the target on image
     if (ImageToPlot == 1)
        ISAR_image_linear = circshift(ISAR_image_linear, [0 -20]);       
   elseif (ImageToPlot == 2)
        ISAR_image_linear = circshift(ISAR_image_linear, [0 25]);       
   elseif (ImageToPlot == 3)
        ISAR_image_linear = circshift(ISAR_image_linear, [0 25]);       
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
     
  %displaying ISAR image
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
  

function alignedHRR_profiles=Haywood_RangeAlignment_Modified(HRR_profiles, OrderOfFit)
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
     rshift = exp(-1i*(2*pi*outp(j).*m_1/n));   % determines range shifts required 
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

sel_rangeBins = find(ScattererArr<scattererSelValue); 
no_selRangeBins = length(sel_rangeBins); %number of range bins selected

%creating array only containing range profiles of selected scatterers
for i = 1:size(sel_rangeBins,2)
    sel_rangeBinsArr(:,i) = alignedHRR_profiles(:,sel_rangeBins(i));
end

ref_rb = sel_rangeBinsArr(1,:);

%determining phase difference for each range profile w.r.t. refrence range profile
for rangeProfile = 2:size(sel_rangeBinsArr,1) 
    product = conj(ref_rb).*sel_rangeBinsArr(rangeProfile,:);
    avgScatterers = sum(product)/size(product,2); 
    phaseShift(rangeProfile) = angle(avgScatterers);%extract only the phase of the complex values and use this for compensation
end

%applynig phase correction
CorrectionVector = exp(-1i*phaseShift.');
CompensationMatrix = repmat(CorrectionVector, 1, size(alignedHRR_profiles,2));
phaseCorrected_HRR_profiles = alignedHRR_profiles.*CompensationMatrix;
end

function IE = ImageEntropy(ISARimage)

ISARimage_Norm = abs(ISARimage).^2/(sum(sum(abs(ISARimage).^2)));
IE = -sum(sum(ISARimage_Norm.*log(ISARimage_Norm)));
end

function IC = ImageContrast(ISARimage_linear)

B = abs((ISARimage_linear)).^2; 
C_df= mean(B,1);
C = mean(C_df);
E = sqrt(mean(mean( (B-C).^2))); 
IC = abs(E/C);

end