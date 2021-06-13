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

    
  CPTWL = 128; 
  OrderOfFit = 1;
  
  %Array containing start and stop profiles for focused CSIR images
  FocusedImages = [584 2760 3272 3336 3400 4872 4936;711 2887 3399 3463 3527 4999 5063];
  CSIR_imageNo = [7 41 49 50 51 74 75];
  ImageToPlot = 7;   %7 focused images - change variable from 1 -7 to view  images
  
  %Profiles required for selected image
  StartProfile = FocusedImages(1,ImageToPlot);
  StopProfile = FocusedImages(2,ImageToPlot); 
  
  % parameters set by the user for modified MSA 
  sel_increment = 0.001; 
  min_rb =2;  
  max_rb =18; 

  %applying translational motion compensation techniques and windowing process
  HRRProfiles_Subset = HRR_profiles(StartProfile:StopProfile,:);   
  alignedHRR_profiles = Haywood_RangeAlignment_Modified(HRRProfiles_Subset,OrderOfFit);
  optimal_scattererSelValue = choosing_scattererSelValue(alignedHRR_profiles,sel_increment,min_rb,max_rb);
  phaseCorrected_HRR_profiles = Yuan_Autofocus(alignedHRR_profiles,optimal_scattererSelValue);
  
  WindowFunction = repmat(hamming(CPTWL),1,size(phaseCorrected_HRR_profiles,2));    
  win_x = phaseCorrected_HRR_profiles.*WindowFunction; %windows extracted segment      
  ISAR_image_linear = fftshift(fft(win_x,[],1),1); %performs DFT on windowed segment 
  
  IE = ImageEntropy(ISAR_image_linear);
  IC = ImageContrast(ISAR_image_linear);
  
  %shifts required to keep target at the middle of the image
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
    findM = r==maxV;
    timeD(i,:) = lags(findM);
end

%using polyfit and polyval to obtain non-integer values for time delay
x = 1:N;
y = timeD.';
coeff = polyfit(x,y,OrderOfFit);
outp = polyval(coeff,x);

 for j=1:N   % looping through HRR_profiles
     rshift = exp(-1i*(2*pi*outp(j).*m_1/n));  %determining range shifts required
     alignedHRR_profiles(j,:)=ifft(rshift.*fft(HRR_profiles(j,:)));% applying shift in the frequency domain

 end

end

function  optimal_scattererSelValue = choosing_scattererSelValue(alignedHRR_profiles,sel_increment,min_rb,max_rb)
%Calculate mean and variance of range bins
mean_rb = mean(abs(alignedHRR_profiles),1); 
var_rb  = var(abs(alignedHRR_profiles));
ScattererArr = var_rb./(var_rb+(mean_rb.^2));
over_max_rb =0;
scattererSelValue =0;
prev_IE_selValue =0;
count=1;

%Selection and phase correction prcoess for given scattererSelValue 
while over_max_rb~=1
    sel_rangeBins = find(ScattererArr<scattererSelValue);
    no_sel_rb = length(sel_rangeBins); %calculates number of range bins selected
    
    if(no_sel_rb>=min_rb) %only starts selection process if min rangebins are found
       
        if(no_sel_rb>max_rb) %stops generating selection values if too many rangebins have been selected
            over_max_rb=1;
        else
         %Determine image entopy and storing this and scatterer selection value
           IE_selValue = determingImageQuality(sel_rangeBins,alignedHRR_profiles);
            if (IE_selValue ~= prev_IE_selValue)
                IE_arr(count) = IE_selValue;
                sel_value_arr(count) = scattererSelValue;
                prev_IE_selValue = IE_selValue;
                count = count +1;
            end
        end
    end  
    scattererSelValue = scattererSelValue + sel_increment;  %increasing selection value
end


%choosing the best scattererSelValue that yields the lowest image entropy
minIE_pos = find(IE_arr==min(IE_arr));
optimal_scattererSelValue = sel_value_arr(minIE_pos)
end

function [IE_selValue]  = determingImageQuality(sel_rangeBins,alignedHRR_profiles)
phaseShift=0;
phaseShift(1)=0;

for i = 1:size(sel_rangeBins,2)  %creating array with only the HRRP of the selected range bins 
    sel_rangeBinsArr(:,i) = alignedHRR_profiles(:,sel_rangeBins(i));
end

ref_rb = sel_rangeBinsArr(1,:);

%Yuan autofocus phase correction with range bins selected from the inputted selection value
for rangeProfile = 2:size(sel_rangeBinsArr,1) 
    product = conj(ref_rb).*sel_rangeBinsArr(rangeProfile,:);
    avgScatterers = sum(product)/size(product,2);
    phaseShift(rangeProfile) = angle(avgScatterers);
end

%phase correction for specific scatterer selection value
CorrectionVector = exp(-1i*phaseShift.');
CompensationMatrix = repmat(CorrectionVector, 1, size(alignedHRR_profiles,2));
phaseCorrected_HRR_profiles = alignedHRR_profiles.*CompensationMatrix;

%determining IE for specific scatterer selection value
WindowFunction = repmat(hamming(size(phaseCorrected_HRR_profiles,1)), 1, size(phaseCorrected_HRR_profiles,2));
ISARimage = fftshift(fft(phaseCorrected_HRR_profiles.*WindowFunction, [], 1),1);
IE_selValue = ImageEntropy(ISARimage);

end

%Yuan phase correction using the optimal selection value
function phaseCorrected_HRR_profiles = Yuan_Autofocus(alignedHRR_profiles,scattererSelValue)

phaseShift=0;
phaseShift(1)=0;
mean_rb = mean(abs(alignedHRR_profiles),1); 
var_rb  = var(abs(alignedHRR_profiles));

ScattererArr = (var_rb)./(var_rb+(mean_rb.^2));

sel_rangeBins = find(ScattererArr<scattererSelValue); %condition for selecting dominant scatterers
no_selRangeBins = length(sel_rangeBins) 

%creating array only containing range profiles of selected scatterers
for i = 1:size(sel_rangeBins,2)
    sel_rangeBinsArr(:,i) = alignedHRR_profiles(:,sel_rangeBins(i));
end

ref_rb = sel_rangeBinsArr(1,:);

%determining phase difference for each range profile w.r.t. refrence range profile
for rangeProfile = 2:size(sel_rangeBinsArr,1) 
    product = conj(ref_rb).*sel_rangeBinsArr(rangeProfile,:);
    avgScatterers = sum(product)/size(product,2); 
    phaseShift(rangeProfile) = angle(avgScatterers);
end

%applying phase correction
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

