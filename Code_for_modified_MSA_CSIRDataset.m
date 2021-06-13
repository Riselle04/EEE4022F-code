
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
overlap = 0.5*CPTWL;  
hop = CPTWL - overlap;
noFrames=1+floor(((size(HRR_profiles,1)-StartProfile)-CPTWL)/hop);

% parameters set by the user for modified MSA 
sel_increment = 0.0001; 
min_rb =2;  
max_rb =18; 
  
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
      optimal_scattererSelValue = choosing_scattererSelValue(alignedHRR_profiles,sel_increment,min_rb,max_rb);
      
      %Yuan phase correction using the optimal scattererSelValue
      phaseCorrected_HRR_profiles = Yuan_Autofocus(alignedHRR_profiles,optimal_scattererSelValue);
      
      WindowFunction = repmat(hamming(CPTWL),1,size(phaseCorrected_HRR_profiles,2));    
      win_x = phaseCorrected_HRR_profiles.*WindowFunction; %windows extracted segment      
      ISAR_image_linear = fftshift(fft(win_x,[],1),1); %performs DFT on windowed segment  
      
      %displaying ISAR image
      n=size(ISAR_image_linear,1);  
      f = ((-CPTWL/2:1:(CPTWL/2-1))*RadarDataset.BRF_Hz/n);  
     
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
      pause;   %click enter to proceed                                                   
                                                                
      %determines the next segment interval
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

 for j=1:N   % looping through HRR_profiles to perform range alignment
     rshift = exp(-1i*(2*pi*outp(j).*m_1/n));  %determining range shift required
     alignedHRR_profiles(j,:)=ifft(rshift.*fft(HRR_profiles(j,:)));% applying shift in the frequency domain

 end

end

function  optimal_scattererSelValue = choosing_scattererSelValue(alignedHRR_profiles,sel_increment,min_rb,max_rb)
%Calculate mean and variance of range bins
mean_rb = mean(abs(alignedHRR_profiles),1); 
var_rb  = var(abs(alignedHRR_profiles));
ScattererArr = (var_rb.^2)./((var_rb.^2)+mean_rb.^2);
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
        IE_selValue  = determingImageQuality(sel_rangeBins,alignedHRR_profiles)
            if (IE_selValue ~= prev_IE_selValue)
                IE_arr(count) = IE_selValue;
                sel_value_arr(count) = scattererSelValue;
                prev_IE_selValue = IE_selValue;
                count = count +1;
            end
        end
    end  
    scattererSelValue = scattererSelValue + sel_increment; %increasing selection value
end
fprintf('IE array');disp(IE_arr);
fprintf('selection value array');disp(sel_value_arr);

%choosing the best scattererSelValue that yields the lowest image entropy
minIE_pos = find(IE_arr==min(IE_arr))
optimal_scattererSelValue = min(sel_value_arr(minIE_pos))

end

function IE_selValue  = determingImageQuality(sel_rangeBins,alignedHRR_profiles)
phaseShift=0;
phaseShift(1)=0;

for i = 1:size(sel_rangeBins,2) %creating array with only the HRRP of the selected range bins 
    sel_rangeBinsArr(:,i) = alignedHRR_profiles(:,sel_rangeBins(i));
end

ref_rb = sel_rangeBinsArr(1,:);

%Yuan autofocus phase correction with range bins selected from the chosen selection value
for rangeProfile = 2:size(sel_rangeBinsArr,1)
    product = conj(ref_rb).*sel_rangeBinsArr(rangeProfile,:);
    avgScatterers = sum(product)/size(product,2);
    phaseShift(rangeProfile) = avgScatterers;
end

%phase correction for specific scatterer selection value
phaseShift = phaseShift.';
CompensationMatrix = repmat(phaseShift, 1, size(alignedHRR_profiles,2));
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

ScattererArr = (var_rb.^2)./((var_rb.^2)+mean_rb.^2);

sel_rangeBins = find(ScattererArr<scattererSelValue); %condition for selecting dominant scatterers
no_selRangeBins = length(sel_rangeBins);

%creating array only containing range profiles of selected scatterers
for i = 1:size(sel_rangeBins,2)
    sel_rangeBinsArr(:,i) = alignedHRR_profiles(:,sel_rangeBins(i));
end

ref_rb = sel_rangeBinsArr(1,:);

%determining phase difference for each range profile w.r.t. refrence range profile
for rangeProfile = 2:size(sel_rangeBinsArr,1) 
    product = conj(ref_rb).*sel_rangeBinsArr(rangeProfile,:);
    avgScatterers = sum(product)/size(product,2); 
    phaseShift(rangeProfile) = avgScatterers;
end

%applying phase correction
phaseShift = phaseShift.';
CompensationMatrix = repmat(phaseShift, 1, size(alignedHRR_profiles,2));
phaseCorrected_HRR_profiles = alignedHRR_profiles.*CompensationMatrix;
end

function IE = ImageEntropy(ISARimage)

ISARimage_Norm = abs(ISARimage).^2/(sum(sum(abs(ISARimage).^2)));
IE = -sum(sum(ISARimage_Norm.*log(ISARimage_Norm)));

end