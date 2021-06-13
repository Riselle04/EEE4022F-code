clear all;
close all;

c = 299792458;

% Load CSIR dataset for Esperance vessel
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
    
 
  CPTWL = 128; % in samples
  StartProfile = 200;  
  StopProfile = StartProfile + CPTWL -1;
  
  OrderOfFit = 1;
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
     
      %performs translational motion compensation tehcniques and windowing on data subset
      HRRProfiles_Subset = HRR_profiles(StartProfile:StopProfile,:); 
      alignedHRR_profiles= Haywood_RangeAlignment(HRRProfiles_Subset, OrderOfFit); 
      optimal_ref_rb = Haywood_Autofocus_Modified(alignedHRR_profiles); %finds best reference range bin
      phaseCorrected_HRR_profiles= phaseCorrection(optimal_ref_rb,alignedHRR_profiles); %DSA phase correction performed using optimal reference range bin
    
      WindowFunction = repmat(hamming(CPTWL),1,size(phaseCorrected_HRR_profiles,2));    
      win_x = phaseCorrected_HRR_profiles.*WindowFunction; %windows extracted segment      
      ISAR_image_linear = fftshift(fft(win_x,[],1),1); %performs DFT on windowed segment to generate ISAR image 
      
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
      ylabel('Doppler Frequency (Hz)');  
      colorbar;
      drawnow; 
      pause;  %press enter to continue                                                      
      
      % determines next profile interval to be extracted 
      StartProfile = StartProfile + hop;
      StopProfile = StopProfile + hop;
 end


 
function alignedHRR_profiles = Haywood_RangeAlignment(HRR_profiles, OrderOfFit)
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

%using polyfit and polyval to obtain non-integer values for time delay
x = 1:N;
y = timeD.';
coeff = polyfit(x,y,OrderOfFit);
outp = polyval(coeff,x);

 for j=1:N   % looping through HRR_profiles to perform range alignment
     rshift = exp(-1i*(2*pi*outp(j).*m_1/n)); %calculates required range shift
     alignedHRR_profiles(j,:)=ifft(rshift.*fft(HRR_profiles(j,:)));% applying shift in the frequency domain

 end

end

function optimal_ref_rb = Haywood_Autofocus_Modified(alignedHRR_profiles)
 stop =0;
 threshold_p =0; 
 prev_rb =0;
 count=1;
    while stop ==0  %generates array of possible reference range bins
        
        ref_rb = reference_rangebin(threshold_p,alignedHRR_profiles); %uses original DSA and additional power threshold criteria
            
           if (ref_rb~=prev_rb) && (ref_rb~=0) %prevents repeated range bins from being stored
                rb_arr(count)=ref_rb;
                threshold_vector(count) = threshold_p;  
                count=count+1;
                prev_rb = ref_rb;
            end
            
            if(ref_rb ==0) %stops loop when no reference range bin is found for a given power threshold
               stop =1;
            end
       
        threshold_p =threshold_p+0.01;  %increasing power threshold
    end  
    fprintf('Possible reference range bins: ');disp(rb_arr);
    fprintf('Threshold values: ');disp(threshold_vector);    
    
     for i=1:size(rb_arr,2)  %phase corrects and determines image entropy for all reference range bins found
         phaseCorrected_HRR_profiles = phaseCorrection(rb_arr(i),alignedHRR_profiles);
         WindowFunction = repmat(hamming(size(phaseCorrected_HRR_profiles,1)), 1, size(phaseCorrected_HRR_profiles,2));
         ISARimage = fftshift(fft(phaseCorrected_HRR_profiles.*WindowFunction, [], 1),1);
         IE(i) = ImageEntropy(ISARimage);
      
     end
     %Determines the dominant scatterer/reference range bin as one that yields the lowest image entropy
     fprintf('Image Entropy: ');disp(IE);
     pos_minIE = find(IE==min(IE));
     optimal_ref_rb = rb_arr(pos_minIE); 
     
     %comparing ISAR image entropy before autofocus and with modified DSA
     WindowFunction = repmat(hamming(size(phaseCorrected_HRR_profiles,1)), 1, size(phaseCorrected_HRR_profiles,2));
     ISARimage = fftshift(fft(alignedHRR_profiles.*WindowFunction, [], 1),1);
     IE_NoAutofocus = ImageEntropy(ISARimage);
     fprintf('Image Entropy with no autofocus: ');disp(IE_NoAutofocus);   
     
     
end
    
function ref_rb = reference_rangebin(threshold_p,alignedHRR_profiles)
n = size(alignedHRR_profiles,2); 
var_arr= var(abs(alignedHRR_profiles));  
lowest=max(var_arr);
ref_rb =0;
power_arr = abs(sum(alignedHRR_profiles.^2,1)); %vector containing the power for each range bin
max_power  = max(power_arr);

%calculation for reference range bin

    for m=1:n %loops through each range bin
        power_m = power_arr(m); %extracts the power for range bin m
        overal_p = (sum(power_arr)-power_m)/(n-1); % calculates total average power for the rest of the range bins excluding the current one
       
         if (var_arr(m)<lowest)&&( power_m >overal_p)% checks that DS criteria is met

             if(power_m>threshold_p*max_power) %ensures weaker scatterer is not picked
              lowest = var_arr(m);
              ref_rb =m; %reference range bin chosen as one that meets the above criteria and has the lowest variance
             end   
         end 
    end 
    
end

function phaseCorrected_HRR_profiles= phaseCorrection (ref_rb,alignedHRR_profiles)
 %performs phase correction for optimal reference range bin
 n = size(alignedHRR_profiles,2); 
 conj_ref_profiles = conj(alignedHRR_profiles(:,ref_rb)); 
 CompensationMatrix = repmat(conj_ref_profiles,1,n);    
 phaseCorrected_HRR_profiles = alignedHRR_profiles.*CompensationMatrix; 

end

function IE = ImageEntropy(ISARimage)
%calculates image entropy using given equation
ISARimage_Norm = abs(ISARimage).^2/(sum(sum(abs(ISARimage).^2)));
IE = -sum(sum(ISARimage_Norm.*log(ISARimage_Norm)));

end
