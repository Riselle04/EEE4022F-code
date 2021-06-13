 %Notes: The steps commented in the code refers to the steps used for the Dominant Scatterer
 %Selection algorithm in the System Design section of report
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
  ImageToPlot =1;   %7 focused images - change variable from 1 -7 to view  images
  
  % extracts profiles needed for selected image
  StartProfile = FocusedImages(1,ImageToPlot);
  StopProfile = FocusedImages(2,ImageToPlot);     
 
  %performs translational motion compensation tehcniques and windowing on
  %data subset
  HRRProfiles_Subset = HRR_profiles(StartProfile:StopProfile,:); 
  alignedHRR_profiles= Haywood_RangeAlignment(HRRProfiles_Subset, OrderOfFit); 
  optimal_ref_rb = Haywood_Autofocus_Modified(alignedHRR_profiles); %finds best/optimal reference range bin
  phaseCorrected_HRR_profiles= phaseCorrection(optimal_ref_rb,alignedHRR_profiles);

  WindowFunction = repmat(hamming(CPTWL),1,size(phaseCorrected_HRR_profiles,2));    
  win_x = phaseCorrected_HRR_profiles.*WindowFunction; %windows extracted segment      
  ISAR_image_linear = fftshift(fft(win_x,[],1),1); %performs DFT on windowed segment  

 
   %shifts to move target to center of the image
     if (ImageToPlot == 1)
        ISAR_image_linear = circshift(ISAR_image_linear, [0 -20]);       
   elseif (ImageToPlot == 2)
        ISAR_image_linear = circshift(ISAR_image_linear, [0 25]);      
   elseif (ImageToPlot == 3)
        ISAR_image_linear = circshift(ISAR_image_linear, [0 20]);       
   elseif (ImageToPlot == 4)
       ISAR_image_linear = circshift(ISAR_image_linear, [0 20]);       
   elseif (ImageToPlot == 5)
        ISAR_image_linear = circshift(ISAR_image_linear, [0 15]);       
   elseif (ImageToPlot == 6)
        ISAR_image_linear = circshift(ISAR_image_linear, [0 -5]);       
   elseif (ImageToPlot == 7)
       ISAR_image_linear = circshift(ISAR_image_linear, [0 -5]);       
   elseif (ImageToPlot == 8)   
       ISAR_image_linear = circshift(ISAR_image_linear, [0 -5]);       
     end
   
  %plotting range-dopler map
  n=size(ISAR_image_linear,1);  
  f = ((-CPTWL/2:1:(CPTWL/2-1))*RadarDataset.BRF_Hz/n); 
 
  figure
  max_ISAR_dB = max(max(20*log10(abs(ISAR_image_linear)))); 
  ISAR_dB = 20*log10(abs(ISAR_image_linear));                   
  ISAR_dB  = ISAR_dB  - max_ISAR_dB;      %normalise maximum of ISAR image to 0 dB
  DynamicRange_dB = 30;                                          
  clims = [ -DynamicRange_dB  0];                              
  imagesc(xaxis, f, ISAR_dB ,clims);     
  colormap('jet');                                              
  xlabel('Range (m)');
  ylabel('Doppler frequency (Hz)');                           
  colorbar;

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

 for j=1:N   % performs range alignment
     rshift = exp(-1i*(2*pi*outp(j).*m_1/n)); %calculates required range shift
     alignedHRR_profiles(j,:)=ifft(rshift.*fft(HRR_profiles(j,:)));% applying shift in the frequency domain

 end

end

function optimal_ref_rb = Haywood_Autofocus_Modified(alignedHRR_profiles)
 stop =0;
 threshold_p =0; 
 prev_rb =0;
 count=1;
    while stop ==0  %Step 1-3: generates array of eligible reference range bins
        
        ref_rb = reference_rangebin(threshold_p,alignedHRR_profiles); %uses original DSA and additional power threshold criteria
            
           if (ref_rb~=prev_rb) && (ref_rb~=0)  %prevents repeated range bins from being stored
                rb_arr(count)=ref_rb;
                threshold_vector(count) = threshold_p; 
                count=count+1;
                prev_rb = ref_rb;
            end
            
            if(ref_rb ==0) %stops loop when no reference range bin is found for a given power threshold
               stop =1;
            end
       
        threshold_p =threshold_p+0.1; %increasing power threshold
    end  
    fprintf('Possible reference range bins: ');disp(rb_arr);
    fprintf('Threshold values: ');disp(threshold_vector);    
     
    %Steps 4-5: phase corrects and determines image entropy for all reference range bins found
     for i=1:size(rb_arr,2)  
         phaseCorrected_HRR_profiles = phaseCorrection(rb_arr(i),alignedHRR_profiles);
         WindowFunction = repmat(hamming(size(phaseCorrected_HRR_profiles,1)), 1, size(phaseCorrected_HRR_profiles,2));
         ISARimage = fftshift(fft(phaseCorrected_HRR_profiles.*WindowFunction, [], 1),1);
         IE(i) = ImageEntropy(ISARimage);
      
     end
     
     %Step 6: Determines the dominant scatterer/reference range bin as one
     %that yields the lowest image entropy
     pos_minIE = find(IE==min(IE));
     optimal_ref_rb = rb_arr(pos_minIE); 
     
     %comparing ISAR image entropy before autofocus and with modified DSA
     WindowFunction = repmat(hamming(size(phaseCorrected_HRR_profiles,1)), 1, size(phaseCorrected_HRR_profiles,2));
     ISARimage = fftshift(fft(alignedHRR_profiles.*WindowFunction, [], 1),1);
     IE_NoAutofocus = ImageEntropy(ISARimage);
     fprintf('Image Entropy with no autofocus: ');disp(IE_NoAutofocus); 
     fprintf('Image Entropy after DS autofocus: ');disp(IE(pos_minIE));
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
       
         if (var_arr(m)<lowest)&&( power_m > overal_p)% checks that DS criteria is met
             if(power_m>threshold_p*max_power) % additional criteria to ensure a weaker scatterer is not picked
              lowest = var_arr(m);
              ref_rb =m; %reference range bin chosen as one that meets the above criteria and has the lowest variance
             end   
         end 
    end 
    
end

function phaseCorrected_HRR_profiles= phaseCorrection (ref_rb,alignedHRR_profiles)
 %Step 7: performs phase correction for optimal reference range bin
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

