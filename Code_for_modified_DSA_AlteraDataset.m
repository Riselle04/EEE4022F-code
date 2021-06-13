% Steps commented refer to steps for Dominant Scatterer Selection Algorithm in System Design section of report 
close all;
clear all; 

StartProfile=1;
StopProfile = 636; 

% Load and process Spanish dataset of a ship at sea
Dataset = load('Altera_CWLFM_Dataset.mat');
HRR_profiles = Dataset.HRR_profiles(StartProfile: StopProfile, :);    

ProfileRepetitionFreq = Dataset.ProfileRepitionFreq;
Range_axis = Dataset.Range_axis; 

% Plot all HRR profiles
NumProfiles_axis = 1:size(HRR_profiles,1);  
figure; imagesc(Range_axis, NumProfiles_axis, 20*log10(abs(HRR_profiles))); 
colorbar
xlabel('Range (m)');
ylabel('Profile number');
title('Unaligned range profiles');
axis xy;
colormap('jet');
 
%Obtain ISAR image generation parameters
OrderOfFit = 3;
window_size = 256;
overlap = 0.5*window_size;  

 hop = window_size - overlap;
 noFrames=1+floor((size(HRR_profiles,1)-window_size)/hop); 
 nextSeg=0;
 
 for i=1:noFrames
      %obtaining data subset
      HRR_profilesToAlign = HRR_profiles((1+nextSeg):(window_size+nextSeg), :); 
      
      %applying translational motion compensation techniques to data subset
      alignedHRR_profiles= Haywood_RangeAlignment(HRR_profilesToAlign, OrderOfFit); 
      optimal_ref_rb = Haywood_Autofocus_Modified(alignedHRR_profiles); %finds best/optimal reference range bin
      phaseCorrected_HRR_profiles= phaseCorrection(optimal_ref_rb,alignedHRR_profiles); %DSA phase correction done using optimal reference range bin
      
      %windowing and ISAR image generation process
      WindowFunction = repmat(hamming(window_size),1,size(phaseCorrected_HRR_profiles,2)); 
      windowed_HRRP = phaseCorrected_HRR_profiles.*WindowFunction; %windows extracted segment  
      ISAR_image_linear = fftshift(fft(windowed_HRRP,[],1),1); %performs DFT on windowed segment  
    
      %displaying ISAR image
      n=size(ISAR_image_linear,1);  
      f = ((-window_size/2:1:(window_size/2-1))*ProfileRepetitionFreq/n);  
     
      figure
      imagesc(Range_axis, f, 20*log10(abs(ISAR_image_linear))); 
      colormap('jet');
      xlabel('Range (m)'); 
      ylabel('Doppler Frequency (Hz)');  
      colorbar;
      drawnow; 
       
      nextSeg = nextSeg+hop; %determining next data subset interval
 end
 

% -------------------------------------------------------------------------------------------------------
      
function alignedHRR_profiles = Haywood_RangeAlignment(HRR_profiles, OrderOfFit)

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
     rshift = exp(-1i*(2*pi*outp(j).*m_1/n)); % determines range shifts required
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
            
           if (ref_rb~=prev_rb) && (ref_rb~=0) % prevents repeated range bins from being stored
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
        power_m = power_arr(m);
        overal_p = (sum(power_arr)-power_m)/(n-1); % calculates total average power for the rest of the range bins excluding the current one
       
         if (var_arr(m)<lowest)&&( power_m > overal_p)% checks that DS criteria is met
             if(power_m>threshold_p*max_power) %ensures weaker scatterer is not picked
              lowest = var_arr(m);
              ref_rb =m; %reference range bin chosen as one that meets the above criteria and has the lowest variance
             end   
         end %if statement checking DS criteria
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

function IC = ImageContrast(ISARimage_linear)
%calculates image contrast using given equation
B = abs((ISARimage_linear)).^2; 
C_df= mean(B,1);
C = mean(C_df);
E = sqrt(mean(mean( (B-C).^2))); 
IC = abs(E/C);
end
