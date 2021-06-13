%The steps referred to in the comments are for the modified MSA technique
%mentioned in the System Design section of the report
close all;
clear all; 

% Load and process Spanish dataset of a ship at sea
Dataset = load('Altera_CWLFM_Dataset.mat');
HRR_profiles = Dataset.HRR_profiles(1:636,:);
ProfileRepetitionFreq = Dataset.ProfileRepitionFreq;
Range_axis = Dataset.Range_axis; 

%Obtain ISAR image generation parameters
StartProfile = 1;  
StopProfile = 256;
CPTWL = StopProfile-StartProfile+1;
OrderOfFit =3;
overlap = 0.5*CPTWL;  
hop = CPTWL - overlap;
noFrames=1+floor((size(HRR_profiles,1)-CPTWL)/hop); 
nextSeg=0;

%Step 1: Selection of required parameters for modified MSA
sel_increment = 0.0001; 
min_rb =2;  
max_rb =18; 
 
 for i=1:noFrames
    
      %applying translational motion compensation and windowing techniques to data subset
      HRR_profilesToAlign = HRR_profiles((1+nextSeg):(CPTWL+nextSeg), :); 
      alignedHRR_profiles= Haywood_RangeAlignment_Modified(HRR_profilesToAlign, OrderOfFit);  
      optimal_scattererSelValue = choosing_scattererSelValue(alignedHRR_profiles,sel_increment,min_rb,max_rb);
      
      %Yuan phase correction using the optimal scattererSelValue
      phaseCorrected_HRR_profiles = Yuan_Autofocus(alignedHRR_profiles,optimal_scattererSelValue); 
    
      WindowFunction = repmat(hamming(size(phaseCorrected_HRR_profiles,1)), 1, size(phaseCorrected_HRR_profiles,2)); %windows extracted segment  
      ISAR_image_linear = fftshift(fft(phaseCorrected_HRR_profiles.*WindowFunction, [], 1),1); %performs DFT on windowed segment 
     
      %displaying ISAR image
      n=size(ISAR_image_linear,1);   
      f = ((-CPTWL/2:1:(CPTWL/2-1))*ProfileRepetitionFreq/n);  

      figure
      max_ISAR_dB = max(max(20*log10(abs(ISAR_image_linear))));      
      DynamicRange_dB = 50;                                          
      clims = [(max_ISAR_dB - DynamicRange_dB) max_ISAR_dB];        
      imagesc(Range_axis, f, 20*log10(abs(ISAR_image_linear)),clims);     
      colormap('jet');      
      xlabel('Range (m)'); 
      ylabel('Doppler Frequency (Hz)');  
      colorbar;
      drawnow; 
      
      
      IE = ImageEntropy(ISAR_image_linear);
      IC = ImageContrast(ISAR_image_linear);
      fprintf('IE for ISAR image %i: \n',IE);
      fprintf('IC for ISAR image %i: \n', IC);
      
      nextSeg = nextSeg+hop;  %determines the next segment interval
 end
     

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

 for j=1:N   % looping through HRR_profiles to perform range alignment
     rshift = exp(-1i*(2*pi*outp(j).*m_1/n));  %determining range shift required
     alignedHRR_profiles(j,:)=ifft(rshift.*fft(HRR_profiles(j,:)));% applying shift in the frequency domain

 end

end

function  optimal_scattererSelValue = choosing_scattererSelValue(alignedHRR_profiles,sel_increment,min_rb,max_rb)
%Step 2: Calculate mean and variance of range bins
mean_rb = mean(abs(alignedHRR_profiles),1); 
var_rb  = var(abs(alignedHRR_profiles));
ScattererArr = var_rb./(var_rb+(mean_rb.^2));
over_max_rb =0;
scattererSelValue =0;
prev_IE_selValue =0;
count=1;

%Steps 3-5: Selection and phase correction prcoess for given scattererSelValue 
while over_max_rb~=1 
    sel_rangeBins = find(ScattererArr<scattererSelValue);
    no_sel_rb = length(sel_rangeBins);
    if(no_sel_rb>=min_rb) %only starts selection process if min range bins are found
       
        if(no_sel_rb>max_rb) %stops generating selection values if max number range bins have been selected
            over_max_rb=1;
        
        % Step 6: Determine image entopy and storing this and scatterer
        % selection value
        else 
             IE_selValue  = determingImageQuality(sel_rangeBins,alignedHRR_profiles);
            if (IE_selValue ~= prev_IE_selValue) %storing unique range bin values (i.e. no repetition)
                IE_arr(count) = IE_selValue;
                sel_value_arr(count) = scattererSelValue;
                no_sel_rb_arr(count) = no_sel_rb;
                prev_IE_selValue = IE_selValue;
                count = count +1;
            end
        end
    end  
    scattererSelValue = scattererSelValue + sel_increment; %Step 7: increasing selection value
end

%Step 9: choosing the best scattererSelValue that yields the lowest image
%entropy
minIE = min(IE_arr);
optimal_scattererSelValue = min(sel_value_arr(find(IE_arr==min(minIE))))
end

function [IE_selValue]  = determingImageQuality(sel_rangeBins,alignedHRR_profiles) 
phaseShift=0;
phaseShift(1)=0;

for i = 1:size(sel_rangeBins,2) %creating array with only the HRRP of the selected range bins 
    sel_rangeBinsArr(:,i) = alignedHRR_profiles(:,sel_rangeBins(i));
end

ref_rb = sel_rangeBinsArr(1,:);

%Yuan autofocus phase correction with range bins selected from the specific selection value
for rangeProfile = 2:size(sel_rangeBinsArr,1) 
    product = conj(ref_rb).*sel_rangeBinsArr(rangeProfile,:);
    avgScatterers = sum(product)/size(product,2);
    phaseShift(rangeProfile) = angle(avgScatterers);
end

CorrectionVector = exp(-1i*phaseShift.');
CompensationMatrix = repmat(CorrectionVector, 1, size(alignedHRR_profiles,2));
phaseCorrected_HRR_profiles = alignedHRR_profiles.*CompensationMatrix;

%determining IE for specific scatterer selection value
WindowFunction = repmat(hamming(size(phaseCorrected_HRR_profiles,1)), 1, size(phaseCorrected_HRR_profiles,2));
ISARimage = fftshift(fft(phaseCorrected_HRR_profiles.*WindowFunction, [], 1),1);
IE_selValue = ImageEntropy(ISARimage);


end

%Step 9: Yuan phase correction using optimal selection value
function phaseCorrected_HRR_profiles = Yuan_Autofocus(alignedHRR_profiles,scattererSelValue)

phaseShift=0;
phaseShift(1)=0;
mean_rb = mean(abs(alignedHRR_profiles),1); 
var_rb  = var(abs(alignedHRR_profiles));

ScattererArr = var_rb./(var_rb+(mean_rb.^2));

sel_rangeBins = find(ScattererArr<scattererSelValue); %condition for selecting dominant scatterers
no_selRangeBins = length(sel_rangeBins); %calculates number of selected range bins for selection value

%creating array only containing range profiles of selected scatterers
for i = 1:size(sel_rangeBins,2)
    sel_rangeBinsArr(:,i) = alignedHRR_profiles(:,sel_rangeBins(i));
end

ref_rb = sel_rangeBinsArr(1,:);

%determining phase difference for each range profile w.r.t. refrence range profile
for rangeProfile = 2:size(sel_rangeBinsArr,1) %looping through for all N range profiles
    product = conj(ref_rb).*sel_rangeBinsArr(rangeProfile,:);
    avgScatterers = sum(product)/size(product,2); 
    phaseShift(rangeProfile) = angle(avgScatterers);
end

%applying phase correction
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