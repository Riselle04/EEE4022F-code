%This file generates all the comparison plots discussed in the Comparison
%of Algorithms section in the report (Chapter 8)
clear all
close all

%Spanish dataset comparison 
Spanish_ISARimages = [1 2 3]; 
noAuto_Spanish_IE = [5.8294 5.5362 5.9336];
noAuto_Spanish_IC = [50.4345 73.2605 51.5562];

DSA_Spanish_IE = [6.0425 6.8056 6.5152];
DSA_Spanish_IC = [49.3827 45.3526 48.3100];

Yuan_Spanish_IE = [6.6747 6.0144 6.0372];
Yuan_Spanish_IC = [38.7201 71.3585 63.9173];

modDSA_Spanish_IE = [5.8254 5.5558 5.8353];
modDSA_Spanish_IC = [51.3751 70.8210 58.6858];

modYuan_Spanish_IE = [5.8177 5.4972 5.8231];
modYuan_Spanish_IC = [51.6238 81.1219 66.9051];


%IE Spanish data plots
%IE comparison for all algorithms
figure
plot(Spanish_ISARimages, noAuto_Spanish_IE, 'r');
hold on 
plot(Spanish_ISARimages, DSA_Spanish_IE, 'b');
hold on
plot(Spanish_ISARimages,Yuan_Spanish_IE,'g');
hold on 
plot(Spanish_ISARimages, modDSA_Spanish_IE,'m');
hold on
plot(Spanish_ISARimages,modYuan_Spanish_IE,'c');
hold off
xlabel('ISAR image');
ylabel('Image Entropy');
legend('No Autofocus','DSA','MSA','Modified DSA', 'Modified MSA');
ax = gca;
ax.XTick = unique( round(ax.XTick) );

%IE for Yuan MSA and Haywood DSA
figure
plot(Spanish_ISARimages, noAuto_Spanish_IE, 'r');
hold on
plot(Spanish_ISARimages, DSA_Spanish_IE, 'b');
hold on
plot(Spanish_ISARimages,Yuan_Spanish_IE,'g');
hold off
xlabel('ISAR image');
ylabel('Image Entropy');
legend('No Autofocus','DSA','MSA');
ax = gca;
ax.XTick = unique( round(ax.XTick) );


%IE for modMSA and modDSA
figure
plot(Spanish_ISARimages, noAuto_Spanish_IE, 'r');
hold on
plot(Spanish_ISARimages, modDSA_Spanish_IE, 'b');
hold on
plot(Spanish_ISARimages,modYuan_Spanish_IE,'g');
hold off
xlabel('ISAR image');
ylabel('Image Entropy');
legend('No Autofocus','Modified DSA','Modified MSA');
ax = gca;
ax.XTick = unique( round(ax.XTick) );



%-----------------------------------------------------------------------------------------------------------------

%IC Spanish data plots

%IC comparison for all algorithms
figure
plot(Spanish_ISARimages, noAuto_Spanish_IC, 'r');
hold on 
plot(Spanish_ISARimages, DSA_Spanish_IC, 'b');
hold on
plot(Spanish_ISARimages,Yuan_Spanish_IC,'g');
hold on 
plot(Spanish_ISARimages, modDSA_Spanish_IC,'m');
hold on
plot(Spanish_ISARimages,modYuan_Spanish_IC,'c');
hold off
xlabel('ISAR image');
ylabel('Image Contrast');
legend('No Autofocus','DSA','MSA','Modified DSA', 'Modified MSA');
ax = gca;
ax.XTick = unique( round(ax.XTick) );

%IC for Yuan MSA and Haywood DSA
figure
plot(Spanish_ISARimages, noAuto_Spanish_IC, 'r');
hold on
plot(Spanish_ISARimages, DSA_Spanish_IC, 'b');
hold on
plot(Spanish_ISARimages,Yuan_Spanish_IC,'g');
hold off
xlabel('ISAR image');
ylabel('Image Contrast');
legend('No Autofocus','DSA autofocus','MSA autofocus');
ax = gca;
ax.XTick = unique( round(ax.XTick) );


%IC for modMSA and modDSA
figure
plot(Spanish_ISARimages, noAuto_Spanish_IC, 'r');
hold on
plot(Spanish_ISARimages, modDSA_Spanish_IC, 'b');
hold on
plot(Spanish_ISARimages,modYuan_Spanish_IC,'g');
hold off
xlabel('ISAR image');
ylabel('Image Contrast');
legend('No Autofocus','Modified DSA','Modified MSA');
ax = gca;
ax.XTick = unique( round(ax.XTick) );




%--------------------------------------------------------------------------------------------------------------------------


%CSIR dataset

CSIR_ISARimages = [1 2 3 4 5 6 7]; 
noAuto_CSIR_IE = [5.9242 5.1557 5.5470 6.0481 6.0656 6.0448 6.0132];
noAuto_CSIR_IC = [6.7027 10.6515 10.6169 6.6226 6.5958 6.7240 6.7498];

DSA_CSIR_IE = [5.6715 5.2063 5.4811 5.6596 5.3380 5.3505 5.1955];
DSA_CSIR_IC = [7.9858 10.4235 11.6434 8.9858 10.0338 12.8454 12.6805];

Yuan_CSIR_IE = [5.3244 5.2648 5.9592 5.5618 5.7258 5.3336 5.6113];
Yuan_CSIR_IC = [10.0633 9.7532 9.9283 11.9138 10.3798 13.1761 11.5400];

modDSA_CSIR_IE = [5.4958 4.6587 5.4811 5.3598 5.2639 5.3505 5.1955];
modDSA_CSIR_IC = [8.6913 16.0796 11.6434 13.2367 14.2504 12.8454 12.6805];

modYuan_CSIR_IE = [5.3126 5.1322 5.4749 5.5169 5.5197 5.3110 5.1742];
modYuan_CSIR_IC = [10.2218 10.8939 11.5389 11.7271 11.1589 13.3557 13.0271];


%IE CSIR data plots
%IE comparison for all algorithms
figure
plot(CSIR_ISARimages, noAuto_CSIR_IE, 'r');
hold on 
plot(CSIR_ISARimages, DSA_CSIR_IE, 'b');
hold on
plot(CSIR_ISARimages,Yuan_CSIR_IE,'g');
hold on 
plot(CSIR_ISARimages, modDSA_CSIR_IE,'m');
hold on
plot(CSIR_ISARimages,modYuan_CSIR_IE,'c');
hold off
xlabel('ISAR image');
ylabel('Image Entropy');
legend('No Autofocus','DSA','MSA','Modified DSA', 'Modified MSA');

%IE for Yuan MSA and Haywood DSA
figure
plot(CSIR_ISARimages, noAuto_CSIR_IE, 'r');
hold on
plot(CSIR_ISARimages, DSA_CSIR_IE, 'b');
hold on
plot(CSIR_ISARimages,Yuan_CSIR_IE,'g');
hold off
xlabel('ISAR image');
ylabel('Image Entropy');
legend('No Autofocus','DSA','MSA');


%IE for modMSA and modDSA
figure
plot(CSIR_ISARimages, noAuto_CSIR_IE, 'r');
hold on
plot(CSIR_ISARimages, modDSA_CSIR_IE, 'b');
hold on
plot(CSIR_ISARimages,modYuan_CSIR_IE,'g');
hold off
xlabel('ISAR image');
ylabel('Image Entropy');
legend('No Autofocus','Modified DSA','Modified MSA');
%-------------------------------------------------------------------------------------------------------------------------

%IC plots

%IC comparison for all algorithms
figure
plot(CSIR_ISARimages, noAuto_CSIR_IC, 'r');
hold on 
plot(CSIR_ISARimages, DSA_CSIR_IC, 'b');
hold on
plot(CSIR_ISARimages,Yuan_CSIR_IC,'g');
hold on 
plot(CSIR_ISARimages, modDSA_CSIR_IC,'m');
hold on
plot(CSIR_ISARimages,modYuan_CSIR_IC,'c');
hold off
xlabel('ISAR image');
ylabel('Image Contrast');
legend('No Autofocus','DSA','MSA','Modified DSA', 'Modified MSA');

%IC for Yuan MSA and Haywood DSA
figure
plot(CSIR_ISARimages, noAuto_CSIR_IC, 'r');
hold on
plot(CSIR_ISARimages, DSA_CSIR_IC, 'b');
hold on
plot(CSIR_ISARimages,Yuan_CSIR_IC,'g');
hold off
xlabel('ISAR image');
ylabel('Image Contrast');
legend('No Autofocus','DSA','MSA');


%IC for modMSA and modDSA
figure
plot(CSIR_ISARimages, noAuto_CSIR_IC, 'r');
hold on
plot(CSIR_ISARimages, modDSA_CSIR_IC, 'b');
hold on
plot(CSIR_ISARimages,modYuan_CSIR_IC,'g');
hold off
xlabel('ISAR image');
ylabel('Image Contrast');
legend('No Autofocus','Modified DSA','Modified MSA');
