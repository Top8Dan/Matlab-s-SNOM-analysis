close all
clear all

% eqns from Vibrational modes in amorphous silicon dioxide Marta Klanjs\ek Gunde 

%% parameters
Radius = 35e-9; %tip radius
Tmin = 10e-9; %tip minimum distance
Tamp = 30e-9; %tip amplitude position
freq = 1; %tip frequency
T = 1/freq; %tip period

%finite Dipole parameters
g = 0.7*exp(0.06*1i);
L = 300e-9;
%demodulation number
n = 3;

%plotting parameters
Resolution = 100; %how many datapoints on graph
freq0  = 800; %start of spectra (wavenumber (cm))
freqMax = 2000; %end of spectra (wavenumber (cm))

%make zeros
Sn = zeros(Resolution,1);
SnAu = zeros(Resolution,1);

%frequency array (cm-1)
freqX = linspace(freq0,freqMax,Resolution);



%layer parameters
Es_Ref = 11.7 +0i; % Si

%sample Dielectic
Es_1 = 1; %air
Es_2 = Es_Ref; %DrudeLorentzDielectricEH(freqX); %SiO2
Es_3 = LorentzOscillator(freqX,11.7,1e29,1400,5e13); 
Es_4 = Es_Ref;
%sample Thickness
d2 = 5e-9;
d3 = 5e-9;


%Ref Dielectric
Es_Ref1 = 1; %air
Es_Ref2 = 3.9 - 0i;
Es_Ref3 = Es_Ref;
Es_Ref4 = Es_Ref;
%Ref Thickness
dRef2 = 0;
dRef3 = 0;


%% handling


%calculate FDM output for the input dielectric
W_0 = 1.31*Radius*L/(L+2*Radius);
W_1 = 0.5*Radius;
iteration = 0;
Ref = ML_n_3layers(Tmin, Tamp, freq, Radius, L, g, W_0, W_1, Es_Ref1, Es_Ref2, Es_Ref3, Es_Ref4, dRef2, dRef3,n);
for z = 1:Resolution
    S0(z) = ML_n_3layers(Tmin, Tamp, freq, Radius, L, g, W_0, W_1, Es_1, Es_2, Es_3(z), Es_4, 0, d3, n);
    S1(z) = ML_n_3layers(Tmin, Tamp, freq, Radius, L, g, W_0, W_1, Es_1, Es_2, Es_3(z), Es_4, 1e-9, d3, n);
    S2(z) = ML_n_3layers(Tmin, Tamp, freq, Radius, L, g, W_0, W_1, Es_1, Es_2, Es_3(z), Es_4, 3e-9, d3, n);
    S3(z) = ML_n_3layers(Tmin, Tamp, freq, Radius, L, g, W_0, W_1, Es_1, Es_2, Es_3(z), Es_4, 5e-9, d3, n);
    S4(z) = ML_n_3layers(Tmin, Tamp, freq, Radius, L, g, W_0, W_1, Es_1, Es_2, Es_3(z), Es_4, 10e-9, d3,  n);

    SRef(z) = Ref;
    iteration = iteration+1
end

NormS0 = S0./SRef;
NormS1 = S1./SRef;
NormS3 = S2./SRef;
NormS6 = S3./SRef;
NormS10 = S4./SRef;

%% plotting

Es = Es_3;

figure(1)
hold on
plot(freqX,real(Es),'linewidth',2,'Color','b')
ylabel('Real(\epsilon)')
yyaxis right
ylabel('Image(\epsilon)')
plot(freqX,imag(Es),'linewidth',2,'Color','r')
xlabel('Wavenumber (cm-1)')

legend('Re(\epsilon)','Im(\epsilon)')
title('Input dielectric')

% figure(2)
% hold on
% plot(freqX,real(NormS))
% yyaxis right
% plot(freqX,imag(NormS))
% xlabel('Wavenumber (cm-1)')
% ylabel('Amp Phase')
% legend('Real','Imag')
% title('Output Multilayer FDM')


figure(3)
hold on
plot(freqX,real(NormS0))
plot(freqX,real(NormS1))
plot(freqX,real(NormS3))
plot(freqX,real(NormS6))
plot(freqX,real(NormS10))
xlabel('Wavenumber (cm^{-1})')
ylabel('Amplitude (arb.)')
legend('0nm','1nm','3nm','5nm','10nm')
title('Varying Thickness of Si over Oscillator / Si')


figure(4)
hold on
plot(freqX,angle(NormS0))
plot(freqX,angle(NormS1))
plot(freqX,angle(NormS3))
plot(freqX,angle(NormS6))
plot(freqX,angle(NormS10))
xlabel('Wavenumber (cm^{-1})')
ylabel('phase (\phi)')
legend('0nm','1nm','3nm','5nm','10nm')
title('Varying Thickness of Si over Oscillator / Si')
