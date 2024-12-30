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
freqMax = 2100; %end of spectra (wavenumber (cm))

%make zeros
Sn = zeros(Resolution,1);
SnAu = zeros(Resolution,1);

%frequency array (cm-1)
freqX = linspace(freq0,freqMax,Resolution);

%%SiO2 values Marta Klanjs\ek Gunde
% Es_inf = 1.3; %high frequency permittivity
% OmegaLv = 1073;
% yLv = 90;
% OmegaTv = 1165;
% yTv = 80;
% DrudeLorentzDielectricEH(freqX,1.3,1072,90,1130,80); %SiO2;

%%Bi2Te3 values 
% %Bi2Te3 - https://pubs.acs.org/doi/10.1021/acs.nanolett.8b03008
% DrudeFunction(freqX,30,5e19,150,0.6);
% DiracDielectric(freqX,1.25e13,5e5); %dirac carriers on surface

%layer parameters
Es_Ref = 11.7 +0i; % Si

var = [5e11,1e12,5e12,1e13,5e13,1e14];

for i = 1:size(var,2)

%sample Dielectic
Es_1 = 1; %air
% Es_2 = DrudeLorentzDielectricEH(freqX,1.1,1250,100,1430,110); %Surface mode;
Es_2 = LorentzOscillator(freqX,10,50e27,1250,var(i));
Es_3 = DrudeFunction(freqX,17.4,3.72e18,150,0.6); %Bi2Te3 - https://pubs.acs.org/doi/10.1021/acs.nanolett.8b03008
Es_4 = zeros(Resolution,1) + Es_Ref;
%sample Thickness
d2 = 7e-9;
d3 = 50e-9;


%Ref Dielectric
Es_Ref1 = 1; %air
Es_Ref2 = DrudeLorentzDielectricEH(freqX,1.3,1072,90,1130,80); %SiO2;
Es_Ref3 = DrudeFunction(freqX,17.4,3.72e18,150,0.6); %Bi2Te3 - https://pubs.acs.org/doi/10.1021/acs.nanolett.8b03008
Es_Ref4 = zeros(Resolution,1) + Es_Ref;
%Ref Thickness
dRef2 = 5e-9;
dRef3 = 0;



%% handling


%calculate FDM output for the input dielectric
W_0 = 1.31*Radius*L/(L+2*Radius);
W_1 = 0.5*Radius;
iteration = 0;

for z = 1:Resolution

    S0(z) = ML_n_3layers(Tmin, Tamp, freq, Radius, L, g, W_0, W_1, Es_1, Es_2(z), Es_3(z), Es_4(z), d2, d3,n);

    SRef(z) = ML_n_3layers(Tmin, Tamp, freq, Radius, L, g, W_0, W_1, Es_Ref1, Es_Ref2(z), Es_Ref3(z), Es_Ref4(z), dRef2, dRef3,n);

    iteration = iteration+1
end

NormS0 = S0./SRef;

NormArray{i} = NormS0;
Es_2Array{i} = Es_2;
end


% Create legend labels from the double array
legendLabels = arrayfun(@(x) sprintf('%.1e', x), var, 'UniformOutput', false);

%% plotting

figure(1)
hold on
for i = 1:size(var,2)
plot(freqX,real(Es_2Array{i}),'linewidth',2)
end
xlabel('Wavenumber (cm-1)')
ylabel('Real Dielectric')
legend(legendLabels)
title('Input dielectric Es_2')
ax = gca;
ax.FontSize = 14;
xlim([freq0 freqMax])
hold off

figure(2)
hold on
for i = 1:size(var,2)
plot(freqX,imag(Es_2Array{i}),'linewidth',2)
end
xlabel('Wavenumber (cm-1)')
ylabel('Imag Dielectric')
legend(legendLabels)
title('Input dielectric Es_2')
ax = gca;
ax.FontSize = 14;
xlim([freq0 freqMax])
hold off


figure(5)
hold on
for i = 1:size(var,2)
plot(freqX,real(NormArray{i}),'linewidth',2)
end
ylabel('Amplitude (AU)')
title('FDM Output (Changing: gamma)')
xlabel('Wavenumber (cm^{-1})')
ax = gca;
ax.FontSize = 14;
xlim([freq0 freqMax])
legend(legendLabels)
hold off

figure(6)
hold on
for i = 1:size(var,2)
plot(freqX,angle(NormArray{i}),'linewidth',2)
end
ylabel('Phase (\phi)')
title('FDM Output (Changing: gamma)')
xlabel('Wavenumber (cm^{-1})')
ax = gca;
ax.FontSize = 14;
xlim([freq0 freqMax])
legend(legendLabels)
hold off