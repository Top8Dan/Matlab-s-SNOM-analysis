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
n = 4;

%plotting parameters
Resolution = 100; %how many datapoints on graph
freq0  = 1100; %start of spectra (wavenumber (cm))
freqMax = 1600; %end of spectra (wavenumber (cm))

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



%sample Dielectic
Es_1 = 1; %air
% Es_2 = DrudeLorentzDielectricEH(freqX,1.1,1250,100,1430,110); %Surface mode;
Es_2 = LorentzOscillator(freqX,50,5e28,1550,10e12);
Es_3 = DrudeFunction(freqX,17.4,3.72e18,150,0.6); %Bi2Te3 - https://pubs.acs.org/doi/10.1021/acs.nanolett.8b03008
Es_4 = zeros(Resolution,1) + Es_Ref;
%sample Thickness
d2 = 7e-9;
d3 = 50e-9;


%Ref Dielectric
Es_Ref1 = 1; %air
Es_Ref2 = DrudeLorentzDielectricEH(freqX,5,1150,90,1250,80); %SiO2;
Es_Ref3 = zeros(Resolution,1) + Es_Ref;
Es_Ref4 = zeros(Resolution,1) + Es_Ref;
%Ref Thickness
dRef2 = 0;
dRef3 = 1;



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




%% plotting

figure(1)
hold on
plot(freqX,real(Es_2),'linewidth',2,'Color','b')
plot(freqX,imag(Es_2),'linewidth',2,'Color','r')
plot(freqX,real(Es_Ref2),'--','linewidth',2,'Color','b')
plot(freqX,imag(Es_Ref2),'--','linewidth',2,'Color','r')
xlabel('Wavenumber (cm-1)')
ylabel('Dielectric')
legend('Real (SS)','Imag (SS)','Real (SiO2)','Imag (SiO2)')
title('Input dielectric Es_2')
ax = gca;
ax.FontSize = 14;
xlim([freq0 freqMax])

figure(2)
hold on
plot(freqX,real(Es_3),'linewidth',2,'Color','b')
plot(freqX,imag(Es_3),'linewidth',2,'Color','r')
xlabel('Wavenumber (cm-1)')
ylabel('dielectric')
legend('Real (Bi2Te3)','Imag (Bi2Te3)')
title('Input dielectric Es_3')
ax = gca;
ax.FontSize = 14;
xlim([freq0 freqMax])




figure(5)
hold on
plot(freqX,real(NormS0),'linewidth',2,'Color','b')
ylabel('Amplitude (AU)')
yyaxis right
plot(freqX,angle(NormS0),'linewidth',2,'Color','r')
title('FDM Output')
ylabel('Phase (\phi)')
xlabel('Wavenumber (cm^{-1})')
ax = gca;
ax.FontSize = 14;
xlim([freq0 freqMax])
hold off
