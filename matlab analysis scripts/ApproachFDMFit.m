clear all
close all

%code to fit FDM model to measured approach curve
%goal is to be able to extract tapping amplitude and tip properties

%% import data and create arrays

%IMPORT HERE
%% Set up the Import Options and import the data
opts = delimitedTextImportOptions("NumVariables", 38);

% Specify range and delimiter
opts.DataLines = [2, Inf];
opts.Delimiter = "\t";

% Specify column names and types
opts.VariableNames = ["Voltage", "NF21", "Voltage1", "NF31", "Voltage2", "NF41", "Voltage3", "Amplitude2Dn", "Voltage4", "Amplitude2Up", "Voltage5", "NF22Dn", "Voltage6", "NF22Up", "Voltage7", "NF32Dn", "Voltage8", "NF32Up", "Voltage9", "NF42Dn", "Voltage10", "NF42Up", "Voltage11", "Amplitude3Dn", "Voltage12", "Amplitude3Up", "Voltage13", "NF23Dn", "Voltage14", "NF23Up", "Voltage15", "NF33Dn", "Voltage16", "NF33Up", "Voltage17", "NF43Dn", "Voltage18", "NF43Up"];
opts.VariableTypes = ["double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double"];

% Specify file level properties
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";

% Import the data
AC = readtable("\\nask.man.ac.uk\home$\Downloads\1550 AC AU.txt", opts);


%% Clear temporary variables
clear opts

%assign and remove NaN values
XData = AC.Voltage4;
XData = XData(~isnan(XData));

%oscillation amplitude
OscAmplitude = AC.Amplitude2Up;
OscAmplitude = OscAmplitude(~isnan(OscAmplitude));
%Signals
O2A = AC.NF22Up;
O2A = O2A(~isnan(O2A));
O3A = AC.NF32Up;
O3A = O3A(~isnan(O3A));
O4A = AC.NF42Up;
O4A = O4A(~isnan(O4A));


%% constants/variables
% E_s = Dielectric_Au(1550); %sample dielectric

E_s = 11.7 +0.1i;

E_t = 1; %tip dielectric (not used)
Tmin = 1e-9; %tip minimum distance
freq = 1; %tip frequency
T = 1/freq; %tip period
g = 0.7*exp(0.06*1i); %g factor

%% normalise data

%detect when the tip goes into contact mode (use figure 1 for this)
n = 300;
OscContact = n;

figure(1)
hold on
plot(XData,OscAmplitude)
xline(XData(OscContact))
title('Contact Point')
ylabel('Tapping Amplitude (mV)')
xlabel('Height (nm)')
hold off
% Shift data to start at contact mode behaviour
XData = XData(OscContact:size(XData));
OscAmplitude = OscAmplitude(OscContact:size(OscAmplitude));
O2A = O2A(OscContact:size(O2A));
O3A = O3A(OscContact:size(O3A));
O4A = O4A(OscContact:size(O4A));


% Scale to value at tip equiblibrium
O2A = O2A/O2A(1);
O3A = (O3A/O3A(1));
O4A = O4A/O4A(1);

%shift hard contact mode to 0nm tip distance
XData = XData - XData(1);

%background shift
%we are assuming the last quarter of the curve is flat here
platau = O3A(size(XData,1)*0.75:size(XData,1));
yshift = mean(platau);

O3A = O3A-yshift;

%% fitting
%3nd Harmonic fitting
[paramsOut, resnorm, residuals, exitflag, output, lambda, jacobian] = FDMfit(XData, O3A, E_s, E_t, Tmin, freq,3,g);
L = paramsOut(1);
R = paramsOut(2);
TAmp = paramsOut(3);

Beta = (E_s - 1)./(E_s +1);



%% 95% confidence interval calculation
%nm conversion
nmParams = paramsOut / 1e-9;
CI = nlparci(nmParams,residuals,'Jacobian',jacobian);

%outputs
Length = nmParams(1)
Length95CI = CI(1,:)

Radius = nmParams(2)
Radius95CI = CI(2,:)

TappingAmplitude = nmParams(3)
TappingAmplitude95CI = CI(3,:)


%plot from there and hope.

%produce fitted curve
for c = 1:size(XData,1)

    Tminplot = XData(c)*1e-9;
    
    O3AFitEqn = @(t) FDMAndreas(Beta,R,Tminplot,TAmp, freq, 3, t, g, L);
    O3AFit(c) = abs(integral(O3AFitEqn, -T/2, +T/2,'arrayvalued', true));

    O3ALowerFitRadiusEqn = @(t) FDMAndreas(Beta,Radius95CI(1)*1e-9,Tminplot,TAmp,freq,3,t,g,L);
    O3ALowerFitRadius(c) = abs(integral(O3ALowerFitRadiusEqn, -T/2, +T/2,'arrayvalued', true));

    O3AHigherFitRadiusEqn = @(t) FDMAndreas(Beta,Radius95CI(2)*1e-9,Tminplot,TAmp,freq,3,t,g,L);
    O3AHigherFitRadius(c) = abs(integral(O3AHigherFitRadiusEqn, -T/2, +T/2,'arrayvalued', true));

    O3ALowerFitTAmpEqn = @(t) FDMAndreas(Beta,R,Tminplot,TappingAmplitude95CI(1)*1e-9, freq, 3, t, g, L);
    O3ALowerFitTAmp(c) = abs(integral(O3ALowerFitTAmpEqn, -T/2, +T/2,'arrayvalued', true));

    O3AHigherFitTAmpEqn = @(t) FDMAndreas(Beta,R,Tminplot,TappingAmplitude95CI(2)*1e-9, freq, 3, t, g, L);
    O3AHigherFitTAmp(c) = abs(integral(O3AHigherFitTAmpEqn, -T/2, +T/2,'arrayvalued', true));

    O3ALowerFitLEqn = @(t) FDMAndreas(Beta,R,Tminplot,TAmp, freq, 3, t, g, Length95CI(1)*1e-9);
    O3ALowerFitL(c) = abs(integral(O3ALowerFitLEqn, -T/2, +T/2,'arrayvalued', true));

    O3AHigherFitLEqn = @(t) FDMAndreas(Beta,R,Tminplot,TAmp, freq, 3, t, g, Length95CI(2)*1e-9);
    O3AHigherFitL(c) = abs(integral(O3AHigherFitLEqn, -T/2, +T/2,'arrayvalued', true));
end

%normalise to 1 and switch row to column 
O3AFit = O3AFit.';
O3AFit =(O3AFit/O3AFit(1));
O3ALowerFitRadius = (O3ALowerFitRadius.'/O3ALowerFitRadius(1));
O3AHigherFitRadius = (O3AHigherFitRadius.'/O3AHigherFitRadius(1));
O3ALowerFitTAmp = (O3ALowerFitTAmp.'/O3ALowerFitTAmp(1));
O3AHigherFitTAmp = (O3AHigherFitTAmp.'/O3AHigherFitTAmp(1));
O3ALowerFitL = (O3ALowerFitL.'/O3ALowerFitL(1));
O3AHigherFitL = (O3AHigherFitL.'/O3AHigherFitL(1));


%% plotting

figure(2)
hold on
title('Approach curve data 3rd Harmonic + fitted approach curve')
ylabel('Signal Amplitude (Norm)')
xlabel('Height (nm)')
plot(XData,O3A);
plot(XData,O3AFit,'lineWidth',3);
legend('3rd Harmonic','fitted FDM 3rd Harmonic')
%plot(XData,O3A);
%plot(XData,O4A);
hold off

figure(3)
hold on
title('Radius 95% CI Bounds')
plot(XData,O3ALowerFitRadius)
plot(XData,O3AFit)
plot(XData,O3AHigherFitRadius)
ylabel('Signal Amplitude (Norm)')
xlabel('Height (nm)')
legend('Lower Radius CI','Fitted Radius','Higher Radius CI')
hold off

figure(4)
hold on
title('Tapping Amplitude 95% CI Bounds')
plot(XData,O3ALowerFitTAmp)
plot(XData,O3AFit)
plot(XData,O3AHigherFitTAmp)
ylabel('Signal Amplitude (Norm)')
xlabel('Height (nm)')
legend('Lower TappingAmp CI','Fitted TappingAmp','Higher TappingAmp CI')
hold off

figure(5)
hold on
title('Length 95% CI Bounds')
plot(XData,O3ALowerFitL)
plot(XData,O3AFit)
plot(XData,O3AHigherFitL)
xlabel('Signal Amplitude (Norm)')
ylabel('Height (nm)')
legend('Lower Length CI','Fitted Length','Higher Length CI')
hold off
