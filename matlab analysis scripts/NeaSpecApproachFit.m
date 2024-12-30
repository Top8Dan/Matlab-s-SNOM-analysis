close all
clear all

%% Set up the Import Options and import the data
opts = delimitedTextImportOptions("NumVariables", 30);

% Specify range and delimiter
opts.DataLines = [29, Inf];
opts.Delimiter = "\t";

% Specify column names and types
opts.VariableNames = ["Row", "Column", "Run", "Depth", "Z", "M0A", "M0P", "M1A", "M1P", "M2A", "M2P", "M3A", "M3P", "M4A", "M4P", "M5A", "M5P", "O0A", "O0P", "O1A", "O1P", "O2A", "O2P", "O3A", "O3P", "O4A", "O4P", "O5A", "O5P", "VarName30"];
opts.VariableTypes = ["double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "string"];

% Specify file level properties
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";

% Specify variable properties
opts = setvaropts(opts, "VarName30", "WhitespaceRule", "preserve");
opts = setvaropts(opts, "VarName30", "EmptyFieldRule", "auto");
opts = setvaropts(opts, ["Row", "Column", "Run", "Depth", "Z", "M0A", "M0P", "M1A", "M1P", "M2A", "M2P", "M3A", "M3P", "M4A", "M4P", "M5A", "M5P", "O0A", "O0P", "O1A", "O1P", "O2A", "O2P", "O3A", "O3P", "O4A", "O4P", "O5A", "O5P"], "ThousandsSeparator", ",");

% Import the data

ApproachSample = readtable("C:\Users\k68558dj\OneDrive - The University of Manchester\DataFiles\TI NW SN123\Bi2Te3 NW July 23\Dan_BiTe_dopedNWs\2023-07-19 1852-Undoped NW\2023-07-19 135258 AC Si_SourceE\2023-07-19 135258 AC Si_SourceE.txt",opts);

clear opts
%% constants/variables
% E_s = Dielectric_Au(1550); %sample dielectric

E_s = 11.7 +0.1i;
Tamp = 40e-9;

E_t = 1; %tip dielectric (not used)
Tmin = 1e-9; %tip minimum distance
freq = 1; %tip frequency
T = 1/freq; %tip period
g = 0.7*exp(0.06*1i); %g factor
%% conversion

%convert Z M to nm

ApproachSample.Z = ApproachSample.Z * 1e9;


%% handling


%% identify the contact point
ContactArea = 30; %0 to N first points
FreeSpaceArea = 700; %N to max points
%fit of contact area
ContactFit = fit(ApproachSample.Z(1:ContactArea),ApproachSample.M1A(1:ContactArea),'poly1');
%fit of free space area
FreeSpaceFit = fit(ApproachSample.Z(FreeSpaceArea:size(ApproachSample.Z),1),ApproachSample.M1A(FreeSpaceArea:size(ApproachSample.Z),1),'poly1');
ContactArray = ContactFit(ApproachSample.Z) - FreeSpaceFit(ApproachSample.Z);

ContactPoint = 1;
while ContactArray(ContactPoint)<0
    ContactPoint = ContactPoint+1;
end

%normalize
O2ANorm = ApproachSample.O2A./ApproachSample.O2A(ContactPoint);
O3ANorm = ApproachSample.O3A./ApproachSample.O3A(ContactPoint);
O4ANorm = ApproachSample.O4A./ApproachSample.O4A(ContactPoint);
O5ANorm = ApproachSample.O5A./ApproachSample.O5A(ContactPoint);
XData = ApproachSample.Z;

%% window

%window of data for fitting

O2DataWindow = 300;
O3DataWindow = 300;
O4DataWindow = 300;
O5DataWindow = 300;

O2ANormWin = O2ANorm(ContactPoint  :O2DataWindow);
XData2 = XData(ContactPoint :O2DataWindow) - XData(ContactPoint);
O3ANormWin = O3ANorm(ContactPoint:O3DataWindow);
XData3 = XData(ContactPoint :O3DataWindow) - XData(ContactPoint);
O4ANormWin = O4ANorm(ContactPoint :O4DataWindow);
XData4 = XData(ContactPoint :O4DataWindow) - XData(ContactPoint);
O5ANormWin = O5ANorm(ContactPoint :O5DataWindow);
XData5 = XData(ContactPoint :O5DataWindow) - XData(ContactPoint);



%% fitting

%3nd Harmonic fitting
[paramsOut, resnorm, residuals, exitflag, output, lambda, jacobian] = FDMfit(XData3, O3ANormWin, E_s, E_t, Tmin,Tamp, freq,3,g);
L = paramsOut(1);
R = paramsOut(2);

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
% 
% TappingAmplitude = nmParams(3)
% TappingAmplitude95CI = CI(3,:)


%plot from there and hope.

%produce fitted curve
for c = 1:size(XData3,1)

    Tminplot = XData3(c)*1e-9;
    
    O3AFitEqn = @(t) FDMHauer(Beta,R,Tminplot,Tamp, freq, 3, t, g, L);
    O3AFit(c) = abs(integral(O3AFitEqn, -T/2, +T/2,'arrayvalued', true));

    O3ALowerFitRadiusEqn = @(t) FDMHauer(Beta,Radius95CI(1)*1e-9,Tminplot,Tamp,freq,3,t,g,L);
    O3ALowerFitRadius(c) = abs(integral(O3ALowerFitRadiusEqn, -T/2, +T/2,'arrayvalued', true));

    O3AHigherFitRadiusEqn = @(t) FDMHauer(Beta,Radius95CI(2)*1e-9,Tminplot,Tamp,freq,3,t,g,L);
    O3AHigherFitRadius(c) = abs(integral(O3AHigherFitRadiusEqn, -T/2, +T/2,'arrayvalued', true));

%     O3ALowerFitTAmpEqn = @(t) FDMHauer(Beta,R,Tminplot,TappingAmplitude95CI(1)*1e-9, freq, 3, t, g, L);
%     O3ALowerFitTAmp(c) = abs(integral(O3ALowerFitTAmpEqn, -T/2, +T/2,'arrayvalued', true));
% 
%     O3AHigherFitTAmpEqn = @(t) FDMHauer(Beta,R,Tminplot,TappingAmplitude95CI(2)*1e-9, freq, 3, t, g, L);
%     O3AHigherFitTAmp(c) = abs(integral(O3AHigherFitTAmpEqn, -T/2, +T/2,'arrayvalued', true));

    O3ALowerFitLEqn = @(t) FDMHauer(Beta,R,Tminplot,Tamp, freq, 3, t, g, Length95CI(1)*1e-9);
    O3ALowerFitL(c) = abs(integral(O3ALowerFitLEqn, -T/2, +T/2,'arrayvalued', true));

    O3AHigherFitLEqn = @(t) FDMHauer(Beta,R,Tminplot,Tamp, freq, 3, t, g, Length95CI(2)*1e-9);
    O3AHigherFitL(c) = abs(integral(O3AHigherFitLEqn, -T/2, +T/2,'arrayvalued', true));
end

%normalise to 1 and switch row to column 
O3AFit = O3AFit.';
O3AFit =(O3AFit/O3AFit(1));
O3ALowerFitRadius = (O3ALowerFitRadius.'/O3ALowerFitRadius(1));
O3AHigherFitRadius = (O3AHigherFitRadius.'/O3AHigherFitRadius(1));
% O3ALowerFitTAmp = (O3ALowerFitTAmp.'/O3ALowerFitTAmp(1));
% O3AHigherFitTAmp = (O3AHigherFitTAmp.'/O3AHigherFitTAmp(1));
O3ALowerFitL = (O3ALowerFitL.'/O3ALowerFitL(1));
O3AHigherFitL = (O3AHigherFitL.'/O3AHigherFitL(1));


%% plotting

figure(2)
hold on
title('Approach curve data 3rd Harmonic + fitted approach curve')
ylabel('Signal Amplitude (Norm)')
xlabel('Height (nm)')
plot(XData3,O3ANormWin);
plot(XData3,O3AFit,'lineWidth',3);
legend('3rd Harmonic','fitted FDM 3rd Harmonic')
%plot(XData,O3A);
%plot(XData,O4A);
hold off

figure(3)
hold on
title('Radius 95% CI Bounds')
plot(XData3,O3ALowerFitRadius)
plot(XData3,O3AFit)
plot(XData3,O3AHigherFitRadius)
ylabel('Signal Amplitude (Norm)')
xlabel('Height (nm)')
legend('Lower Radius CI','Fitted Radius','Higher Radius CI')
hold off

% figure(4)
% hold on
% title('Tapping Amplitude 95% CI Bounds')
% plot(XData3,O3ALowerFitTAmp)
% plot(XData3,O3AFit)
% plot(XData3,O3AHigherFitTAmp)
% ylabel('Signal Amplitude (Norm)')
% xlabel('Height (nm)')
% legend('Lower TappingAmp CI','Fitted TappingAmp','Higher TappingAmp CI')
% hold off

figure(5)
hold on
title('Length 95% CI Bounds')
plot(XData3,O3ALowerFitL)
plot(XData3,O3AFit)
plot(XData3,O3AHigherFitL)
xlabel('Signal Amplitude (Norm)')
ylabel('Height (nm)')
legend('Lower Length CI','Fitted Length','Higher Length CI')
hold off