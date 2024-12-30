%Combination of amplitude 3/4.
close all
clear all

SaveLocation = 'C:\Users\k68558dj\OneDrive - The University of Manchester\DataFiles\TI NW SN123\Bi2Te3 NW July 23 (flood doped)\Undoped\FDM\Source B';

Harmonic3Fig = openfig("C:\Users\k68558dj\OneDrive - The University of Manchester\DataFiles\TI NW SN123\Bi2Te3 NW July 23 (flood doped)\Undoped\FDM\Source B\3rd Harmonic\SourceBAmpPhaseMulti.fig");

Harmonic4Fig = openfig("C:\Users\k68558dj\OneDrive - The University of Manchester\DataFiles\TI NW SN123\Bi2Te3 NW July 23 (flood doped)\Undoped\FDM\Source B\4th Harmonic\SourceBAmpPhaseMulti.fig");

Data3 = findobj(Harmonic3Fig,'-property','YData');
Data4 = findobj(Harmonic4Fig,'-property','YData');


ax = gca;

ax.XLim;

AmpData3X = Data3(2).XData;
AmpData3Y = Data3(2).YData;
PhaseData3Y = Data3(1).YData;


AmpData4X = Data4(2).XData;
AmpData4Y = Data4(2).YData;
PhaseData4Y = Data4(1).YData;

Ratio43Amp = AmpData4Y ./ AmpData3Y;
Ratio43Phase = PhaseData4Y - PhaseData3Y;




figure(3)
hold on
plot(AmpData3X,Ratio43Amp,'linewidth',2,'Color','b')
ylabel('Amplitude (AU)')
yyaxis right
plot(AmpData3X,Ratio43Phase,'linewidth',2,'Color','r')
xlim(ax.XLim)
title('Ratio S4/S3')
ylabel('Phase (\phi)')
xlabel('Wavenumber (cm^{-1})')
ax = gca;
ax.FontSize = 14;
%xlim([900 1400])
hold off


figure(4)
hold on
plot(AmpData3X,AmpData3Y,'b-','linewidth',2)
ylabel('Amplitude (AU)')
plot(AmpData4X,AmpData4Y,'b--','linewidth',2)
xlim(ax.XLim)
%xlim([900 1400])
xlabel('Wavenumber (cm^{-1})')
legend('3rd Harmonic','4th Harmonic')
ax = gca;
ax.FontSize = 14;

figure(5)
hold on
plot(AmpData3X,PhaseData3Y,'r-','linewidth',2)
ylabel('Phase (\phi)')
plot(AmpData4X,PhaseData4Y,'r--','linewidth',2)
xlim(ax.XLim)
%xlim([900 1400])
xlabel('Wavenumber (cm^{-1})')
legend('3rd Harmonic','4th Harmonic')
ax = gca;
ax.FontSize = 14;

saveas(figure(3),strcat(SaveLocation, '/Ratio43.fig'))
saveas(figure(4),strcat(SaveLocation, '/AmpMulti.fig'))
saveas(figure(5),strcat(SaveLocation, '/PhaseMulti.fig'))
