%Combination of amplitude 3/4.
close all
clear all


Amplitude3Fig = openfig("C:\Users\k68558dj\OneDrive - The University of Manchester\NWPaper\Spectra (fig1)\Cylin NW\Wider Window\Source A\3rd Harmonic\SourceAAvg Amp.fig");

Amplitude4Fig = openfig("C:\Users\k68558dj\OneDrive - The University of Manchester\NWPaper\Spectra (fig1)\Cylin NW\Wider Window\Source A\4th Harmonic\SourceAAvg Amp.fig");

AmpData = findobj(Amplitude3Fig,'-property','YData');
PhaseData = findobj(Amplitude4Fig,'-property','YData');

ax = gca;

ax.XLim;

AmpData3X = AmpData.XData;
AmpData3Y = AmpData.YData;

AmpData4X = PhaseData.XData;
AmpData4Y = PhaseData.YData;


figure(3)
hold on
plot(AmpData3X,AmpData3Y,'b-')
ylabel('Amplitude')
plot(AmpData4X,AmpData4Y,'b--')
xlim(ax.XLim)
%xlim([900 1400])
xlabel('Wavenumber (cm^{-1})')
legend('3rd Harmonic','4th Harmonic')

saveas(figure(3),strcat('C:\Users\k68558dj\OneDrive - The University of Manchester\NWPaper\Spectra (fig1)\Cylin NW\Wider Window\Source A' , '/34AmpMulti.fig'))