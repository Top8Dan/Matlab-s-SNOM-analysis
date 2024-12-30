%Combination of amplitude and phase.
close all
clear all


AmplitudeFig3 = openfig("C:\Users\k68558dj\OneDrive - The University of Manchester\DataFiles\Lars-Visit Jul2024\2024-07-19 298 9K spectra\Analysis NFSpectra\3rd Harmonic\SourceAAvg Amp.fig");
AmplitudeFig4 = openfig("C:\Users\k68558dj\OneDrive - The University of Manchester\DataFiles\Lars-Visit Jul2024\2024-07-19 298 9K spectra\Analysis NFSpectra\4th Harmonic\SourceAAvg Amp.fig");

PhaseFig3 = openfig("C:\Users\k68558dj\OneDrive - The University of Manchester\DataFiles\Lars-Visit Jul2024\2024-07-19 298 9K spectra\Analysis NFSpectra\3rd Harmonic\SourceAAvg Angle.fig");
PhaseFig4 = openfig("C:\Users\k68558dj\OneDrive - The University of Manchester\DataFiles\Lars-Visit Jul2024\2024-07-19 298 9K spectra\Analysis NFSpectra\4th Harmonic\SourceAAvg Angle.fig");

AmpData3 = findobj(AmplitudeFig3,'-property','YData');
AmpData4 = findobj(AmplitudeFig4,'-property','YData');

PhaseData3 = findobj(PhaseFig3,'-property','YData');
PhaseData4 = findobj(PhaseFig4,'-property','YData');

ax = gca;

ax.XLim;

AmpDataX3 = AmpData3.XData;
AmpDataY3 = AmpData3.YData;
AmpDataX4 = AmpData4.XData;
AmpDataY4 = AmpData4.YData;

PhaseDataX3 = PhaseData3.XData;
PhaseDataY3 = PhaseData3.YData;
PhaseDataX4 = PhaseData4.XData;
PhaseDataY4 = PhaseData4.YData;

AbsorptionSpectra = AmpDataY3.*sin(PhaseDataY3);


figure(6)
hold on
plot(AmpDataX3,AmpDataY3,'linewidth',2,'Color','b')
ylabel('Amplitude')
yyaxis right
plot(PhaseDataX3,PhaseDataY3,'linewidth',2,'Color','r')
ylabel('Phase (\phi)')
xlim(ax.XLim)
xlabel('Wavenumber (cm^{-1})')
title('Si Norm to Au 3rd Harmonic')
hold off

figure(7)
hold on
plot(AmpDataX4,AmpDataY4,'linewidth',2,'Color','b')
ylabel('Amplitude')
yyaxis right
plot(PhaseDataX4,PhaseDataY4,'linewidth',2,'Color','r')
ylabel('Phase (\phi)')
xlim(ax.XLim)
xlabel('Wavenumber (cm^{-1})')
title('Si Norm to Au 4th Harmonic')
hold off


figure(8)
hold on
plot(AmpDataX3,AmpDataY3,'linewidth',2,'Color','b')
plot(AmpDataX4,AmpDataY4,'linewidth',2,'Color','b','LineStyle','--')
ylabel('Amplitude (arb.)')
yyaxis right
plot(PhaseDataX3,PhaseDataY3,'linewidth',2,'Color','r')
plot(PhaseDataX4,PhaseDataY4,'linewidth',2,'Color','r','LineStyle','--')
ylabel('Phase (\phi)')
xlim(ax.XLim)
xlabel('Wavenumber (cm^{-1})')
title('Bi2Te3 Norm to Si')
legend('Amp 3rd','Amp 4th','Phase 3rd','Phase 4th')
set(gca,'FontSize',14);
set(gca,'YColor','r');
hold off



% SaveFolder = 'C:\Users\k68558dj\OneDrive - The University of Manchester\NWPaper\Spectra (fig1)\Cylin NW\Source B';
% saveas(figure(3),strcat(SaveFolder, '\AmpPhase.fig'))
% saveas(figure(4),strcat(SaveFolder, '\AbsorptionSpectra.fig'))