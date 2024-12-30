%Combination of Phase 3/4.
close all
clear all


Phase3Fig = openfig("C:\Users\k68558dj\OneDrive - The University of Manchester\NWPaper\Spectra (fig1)\Cylin NW\Wider Window\Source E\3rd Harmonic\SourceEAvg Angle.fig");

Phase4Fig = openfig("C:\Users\k68558dj\OneDrive - The University of Manchester\NWPaper\Spectra (fig1)\Cylin NW\Wider Window\Source E\4th Harmonic\SourceEAvg Angle.fig");

Phase3Data = findobj(Phase3Fig,'-property','YData');
Phase4Data = findobj(Phase4Fig,'-property','YData');

ax = gca;

ax.XLim;

PhaseData3X = (Phase3Data.XData);
PhaseData4X = (Phase4Data.XData);

PhaseData3Y = (Phase3Data.YData);
PhaseData4Y = (Phase4Data.YData);

% PhaseData3Y = unwrap(PhaseData3Y);
% PhaseData4Y = unwrap(PhaseData4Y);
% PhaseData3Y = PhaseData3Y - 22;
% PhaseData4Y = PhaseData4Y - 34;

figure(3)
hold on
plot(PhaseData3X,PhaseData3Y,'r-')
ylabel('Phase (\phi)')
plot(PhaseData4X,PhaseData4Y,'r--')
xlim(ax.XLim);
%xlim([900 1400])
xlabel('Wavenumber (cm^{-1})')
legend('3rd Harmonic','4th Harmonic')

saveas(figure(3),strcat('C:\Users\k68558dj\OneDrive - The University of Manchester\NWPaper\Spectra (fig1)\Cylin NW\Wider Window\Source E' , '/34PhaseMulti.fig'))