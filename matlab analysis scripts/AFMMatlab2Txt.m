AFMFig = openfig("C:\Users\k68558dj\OneDrive - The University of Manchester\DataFiles\NPL Measurements\June 2022\nanowire 2\tip2\Images\1140\NormSubstrate\AFM.fig");

AFMData = findobj(AFMFig);

ImageArray = AFMData(5);

Array = ImageArray.CData;

saveas