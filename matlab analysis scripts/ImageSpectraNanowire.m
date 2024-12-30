close all
imageFolder = 'C:\Users\k68558dj\Dropbox (The University of Manchester)\NPL\NPL measurements\nanowire 2\tip2\Images';

N1140 = openfig(strcat(imageFolder, '\1140\NormSubstrate\NF3.fig'));
N1160 = openfig(strcat(imageFolder, '\1160\NormSubstrate\NF3.fig'));
N1250 = openfig(strcat(imageFolder, '\1250\NormSubstrate\NF3.fig'));
N1255 = openfig(strcat(imageFolder, '\1255\NormSubstrate\NF3.fig'));
N1260 = openfig(strcat(imageFolder, '\1260\NormSubstrate\NF3.fig'));
N1270 = openfig(strcat(imageFolder, '\1270\NormSubstrate\NF3.fig'));
N1290 = openfig(strcat(imageFolder, '\1290\NormSubstrate\NF3.fig'));
N1320 = openfig(strcat(imageFolder, '\1320\NormSubstrate\NF3.fig'));
N1340 = openfig(strcat(imageFolder, '\1340\NormSubstrate\NF3.fig'));
N1360 = openfig(strcat(imageFolder, '\1360\NormSubstrate\NF3.fig'));


N1140 = N1140.Children;
N1160 = N1160.Children;
N1250 = N1250.Children;
N1255 = N1255.Children;
N1260 = N1260.Children;
N1270 = N1270.Children;
N1290 = N1290.Children;
N1320 = N1320.Children;
N1340 = N1340.Children;
N1360 = N1360.Children;

N1140 = N1140.Children;
N1160 = N1160.Children;
N1250 = N1250.Children;
N1255 = N1255.Children;
N1260 = N1260.Children;
N1270 = N1270.Children;
N1290 = N1290.Children;
N1320 = N1320.Children;
N1340 = N1340.Children;
N1360 = N1360.Children;

N1140 = N1140(4);
N1160 = N1160(4);
N1250 = N1250(4);
N1255 = N1255(4);
N1260 = N1260(4);
N1270 = N1270(4);
N1290 = N1290(4);
N1320 = N1320(4);
N1340 = N1340(4);
N1360 = N1360(4);

N1140 = N1140.Children;
N1160 = N1160.Children;
N1250 = N1250.Children;
N1255 = N1255.Children;
N1260 = N1260.Children;
N1270 = N1270.Children;
N1290 = N1290.Children;
N1320 = N1320.Children;
N1340 = N1340.Children;
N1360 = N1360.Children;


%calculate nanowire amplitude
XPos1140 = 192;
YPos1140 = 220;

XPos1160 = 200;
YPos1160 = 220;

XPos1250 = 200;
YPos1250 = 160;

XPos1255 = 216;
YPos1255 = 190;

XPos1260 = 150;
YPos1260 = 204;

XPos1270 = 160;
YPos1270 = 210;

XPos1290 = 150;
YPos1290 = 200;

XPos1320 = 180;
YPos1320 = 210;

XPos1340 = 200;
YPos1340 = 220;

XPos1360 = 200;
YPos1360 = 200;

YArray = [N1140.CData(XPos1140,YPos1140),N1160.CData(XPos1160,YPos1160),N1250.CData(XPos1250,YPos1250),...
N1255.CData(XPos1255,YPos1255),N1260.CData(XPos1260,YPos1260),N1270.CData(XPos1270,YPos1270),N1290.CData(XPos1290,YPos1290),...
N1320.CData(XPos1320,YPos1320),N1340.CData(XPos1340,YPos1340),N1360.CData(XPos1360,YPos1360)];

XArray = [1140,1160,1250,1255,1260,1270,1290,1320,1340,1360];

figure(11)
hold on
plot(XArray,YArray, 'LineWidth',2)
plot(XArray,YArray, 'x')
title('Average Amplitude')
xlabel('wavenumber (cm^{-1})')
ylabel('Amplitude average (A.U)')
hold off
