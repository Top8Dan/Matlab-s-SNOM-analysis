close all
imageFolder = 'C:\Users\k68558dj\Dropbox (The University of Manchester)\NPL\anasys NPL-DJ Manchester\ThinFilmImages';

N1280 = openfig(strcat(imageFolder, '\1280\Norm\NF3.fig'));
N1300 = openfig(strcat(imageFolder, '\1300\Norm\NF3.fig'));
N1320 = openfig(strcat(imageFolder, '\1320\Norm\NF3.fig'));
N1520 = openfig(strcat(imageFolder, '\1520\Norm\NF3.fig'));
N1530 = openfig(strcat(imageFolder, '\1530\Norm\NF3.fig'));
N1540 = openfig(strcat(imageFolder, '\1540\Norm\NF3.fig'));
N1550 = openfig(strcat(imageFolder, '\1550\Norm\NF3.fig'));
N1570 = openfig(strcat(imageFolder, '\1570\Norm\NF3.fig'));
N1590 = openfig(strcat(imageFolder, '\1590\Norm\NF3.fig'));


N1280 = N1280.Children;
N1300 = N1300.Children;
N1320 = N1320.Children;
N1520 = N1520.Children;
N1530 = N1530.Children;
N1540 = N1540.Children;
N1550 = N1550.Children;
N1570 = N1570.Children;
N1590 = N1590.Children;


N1280 = N1280.Children;
N1300 = N1300.Children;
N1320 = N1320.Children;
N1520 = N1520.Children;
N1530 = N1530.Children;
N1540 = N1540.Children;
N1550 = N1550.Children;
N1570 = N1570.Children;
N1590 = N1590.Children;


N1280 = N1280(4);
N1300 = N1300(4);
N1320 = N1320(4);
N1520 = N1520(4);
N1530 = N1530(4);
N1540 = N1540(4);
N1550 = N1550(4);
N1570 = N1570(4);
N1590 = N1590(4);

N1280 = N1280.Children;
N1300 = N1300.Children;
N1320 = N1320.Children;
N1520 = N1520.Children;
N1530 = N1530.Children;
N1540 = N1540.Children;
N1550 = N1550.Children;
N1570 = N1570.Children;
N1590 = N1590.Children;

%calculate average
N1280 = mean(mean(N1280.CData,'omitnan'),'omitnan');
N1300 = mean(mean(N1300.CData,'omitnan'),'omitnan');
N1320 = mean(mean(N1320.CData,'omitnan'),'omitnan');
N1520 = mean(mean(N1520.CData,'omitnan'),'omitnan');
N1530 = mean(mean(N1530.CData,'omitnan'),'omitnan');
N1540 = mean(mean(N1540.CData,'omitnan'),'omitnan');
N1550 = mean(mean(N1550.CData,'omitnan'),'omitnan');
N1570 = mean(mean(N1570.CData,'omitnan'),'omitnan');
N1590 = mean(mean(N1590.CData,'omitnan'),'omitnan');

XArray = [1280,1300,1320,1520,1530,1540,1550,1570,1590];

YArray = [N1280,N1300,N1320,N1520,N1530,N1540,N1550,N1570,N1590];

figure(10)
hold on
plot(XArray,YArray, 'LineWidth',2)
plot(XArray,YArray, 'x')
title('Average Amplitude')
xlabel('wavenumber (cm)')
ylabel('Amplitude average (A.U)')
hold off
