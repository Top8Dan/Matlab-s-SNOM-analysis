%multiscan plotter
%idea is that it will plot the AFM and overlay each of the scan positions
%onto it. lets find out
clear all
close all

Harmonic = 3;
Source = 'C';

Save = 'Y'

SaveFile = 'C:\Users\k68558dj\OneDrive - The University of Manchester\DataFiles\TI NW SN123\Dan_BiTe_TargetedDopedNWs\2023-08-08 1869\Source C\Figures - Substrate norm to 1';

%put the gwy file in here
ImageFile = "C:\Users\k68558dj\OneDrive - The University of Manchester\DataFiles\TI NW SN123\Dan_BiTe_TargetedDopedNWs\2023-08-08 1869\Source C\2023-08-08 144253 WL DopedNW Targeted doping Source C WL\2023-08-08 144253 WL DopedNW Targeted doping Source C WL.gwy";

%folder of scans
ScanFile = "C:\Users\k68558dj\OneDrive - The University of Manchester\DataFiles\TI NW SN123\Dan_BiTe_TargetedDopedNWs\2023-08-08 1869\Source C\Scans";

RefFile = "C:\Users\k68558dj\OneDrive - The University of Manchester\DataFiles\TI NW SN123\Dan_BiTe_TargetedDopedNWs\2023-08-08 1869\Source C\Scans\2023-08-08 145731 NF S\2023-08-08 145731 NF S Interferograms.txt.txt";

%in neaspec software the z raw channel will be saved on channel 50
%typically, using a Gwy file converter found: https://uk.mathworks.com/matlabcentral/fileexchange/32893-gwyddion-file-importer
AFM = readgwychannel(ImageFile,50);
AFMImage = AFM.data;

directory = dir(ScanFile);
FileNumber = 1;
%isolate file names
for i = 1:size(directory,1)
    if directory(i).name(1) ~= '.'
        ScanDirectory(FileNumber) = directory(i);
        FileNumber = FileNumber+1;
    end
end

%source limits for plotting
if Source == 'A'
        WindowLeft = 600;
        WindowRight = 1250;

elseif Source == 'B'
        WindowLeft = 700;
        WindowRight = 1450;
    
elseif Source == 'C'
        WindowLeft = 900;
        WindowRight = 1550;

elseif Source == 'D'
        WindowLeft = 1150;
        WindowRight = 1850;
        
elseif Source == 'E'
        WindowLeft = 1200;
        WindowRight = 2100;
end

%extract AFM centre position
    %txt file
    AFMText = fileread(strcat(ImageFile{1}(1:strlength(ImageFile)-4),'.txt'));
    AFMPosition = extractBetween(AFMText,'Scanner Center Position (X, Y):	[µm]	','# Rotation:');
    AFMNameCell = extractBetween(AFMText,'Description:','# Date');
    AFMName = AFMNameCell{1};
    AFMPositionX = str2double(AFMPosition{1}(1:5));
    AFMPositionY = str2double(AFMPosition{1}(6:10));

    %Extract AFM scan size
    AFMSize = extractBetween(AFMText,'Scan Area (X, Y, Z):	[µm]	','	0.000');
    AFMSizeX = str2double(AFMSize{1}(1:5));
    AFMSizeY = str2double(AFMSize{1}(6:10));


%extract positions and data from scans
for i = 1:(size(ScanDirectory,2))
    %read text from file
    ScanText(i).name = ScanDirectory(i).name;
    ScanText(i).data = fileread(strcat(ScanFile,'\',ScanDirectory(i).name,'\',ScanDirectory(i).name,' Interferograms.txt.txt'));

    [ScanText(i).Spectra,ScanText(i).axis,ScanText(i).SpectraRaw] = NFSpectraFunction(strcat(ScanFile,'\',ScanDirectory(i).name,'\',ScanDirectory(i).name,' Interferograms.txt.txt'),RefFile,Harmonic,Source);
    

    %find scan position
    PositionCell(i) = extractBetween(ScanText(i).data,'Scanner Center Position (X, Y):	[µm]	','# Rotation:');
    ScanText(i).positionX = str2double(PositionCell{1,i}(1:5));
    ScanText(i).positionY = str2double(PositionCell{1,i}(6:10));

end


%convert um position to pixel location
PixelX = size(AFMImage,1);
PixelY = size(AFMImage,2);
AFMAxisX = linspace(AFMPositionX - (AFMSizeX/2),AFMPositionX + (AFMSizeX/2),PixelX);
AFMAxisY = linspace(AFMPositionY - (AFMSizeY/2),AFMPositionY + (AFMSizeY/2),PixelY);


%% plotting

figure(1)
hold on
Image = imagesc(AFMImage);
plot1 = plot([ScanText.positionX],[ScanText.positionY],'X');
plot1.MarkerSize = 8;
plot1.MarkerEdgeColor = [1,0,0];
Image.XData = [AFMAxisX(1) AFMAxisX(end)];
Image.YData = [AFMAxisY(1) AFMAxisY(end)];
ylim([AFMAxisY(1) AFMAxisY(end)])
xlim([AFMAxisX(1) AFMAxisX(end)])
set(gca,'ydir','reverse')
dx = 0.1; dy = 0.1;

    %sequential labels
    NameList = linspace(1,size(ScanDirectory,2),size(ScanDirectory,2));
    for i=1:size(ScanDirectory,2)
        text([ScanText(i).positionX]+dx,[ScanText(i).positionY]+dy,num2str(NameList(i)))
    end
    
    % %file labels
    % NameList = vertcat(ScanDirectory.name);
    % NameList = NameList(:,12:17);
    % text([ScanText.positionX]+dx,[ScanText.positionY]+dy,NameList)

xlabel('um')
ylabel('um')
title(AFMName)
hold off


% Create a waterfall plot
offset = 0.1;
%reference
figure(3)
hold on
for i = [1 3 5 7 9]
    axisAll = ScanText(i).axis(1,:);
    SpectraAbs = abs(ScanText(i).Spectra)+((i*offset));
    plot(axisAll,SpectraAbs)

end
xlim([WindowLeft WindowRight]);
title('Ref Spectra (abs)')
ylabel ('amplitude (arb, offset)')
xlabel('wavenumber cm^{-1}')
legend('1','3','5','7','9')
hold off

figure(4)
hold on
for i = [1 3 5 7 9]
    axisAll = ScanText(i).axis(1,:);
    Spectra = angle(ScanText(i).Spectra)+((i*offset));
    plot(axisAll,Spectra)
end
xlim([WindowLeft WindowRight]);
title('Ref Spectra (angle)')
ylabel ('angle (Phi, offset)')
xlabel('wavenumber cm^{-1}')
legend('1','3','5','7','9')
hold off

%sample
figure(5)
hold on
for i = [2 4 6 8 10]
    axisAll = ScanText(i).axis(1,:);
    SpectraAbs = abs(ScanText(i).Spectra)+((i*offset));
    plot(axisAll,SpectraAbs)

end
xlim([WindowLeft WindowRight]);
title('Sample Spectra (abs)')
ylabel ('amplitude (arb, offset)')
xlabel('wavenumber cm^{-1}')
legend('2','4','6','8','10')
hold off

figure(6)
hold on
for i = [2 4 6 8 10]
    axisAll = ScanText(i).axis(1,:);
    Spectra = angle(ScanText(i).Spectra)+((i*offset));
    plot(axisAll,Spectra)
end
xlim([WindowLeft WindowRight]);
title('Sample Spectra (angle)')
ylabel ('angle (Phi, offset)')
xlabel('wavenumber cm^{-1}')
legend('2','4','6','8','10')
hold off

%% no offset plotting

%reference
figure(7)
hold on
for i = [1 3 5 7 9]
    axisAll = ScanText(i).axis(1,:);
    SpectraAbs = abs(ScanText(i).Spectra);
    plot(axisAll,SpectraAbs)

end
xlim([WindowLeft WindowRight]);
title('Ref Spectra (abs)')
ylabel ('amplitude (arb)')
xlabel('wavenumber cm^{-1}')
legend('1','3','5','7','9')
hold off

figure(8)
hold on
for i = [1 3 5 7 9]
    axisAll = ScanText(i).axis(1,:);
    Spectra = angle(ScanText(i).Spectra);
    plot(axisAll,Spectra)
end
xlim([WindowLeft WindowRight]);
title('Ref Spectra (angle)')
ylabel ('angle (Phi)')
xlabel('wavenumber cm^{-1}')
legend('1','3','5','7','9')
hold off

%sample
figure(9)
hold on
for i = [2 4 6 8 10]
    axisAll = ScanText(i).axis(1,:);
    SpectraAbs = abs(ScanText(i).Spectra);
    plot(axisAll,SpectraAbs)

end
xlim([WindowLeft WindowRight]);
title('Sample Spectra (abs)')
ylabel ('amplitude (arb)')
xlabel('wavenumber cm^{-1}')
legend('2','4','6','8','10')
hold off

figure(10)
hold on
for i = [2 4 6 8 10]
    axisAll = ScanText(i).axis(1,:);
    Spectra = angle(ScanText(i).Spectra);
    plot(axisAll,Spectra)
end
xlim([WindowLeft WindowRight]);
title('Sample Spectra (angle)')
ylabel ('angle (Phi)')
xlabel('wavenumber cm^{-1}')
legend('2','4','6','8','10')
hold off

figure(11)
hold on
for i = [1 3 5 7 9]
    axisAll = ScanText(i).axis(1,:);
    Spectra = abs(ScanText(i).SpectraRaw);
    plot(axisAll,Spectra)
end
xlim([WindowLeft WindowRight]);
title('Ref Spectra (abs,raw)')
ylabel ('angle (Phi)')
xlabel('wavenumber cm^{-1}')
legend('2','4','6','8','10')
hold off

figure(12)
hold on
for i = [2 4 6 8 10]
    axisAll = ScanText(i).axis(1,:);
    Spectra = abs(ScanText(i).SpectraRaw);
    plot(axisAll,Spectra)
end
xlim([WindowLeft WindowRight]);
title('Sample Spectra (abs,raw)')
ylabel ('angle (Phi)')
xlabel('wavenumber cm^{-1}')
legend('2','4','6','8','10')
hold off

%% Saving

 if Save == 'Y'
     saveas(figure(1),strcat(SaveFile,'\','Source',Source,'AFM.fig'));
     saveas(figure(3),strcat(SaveFile,'\','Source',Source,'RefSpectraAbs Offset.fig'));
     saveas(figure(4),strcat(SaveFile,'\','Source',Source,'RefSpectraAngle Offset.fig'));
     saveas(figure(5),strcat(SaveFile,'\','Source',Source,'SpectraAbs Offset.fig'));
     saveas(figure(6),strcat(SaveFile,'\','Source',Source,'SpectraAngle Offset.fig'));
     saveas(figure(7),strcat(SaveFile,'\','Source',Source,'RefSpectraAbs.fig'));
     saveas(figure(8),strcat(SaveFile,'\','Source',Source,'RefSpectraAngle.fig'));
     saveas(figure(9),strcat(SaveFile,'\','Source',Source,'SpectraAbs.fig'));
     saveas(figure(10),strcat(SaveFile,'\','Source',Source,'SpectraAngle.fig'));
     saveas(figure(11),strcat(SaveFile,'\','Source',Source,'RawRefSpectra.fig'));
     saveas(figure(12),strcat(SaveFile,'\','Source',Source,'RawSampleSpectra.fig'));
 end