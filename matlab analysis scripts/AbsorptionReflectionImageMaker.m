eq%% image maker, imports in two files, one reflection one absorption and creates an amplitude and phase image.

close all

%folder of reflection
RefFolder = 'C:\Users\k68558dj\OneDrive - The University of Manchester\DataFiles\NPL Measurements\Jan 2023\Doped NW\High Dose\Candidate 1\1270 R';

%folder of Absorption
AbsFolder = 'C:\Users\k68558dj\OneDrive - The University of Manchester\DataFiles\NPL Measurements\Jan 2023\Doped NW\High Dose\Candidate 1\1270 I';

% Import the FIXED data
ReflectionHeight = strcat(RefFolder,'\', '1_1270 R_Height_PrimaryTrace_1.txt');

% Import the  MOVING data
AbsorptionHeight = strcat(AbsFolder,'\', '1_1270 I_Height_PrimaryTrace_1.txt');

%Do you wanna save??
Save = 'N';

SaveFolder = 'C:\Users\k68558dj\OneDrive - The University of Manchester\DataFiles\NPL Measurements\May 2024\Analysis\Doped NW\Cropped\1295';


alignTest(RefFolder,AbsFolder,ReflectionHeight,AbsorptionHeight)

%Size of image (um)
Xlength = 4;
Ylength = 4;




%% Set up the Import Options and import the data
opts = delimitedTextImportOptions("NumVariables", 256);

% Specify range and delimiter
opts.DataLines = [1, Inf];
opts.Delimiter = "\t";

% Specify column names and types
opts.VariableNames = ["e08", "e1", "e2", "e3", "e4", "e5", "e6", "e7", "e8", "e9", "e10", "e11", "e12", "e13", "e14", "e15", "e16", "e17", "e18", "e19", "e20", "e21", "e22", "e23", "e24", "e25", "e26", "e27", "e28", "e29", "e30", "e31", "e32", "e33", "e34", "e35", "e36", "e37", "e38", "e39", "e40", "e41", "e42", "e43", "e44", "e45", "e46", "e47", "e48", "e49", "e50", "e51", "e52", "e53", "e54", "e55", "e56", "e57", "e58", "e59", "e60", "e61", "e62", "e63", "e64", "e65", "e66", "e67", "e68", "e69", "e70", "e71", "e72", "e73", "e74", "e75", "e76", "e77", "e78", "e79", "e80", "e81", "e82", "e83", "e84", "e85", "e86", "e87", "e88", "e89", "e90", "e91", "e92", "e93", "e94", "e95", "e96", "e97", "e98", "e99", "e100", "e101", "e102", "e103", "e104", "e105", "e106", "e107", "e108", "e109", "e110", "e111", "e112", "e113", "e114", "e115", "e116", "e117", "e118", "e119", "e120", "e121", "e122", "e123", "e124", "e125", "e126", "e127", "e128", "e129", "e130", "e131", "e132", "e133", "e134", "e135", "e136", "e137", "e138", "e139", "e140", "e141", "e142", "e143", "e144", "e145", "e146", "e147", "e148", "e149", "e150", "e151", "e152", "e153", "e154", "e155", "e156", "e157", "e158", "e159", "e160", "e161", "e162", "e163", "e164", "e165", "e166", "e167", "e168", "e169", "e170", "e171", "e172", "e173", "e174", "e175", "e176", "e177", "e178", "e179", "e180", "e181", "e182", "e183", "e184", "e185", "e186", "e187", "e188", "e189", "e190", "e191", "e192", "e193", "e194", "e195", "e196", "e197", "e198", "e199", "e200", "e201", "e202", "e203", "e204", "e205", "e206", "e207", "e208", "e209", "e210", "e211", "e212", "e213", "e214", "e215", "e216", "e217", "e218", "e219", "e07", "e220", "e221", "e222", "e223", "e224", "e225", "e226", "e227", "e228", "e229", "e230", "e231", "e232", "e233", "e234", "e235", "e236", "e237", "e238", "e239", "e240", "e241", "e242", "e243", "e244", "e245", "e246", "e247", "e248", "e249", "e250", "e251", "e252", "e253", "e254"];
opts.VariableTypes = ["double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double"];

% Specify file level properties
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";



% Import the reflection data nf2 nf3 nf4
RefFileListShift = dir(strcat(RefFolder,'\shifted reflection data'));
for n = 1:size(RefFileListShift,1)
    %if it a txt file that is NOT a metaData
     if endsWith(RefFileListShift(n).name, '.txt') == 1 

         %if it is AFM
         RefAFM = readtable(ReflectionHeight, opts);

         %if it is NF2 
         if contains(RefFileListShift(n).name,'NF2','IgnoreCase',true) ==1
         RefNF2 = readtable(strcat(RefFolder,'\shifted reflection data','\',RefFileListShift(n).name), opts);
         end
         
         %if it is NF3
         if contains(RefFileListShift(n).name,'NF3','IgnoreCase',true) ==1
         RefNF3 = readtable(strcat(RefFolder,'\shifted reflection data','\',RefFileListShift(n).name), opts);
         end

         %if it is NF4
         if contains(RefFileListShift(n).name,'NF4','IgnoreCase',true) ==1
         RefNF4 = readtable(strcat(RefFolder,'\shifted reflection data','\',RefFileListShift(n).name), opts);
         end

         %if it is MCT Intensity
         if contains(RefFileListShift(n).name,'MCT','IgnoreCase',true) ==1
         RefMCT = readtable(strcat(RefFolder,'\shifted reflection data','\',RefFileListShift(n).name), opts);
         end
     end
end

% Import the absorption data nf2 nf3 nf4
AbsFileListShift = dir(strcat(AbsFolder,'\shifted absorption data'));
for n = 1:size(AbsFileListShift,1)
    %if it a txt file that is NOT a metaData
     if endsWith(AbsFileListShift(n).name, '.txt') == 1  && endsWith(AbsFileListShift(n).name, 'MetaData.txt') == 0 

         %if it is AFM
         AbsAFM = readtable(AbsorptionHeight,opts);

         %if it is NF2 
         if contains(AbsFileListShift(n).name,'NF2','IgnoreCase',true) ==1
         AbsNF2 = readtable(strcat(AbsFolder,'\shifted absorption data','\',AbsFileListShift(n).name), opts);
         end
         
         %if it is NF3
         if contains(AbsFileListShift(n).name,'NF3','IgnoreCase',true) ==1
         AbsNF3 = readtable(strcat(AbsFolder,'\shifted absorption data','\',AbsFileListShift(n).name), opts);
         end

         %if it is NF4
         if contains(AbsFileListShift(n).name,'NF4','IgnoreCase',true) ==1
         AbsNF4 = readtable(strcat(AbsFolder,'\shifted absorption data','\',AbsFileListShift(n).name), opts);
         end

         %if it is MCT Intensity
         if contains(AbsFileListShift(n).name,'MCT','IgnoreCase',true) ==1
         AbsMCT = readtable(strcat(AbsFolder,'\shifted absorption data','\',AbsFileListShift(n).name), opts);
         end
     end
end

%% Clear temporary variables
clear opts

%% handling
RefAFM = table2array(RefAFM);
AbsAFM = table2array(AbsAFM);

RefData2 = table2array(RefNF2);
AbsData2 = table2array(AbsNF2);

RefData3 = table2array(RefNF3);
AbsData3 = table2array(AbsNF3);

RefData4 = table2array(RefNF4);
AbsData4 = table2array(AbsNF4);

RefMCT = table2array(RefMCT);
AbsMCT = table2array(AbsMCT);

%check images are the same size
if size(RefData2)~= size(AbsData2)
    error('images are not the same size!')
end





%for each pixel, perform calculation, can you do this in one? I didnt want
%to fuck about with matrix calculations, its been a long time since Alevel
%maths
AmpData2 = zeros(size(RefData2,1),size(RefData2,2));
PhaseData2 = zeros(size(RefData2,1),size(RefData2,2));

AmpData3 = zeros(size(RefData2,1),size(RefData2,2));
PhaseData3 = zeros(size(RefData2,1),size(RefData2,2));

AmpData4 = zeros(size(RefData2,1),size(RefData2,2));
PhaseData4 = zeros(size(RefData2,1),size(RefData2,2));

MCTData = zeros(size(RefData2,1),size(RefData2,2));

for n = 1:size(RefData2,1)
    for i = 1:size(RefData2,2)
        AmpData2(n,i) = sqrt(RefData2(n,i)^2 + AbsData2(n,i)^2);
        PhaseData2(n,i) = atan(RefData2(n,i)/AbsData2(n,i));

        AmpData3(n,i) = sqrt(RefData3(n,i)^2 + AbsData3(n,i)^2);
        PhaseData3(n,i) = atan(RefData3(n,i)/AbsData3(n,i));

        AmpData4(n,i) = sqrt(RefData4(n,i)^2 + AbsData4(n,i)^2);
        PhaseData4(n,i) = atan(RefData4(n,i)/AbsData4(n,i));

        MCTData(n,i) = sqrt(RefMCT(n,i)^2 + AbsMCT(n,i)^2);
        MCTPhase(n,i) = atan(RefMCT(n,i)/AbsMCT(n,i));

    end
end

% %normalise amplitude data to MCT 
% for n = 1:size(RefData2,1)
%     for i = 1:size(RefData2,2)
% 
%         AmpData2(n,i) = AmpData2(n,i) / MCTData(n,i);
%         PhaseData2(n,i) = PhaseData2(n,i) - MCTPhase(n,i);
% 
%         AmpData3(n,i) = AmpData3(n,i) / MCTData(n,i);
%         PhaseData3(n,i) = PhaseData3(n,i) - MCTPhase(n,i);
% 
%         AmpData4(n,i) = AmpData4(n,i) / MCTData(n,i);
%         PhaseData4(n,i) = PhaseData4(n,i) - MCTPhase(n,i);
% 
% 
%     end
% end

%normalise to substrate
X0 = 1;
XEnd = 70;
Y0 = 1;
YEnd = 100;

AmpAvg2 = mean(AmpData2(X0:XEnd,Y0:YEnd),'all','omitnan');
AmpAvg3 = mean(AmpData3(X0:XEnd,Y0:YEnd),'all','omitnan');
AmpAvg4 = mean(AmpData4(X0:XEnd,Y0:YEnd),'all','omitnan');

PhaseAvg2 = mean(PhaseData2(X0:XEnd,Y0:YEnd),'all','omitnan');
PhaseAvg3 = mean(PhaseData3(X0:XEnd,Y0:YEnd),'all','omitnan');
PhaseAvg4 = mean(PhaseData4(X0:XEnd,Y0:YEnd),'all','omitnan');

% for n = 1:size(RefData2,1)
%     for i = 1:size(RefData2,2)
% 
%         AmpData2(n,i) = AmpData2(n,i) / AmpAvg2;
%         AmpData3(n,i) = AmpData3(n,i) / AmpAvg3;
%         AmpData4(n,i) = AmpData4(n,i) / AmpAvg4;
% 
%         PhaseData2(n,i) = PhaseData2(n,i) - PhaseAvg2;
%         PhaseData3(n,i) = PhaseData3(n,i) - PhaseAvg3;
%         PhaseData4(n,i) = PhaseData4(n,i) - PhaseAvg4;
% 
%     end
% end





%% plots

%calculate scale bar
%size of x axis in nm
XSize = 1000;
ScaleSize = 300;
LabelLength = (size(RefData2,1)/XSize)*ScaleSize;

LabelPosX = 300;
LabelPosY = 100;



figure(2)
subplot(1,2,1)
hold on
AmpImage2 = image(AmpData2,'CDataMapping','scaled');
colormap(gca,'hot')
colorbar
%uncomment if you want colourbars to be automatic
clim(gca,[0 prctile(AmpData2,99,"all")])
%clim(gca,[0,18]) %Here you can set the limits of colourbar
title('Amplitude NF2')
AmpImage2.XData = [0 Xlength];
AmpImage2.YData = [0 Ylength];
ax = gca;
xlabel('um')
ylabel('um')
ax.YDir = 'reverse';
xlim(ax,AmpImage2.XData)
ylim(ax,AmpImage2.YData)
pbaspect([Xlength Ylength 1])



subplot(1,2,2)
hold on
PhaseImage2 = image(PhaseData2,'CDataMapping','scaled');
colormap(gca,'cool')
colorbar
title('Phase NF2')
PhaseImage2.XData = [0 Xlength];
PhaseImage2.YData = [0 Ylength];
ax = gca;
xlabel('um')
ylabel('um')
ax.YDir = 'reverse';
xlim(ax,PhaseImage2.XData)
ylim(ax,PhaseImage2.YData)
pbaspect([Xlength Ylength 1])


figure(3)
subplot(1,2,1)
hold on
AmpImage3 = image(AmpData3,'CDataMapping','scaled');
colormap(gca,'hot')
colorbar
%uncomment if you want colourbars to be automatic
clim(gca,[0 prctile(AmpData3,99,"all")]) 
% clim(gca,[0,15]) %Here you can set the limits of colourbar
title('Amplitude NF3')
AmpImage3.XData = [0 Xlength];
AmpImage3.YData = [0 Ylength];
ax = gca;
xlabel('um')
ylabel('um')
ax.YDir = 'reverse';
xlim(ax,AmpImage3.XData)
ylim(ax,AmpImage3.YData)
pbaspect([Xlength Ylength 1])



subplot(1,2,2)
hold on
PhaseImage3 = image(PhaseData3,'CDataMapping','scaled');
colormap(gca,'cool')
colorbar
title('Phase NF3')
PhaseImage3.XData = [0 Xlength];
PhaseImage3.YData = [0 Ylength];
ax = gca;
xlabel('um')
ylabel('um')
ax.YDir = 'reverse';
xlim(ax,PhaseImage3.XData)
ylim(ax,PhaseImage3.YData)
pbaspect([Xlength Ylength 1])


figure(4)
subplot(1,2,1)
hold on
AmpImage4 = image(AmpData4,'CDataMapping','scaled');
colormap(gca,'hot')
colorbar
%uncomment if you want colourbars to be automatic
clim(gca,[0 prctile(AmpData4,99,"all")])
%clim(gca,[0,13]) %Here you can set the limits of colourbar
title('Amplitude NF4')
AmpImage4.XData = [0 Xlength];
AmpImage4.YData = [0 Ylength];
ax = gca;
xlabel('um')
ylabel('um')
% xline((XEnd/i)*Xlength,'linewidth',2,'Color','b')
% yline((YEnd/i)*Ylength,'linewidth',2,'Color','b')
ax.YDir = 'reverse';
xlim(ax,AmpImage4.XData)
ylim(ax,AmpImage4.YData)
pbaspect([Xlength Ylength 1])



subplot(1,2,2)
hold on
PhaseImage4 = image(PhaseData4,'CDataMapping','scaled');
colormap(gca,'cool')
colorbar
title('Phase NF4')
PhaseImage4.XData = [0 Xlength];
PhaseImage4.YData = [0 Ylength];
ax = gca;
xlabel('um')
ylabel('um')
ax.YDir = 'reverse';
xlim(ax,PhaseImage4.XData)
ylim(ax,PhaseImage4.YData)
pbaspect([Xlength Ylength 1])


figure(5)
tiledlayout(1,2)
nexttile
hold on
RefAFMImage = image(RefAFM,'CDataMapping','scaled');
colormap hot
title('Reflection AFM')
RefAFMImage.XData = [0 Xlength];
RefAFMImage.YData = [0 Ylength];
ax = gca;
xlabel('um')
ylabel('um')
ax.YDir = 'reverse';
xlim(ax,RefAFMImage.XData)
ylim(ax,RefAFMImage.YData)
pbaspect([Xlength Ylength 1])


nexttile
hold on
AbsAFMImage = image(AbsAFM,'CDataMapping','scaled');
title('Absorption AFM')
AbsAFMImage.XData = [0 Xlength];
AbsAFMImage.YData = [0 Ylength];
ax = gca;
xlabel('um')
ylabel('um')
ax.YDir = 'reverse';
xlim(ax,AbsAFMImage.XData)
ylim(ax,AbsAFMImage.YData)
pbaspect([Xlength Ylength 1])

figure(6)
tiledlayout(1,2)
nexttile
hold on
NormImage = image(MCTData,'CDataMapping','scaled');
colorbar
colormap hot
title('MCT intensity Amplitude')
NormImage.XData = [0 Xlength];
NormImage.YData = [0 Ylength];
ax = gca;
xlabel('um')
ylabel('um')
ax.YDir = 'reverse';
xlim(ax,NormImage.XData)
ylim(ax,NormImage.YData)
pbaspect([Xlength Ylength 1])

nexttile
hold on
PhaseImage = image(MCTPhase,'CDataMapping','scaled');
colorbar
colormap hot
title('MCT intensity Amplitude')
PhaseImage.XData = [0 Xlength];
PhaseImage.YData = [0 Ylength];
ax = gca;
xlabel('um')
ylabel('um')
ax.YDir = 'reverse';
xlim(ax,PhaseImage.XData)
ylim(ax,PhaseImage.YData)
pbaspect([Xlength Ylength 1])

if Save == 'Y'
saveas(figure(1),strcat(SaveFolder,'\shift.fig'));
saveas(figure(2),strcat(SaveFolder,'\NF2.fig'));
saveas(figure(3),strcat(SaveFolder,'\NF3.fig'));
saveas(figure(4),strcat(SaveFolder,'\NF4.fig'));
saveas(figure(5),strcat(SaveFolder,'\AFM.fig'));
% saveas(figure(6),strcat(SaveFolder,'\MCT.fig'));
end