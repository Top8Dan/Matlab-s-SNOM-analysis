clear all
close all
%% import the spectra to be windowed
opts = delimitedTextImportOptions("NumVariables", 18);

% Specify range and delimiter
opts.DataLines = [31, Inf];
opts.Delimiter = "\t";

% Specify column names and types
opts.VariableNames = ["R", "C", "Run", "Depth", "M", "O0A", "O0P", "O1A", "O1P", "O2A", "O2P", "O3A", "O3P", "O4A", "O4P", "O5A", "O5P", "VarName18"];
opts.VariableTypes = ["double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "string"];

% Specify file level properties
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";

% Specify variable properties
opts = setvaropts(opts, "VarName18", "WhitespaceRule", "preserve");
opts = setvaropts(opts, "VarName18", "EmptyFieldRule", "auto");
opts = setvaropts(opts, ["R", "C", "Run", "Depth", "M", "O0A", "O0P", "O1A", "O1P", "O2A", "O2P", "O3A", "O3P", "O4A", "O4P", "O5A", "O5P"], "ThousandsSeparator", ",");



%% Variables

% Import the data
NFSpectraInterferogramsRaw = readtable("C:\Users\k68558dj\OneDrive - The University of Manchester\DataFiles\2022-09-30_Tom_MnImplantedBi2Se3\2022-10-04 1785\2022-10-04 172054 NF LS tip0_sourceB_film_lineSpectra\2022-10-04 172054 NF LS tip0_sourceB_film_lineSpectra Interferograms.txt.txt", opts);

% import Ref data
NFSpectraInterferogramsRawRef = readtable("C:\Users\k68558dj\OneDrive - The University of Manchester\DataFiles\2022-09-30_Tom_MnImplantedBi2Se3\2022-10-04 1780\2022-10-04 110133 NF S tip0_sourceB_Au_spectrum\2022-10-04 110133 NF S tip0_sourceB_Au_spectrum Interferograms.txt.txt", opts);

%Folder to Save in
SaveFolder = "C:\Users\k68558dj\OneDrive - The University of Manchester\NWPaper\Spectra (fig1)\Mn Doped Thin Film\All Avg\Source B\3rd Harmonic";

%Save?
Save = "Y";

%Harmonic?
Harmonic = 3;

% Source to be plotted (cm-1)
Source = 'B';

%% Windowing variables
%zero padding
ExtensionFactor = 1;

%Blackman Harris Shape
FlatTopSize = 0.3;
BHSize = 0.2;

%need to take avg of start and end flat lines for complex column (estimated length = 400
%datapoints)
AvgLength = 400;

%% Phase correction variables
%Integration time (s)
IntegrationTime = 100e-3; % 100 ms not 0.1s 

% time gap between scan start times (s)
TimeGap = 1200;

%PhaseCorrection? (Y/N)
PhaseCorrection = 'N';

%Reference measured first? (Y/N)
ReferenceFirst = 'Y';

%% FDM Variables
%parameters
Radius = 35e-9; %tip radius
Tmin = 1e-9; %tip minimum distance
Tamp = 50e-9; %tip amplitude position
freq = 1; %tip frequency
T = 1/freq; %tip period

% finite Dipole parameters
g = 0.7*exp(0.06*1i);
L = 300e-9;

%reference calculation: this is what the fdm model will normalise its
    %fitting to so must match what you are normalising your normal data to.
    %E_sRef = Dielectric_Au((1200));
    E_sRef = 11.7 +0.1i;
    %E_sRef =  4.71 +0.66i;

%% CODE BELOW




    if Source == 'A'
        WindowLeft = 850;
        WindowRight = 1250;

    elseif Source == 'B'
        WindowLeft = 900;
        WindowRight = 1450;
    
    elseif Source == 'C'
        WindowLeft = 1100;
        WindowRight = 1600;

    elseif Source == 'D'
        WindowLeft = 1350;
        WindowRight = 1900;
        
    elseif Source == 'E'
        WindowLeft = 1350;
        WindowRight = 2000;
    end


%% Run twice, one for data once for Ref
for x=1:2
    %second run for Ref
    if x == 2
        NFSpectraInterferogramsRaw = NFSpectraInterferogramsRawRef;
    end

%% seperate runs into a cell array

%remove metadata issues
for i = 1:10
    if NFSpectraInterferogramsRaw.Depth(i) == 0
    NFSpectraInterferogramsRaw = NFSpectraInterferogramsRaw(i:end,:);
    end
end

%For a line scan, convert extra coloumns to extra runs to average across.
RunsPerColoumn = max(NFSpectraInterferogramsRaw.Run(:));
for i = 1: size(NFSpectraInterferogramsRaw,1)
    NFSpectraInterferogramsRaw.Run(i) = NFSpectraInterferogramsRaw.Run(i) + (NFSpectraInterferogramsRaw.C(i)*RunsPerColoumn);
end


%amount of runs
Runs = max(NFSpectraInterferogramsRaw.Run);

%amount of measurements per run
Resolution = 0;
while NFSpectraInterferogramsRaw.Run(Resolution+1) == 0
    Resolution = Resolution + 1;
end



%create cell array
NFSpectraInterferograms = cell(1,Runs+1);

%seperate each measurement into a cell array
for n = 1:Runs+1
    NFSpectraInterferograms(n) = {NFSpectraInterferogramsRaw(((n-1)*Resolution+1):n*Resolution,:)};
end

%add complex value column to table
if Harmonic == 3
for n = 1:Runs+1
    NFSpectraInterferograms{n}.Complex = NFSpectraInterferograms{n}.O3A .* exp(1i.*NFSpectraInterferograms{n}.O3P);
end
end

if Harmonic == 4
for n = 1:Runs+1
    NFSpectraInterferograms{n}.Complex = NFSpectraInterferograms{n}.O4A .* exp(1i.*NFSpectraInterferograms{n}.O4P);
end
end
%% handling


%% zero padding

%how many lengths of zeros to add either side.


MLength = NFSpectraInterferograms{1}.M(end) - NFSpectraInterferograms{1}.M(1); % overall length of intefergoram (should be interferogram distance - 100um checked and it is) 


ZeroArray = zeros(Resolution*ExtensionFactor, size(NFSpectraInterferograms{n},2));
ZeroTable = array2table(ZeroArray,'VariableNames',["R", "C", "Run", "Depth", "M", "O0A", "O0P", "O1A", "O1P", "O2A", "O2P", "O3A", "O3P", "O4A", "O4P", "O5A", "O5P", "VarName18", "Complex"]);

% define distance step
dX = median(diff(NFSpectraInterferograms{1}.M));



%zero padding
for n = 1:Runs+1
    %padding for complex
    StartAvgComplex(n) = mean(NFSpectraInterferograms{n}.Complex(1:AvgLength));
    EndAvgComplex(n) = mean(NFSpectraInterferograms{n}.Complex(end-AvgLength:end));

   

    %padding for position
    StartPositionArray(n,:) = linspace(NFSpectraInterferograms{n}.M(1)-(dX*ExtensionFactor*Resolution),NFSpectraInterferograms{n}.M(1)-dX,ExtensionFactor*Resolution);
    EndPositionArray(n,:) = linspace(NFSpectraInterferograms{n}.M(end)+dX,NFSpectraInterferograms{n}.M(end)+(dX*ExtensionFactor*Resolution),ExtensionFactor*Resolution);

    StartTable = ZeroTable;
    EndTable = ZeroTable;

    for i = 1:size(ZeroTable,1)
        StartTable.Complex(i) = StartAvgComplex(n);%complex
        EndTable.Complex(i) = EndAvgComplex(n);
        StartTable.M(i) = StartPositionArray(n,i); %position
        EndTable.M(i) = EndPositionArray(n,i);
        
    end
    NFSpectraInterferograms{n} = [StartTable; NFSpectraInterferograms{n}; EndTable];
end


%% set mirror position = 0 at max intensity
for n = 1:Runs+1
    %find max position
   [Max,iMax] = max(abs(NFSpectraInterferograms{n}.Complex));
   %assign M = 0 at max position
   NFSpectraInterferograms{n}.M = NFSpectraInterferograms{n}.M - NFSpectraInterferograms{n}.M(iMax);
end

%% Windowing

%% blackman Harris windowing

% %what percentage (from the middle of the signal) do you want to window?
% WindowPercent = 0.9;
% 
% BHWindow = blackmanharris(round(runsize*WindowPercent));
% 
% %add surrounding zeros
% BHZeros = zeros(round((runsize*(1-WindowPercent))/2),1);
% BHWindow = [BHZeros; BHWindow; BHZeros];
% 


%% flat top BH Windowing

% FlatTopSize = 0.80;
% 
% FT = ones(round((Resolution*FlatTopSize)),1);
% BH = blackmanharris(Resolution - size(FT,1)-1);
% 
% Window  = [BH(1:round(size(BH)/2)); FT ; BH(round(size(BH)/2):end)];
% 
% %add zero padding
% Window = [zeros((Resolution*ExtensionFactor),1); Window; zeros((Resolution*ExtensionFactor),1)];


%% Offset BH Windowing
% 
% 
% FT = ones(round((Resolution*FlatTopSize)),1);
% BH = blackmanharris(Resolution*BHSize);
% 
% Window  = [BH(1:round(size(BH)/2)); FT ; BH(round(size(BH)/2):end)];
% 
% WindowZeros = zeros(Resolution - length(Window),1);
% 
% 
% WindowPosition = iMax - (Resolution*ExtensionFactor) - (round(length(Window)/2));
% 
% Window = [WindowZeros(1:WindowPosition,1); Window; WindowZeros(WindowPosition+1:end,1)];
% 


%cut to size if BH goes over Right edge
FT = ones(round((Resolution*FlatTopSize)),1);
BH = blackmanharris(Resolution*BHSize);

Window  = [BH(1:round(size(BH)/2)); FT ; BH(round(size(BH)/2):end)];

WindowZeros = zeros(length(NFSpectraInterferograms{1}.Complex) - length(Window),1);


%% THIS IS REPORTING TOO BIG
WindowPosition = iMax -(round(length(Window)/2));



Window = [WindowZeros(1:WindowPosition,1); Window; WindowZeros(WindowPosition+1:end,1)];

 %if window clips right hand side, smooth it out.
if ((size(FT)/2) + WindowPosition -(Resolution*ExtensionFactor) ) > Resolution
    Window((Resolution*(ExtensionFactor+1))-round(size(BH)/2):(Resolution*(ExtensionFactor+1))) = BH(round(size(BH)/2):end);
end

%add zero padding
%Window = [zeros((Resolution*ExtensionFactor),1); Window; zeros((Resolution*ExtensionFactor),1)];



for n = 1:Runs+1
    %avg of each 
    ZeroShift(n) = (StartAvgComplex(n)+EndAvgComplex(n))/2;

    % apply window
    NFSpectraInterferograms{n}.Complex = NFSpectraInterferograms{n}.Complex-ZeroShift(n);%bring complex data centralised on zero for windowing 
    NFSpectraInterferograms{n}.Complex = NFSpectraInterferograms{n}.Complex.*Window;%apply window
%     NFSpectraInterferograms{n}.Complex = NFSpectraInterferograms{n}.Complex+ZeroShift(n);%bring complex data back to value
    
end




%% fft function
for n = 1:Runs+1
    %perform fft
    NFSpectraInterferograms{n}.FFT = fft(NFSpectraInterferograms{n}.Complex);
    NFSpectraInterferograms{n}.FFT = fftshift(NFSpectraInterferograms{n}.FFT);

    %creating frequency axis for FFT (note: this doesn't need to be a loop
    %as it should be the same for all interferograms - so you can do this
    %outside loop save time - Not the Same sadly :(
    t(n,:) = 2.*(NFSpectraInterferograms{n}.M)./3e8;  % in seconds
    Npts = length(t(n,:));
    doubleBW(n) = Npts/(t(n,end)-t(1,1));
    Freq(n,:) = linspace(-doubleBW(n)/2, doubleBW(n)/2, Npts); %in Hz

    wavenumberAxis(n,:) = Freq(n,:)/3e8;

    wavenumberAxisCm(n,:) = wavenumberAxis(n,:)/100;

    %find window values
    for i = 1:size(wavenumberAxisCm(n,:),2)
        if wavenumberAxisCm(n,i) < WindowLeft
            FFTWindow(1) = i;
        end
        if wavenumberAxisCm(n,i) < WindowRight
            FFTWindow(2) = i;
        end
    end

%extracting amplitude and phase
%     Sn(n,:) = abs(NFSpectraInterferograms{n}.FFT);
%     phi(n,:) = angle(NFSpectraInterferograms{n}.FFT);

%     %output data file or ref file
%     if x == 1
%          FFTComplex(n,:) = Sn(n,:).*exp(1i*phi(n,:)); %back to complex  
%     else
%         FFTComplexRef(n,:) = Sn(n,:).*exp(1i*phi(n,:)); %back to complex
%     end
% 
% end

    %output data file or ref file
    if x == 1
         FFTComplex(n,:) = NFSpectraInterferograms{n}.FFT; %back to complex  
    else
        FFTComplexRef(n,:) = NFSpectraInterferograms{n}.FFT; %back to complex
    end

end


%end of Mega loop
end


%average data
for n =1:size(FFTComplex,2)
    AvgSpectra(n) = mean(FFTComplex(:,n));
    AvgRef(n) = mean(FFTComplexRef(:,n));
end

%% Phase Correction

for n=1:Runs+1 %create phase drift plots
%Fit method
    %normalise scan to avg scan
    FFTComplexNorm(n,:) = FFTComplex(n,:) ./ AvgSpectra(1,:);
    FFTComplexNormRef(n,:) = FFTComplexRef(n,:) ./ AvgRef(1,:);

    %fit straight line
    FFTWLPFit{n} = fit(wavenumberAxis(n,FFTWindow(1):FFTWindow(2)).',(angle(FFTComplexNorm(n,FFTWindow(1):FFTWindow(2)).')),'poly1');
    FFTWLPFitRef{n} = fit(wavenumberAxis(n,FFTWindow(1):FFTWindow(2)).',(angle(FFTComplexNormRef(n,FFTWindow(1):FFTWindow(2)).')),'poly1');
    %convert gradient to WLP
    Gradient(n,:) = coeffvalues(FFTWLPFit{n});
    GradientRef(n,:) = coeffvalues(FFTWLPFitRef{n});
    dWLP(n) = Gradient(n,1)/(2*pi);
    dWLPRef(n) = GradientRef(n,1)/(2*pi); 

end

%% fit WLP drift
%create t axis

RefX = 1:Runs+1;
SampleX = 1:Runs+1;

%% fit of WLP against time
if ReferenceFirst == 'Y'

    RefXt = RefX*IntegrationTime*Resolution; % calculate length of scan from integration time
    TimeGap = TimeGap - RefXt(end);
    SampleXt = SampleX*IntegrationTime*Resolution+TimeGap+RefXt(end); % calculate length of scan from integration time
    tAxis = [RefXt SampleXt]; %universal t axis taking into account time between scans
    
    RefWLPFitFit = fit(RefXt(:),dWLPRef(:),'poly1'); %perform fit of drift
    SampleWLPFitFit = fit(SampleXt(:),dWLP(:),'poly1'); 
    
    RefWLPShift = 0;
    SampleWLPShift = RefWLPFitFit(SampleXt(1)) - SampleWLPFitFit(SampleXt(1)); 

elseif ReferenceFirst == 'N'

    SampleXt = SampleX*Resolution*IntegrationTime; % calculate length of scan from integration time
    TimeGap = TimeGap - SampleXt(end);
    RefXt = RefX*IntegrationTime*Resolution+TimeGap+SampleXt(end); % calculate length of scan from integration time
    tAxis = [SampleXt RefXt];  %universal t axis taking into account time between scans

    RefWLPFitFit = fit(RefXt(:),dWLPRef(:),'poly1'); %perform fit of drift
    SampleWLPFitFit = fit(SampleXt(:),dWLP(:),'poly1');

    RefWLPShift = SampleWLPFitFit(RefXt(1)) - RefWLPFitFit(RefXt(1));
    SampleWLPShift = 0;    
end

RefWLPCorrectionFit = RefWLPFitFit(RefXt) + RefWLPShift;
SampleWLPCorrectionFit = SampleWLPFitFit(SampleXt) + SampleWLPShift;

 

%apply phase correction
if PhaseCorrection == 'Y'
    for n = 1:Runs+1
        %convert from WLP to phase

        SamplePhaseCorrection(n,:) = SampleWLPCorrectionFit(n).*2.*pi.*wavenumberAxis(n,:);
        RefPhaseCorrection(n,:) = RefWLPCorrectionFit(n).*2.*pi.*wavenumberAxis(n,:);

        FFTComplex(n,:) = FFTComplex(n,:).*exp(-1i*SamplePhaseCorrection(n,:));
        FFTComplexRef(n,:) = FFTComplexRef(n,:).*exp(-1i*RefPhaseCorrection(n,:));

    end
end

%Normalise data by Ref
for n=1:Runs+1
NormalisedSpectra(n,:) = FFTComplex(n,:) ./ FFTComplexRef(n,:);
end

% avg of all runs
for n =1:size(NormalisedSpectra,2)
    
    NormalisedSpectraAvg(n) = mean(NormalisedSpectra(:,n));
    FFTComplexAvg(n) = mean(FFTComplex(:,n));

end

%% fitting of FDM



% discard anything that is noise outside of the range of the source being used.
range = wavenumberAxisCm(2,:) >= WindowLeft & wavenumberAxisCm(2,:) <= WindowRight;

NormalisedSpectraAvgWindow = NormalisedSpectraAvg(range);
wavenumberAxisCmWindow = wavenumberAxisCm(2,range);

%% fitting

resolutionwindow = size(wavenumberAxisCmWindow,2);

for z = 1:resolutionwindow
    ParamsOut = FDMBetafit(wavenumberAxisCmWindow(z), NormalisedSpectraAvgWindow(z), Tmin,Tamp,freq,Harmonic,g,Radius,L,E_sRef);
    BetaFit(z) = ParamsOut(1);
    E_sRefA(z) = E_sRef;
end

%Calculate epsilon
%E_s = (BetaFit+1)./(1-BetaFit);

%new eqn
E_s = (1-(BetaFit+1))./(BetaFit-1);


%% noise plotting
for n = 1:Runs+1 
    LogRefRaw(n,:) = log(abs(FFTComplexRef(n,:)));
end
LogRefRawAvg = mean(LogRefRaw(:,:));


%% plotting
NoPaddingL = (Resolution*ExtensionFactor+1);
NoPaddingR = Resolution*(ExtensionFactor+1);
figure(1)
    %plot(Window)
    hold on
    for n = 1:Runs+1
        plot(NFSpectraInterferograms{n}.M(NoPaddingL:NoPaddingR)./1e-6,real(NFSpectraInterferograms{n}.Complex(NoPaddingL:NoPaddingR)./max(NFSpectraInterferograms{n}.Complex(NoPaddingL:NoPaddingR))))
    end

    plot(NFSpectraInterferograms{n}.M(NoPaddingL:NoPaddingR)./1e-6,Window(NoPaddingL:NoPaddingR))

    title(strcat('windowing',' FT',num2str(FlatTopSize),' BH',num2str(BHSize)))
    xlabel('um')
    ylabel('real amp')
    legend
    %xlim([Resolution*ExtensionFactor,Resolution*(ExtensionFactor+1)])
    hold off

%Amplitude Plot
figure(2)
    for n =1:Runs+1
        hold on
        plot(wavenumberAxisCm(n,:),abs(FFTComplex(n,:)))
        
    end
    title('Absolute (Sample)')
    xlabel('wavenumber (cm-1)')
    ylabel('arb')
    xlim([500,2600])
    hold off

%phase plot
figure(3)
    for n =1:Runs+1
        hold on
        plot(wavenumberAxisCm(n,:),angle(FFTComplex(n,:)))
        
    end
    title('phase (Sample)')
    xlabel('wavenumber (cm-1)')
    ylabel('\phi')
    xlim([800,2600])
    hold off

figure(4)
    hold on
    plot(wavenumberAxisCm(2,:),abs(NormalisedSpectraAvg))
    title('Absolute (Norm,avg)')
    if PhaseCorrection == 'Y'
        title('Absolute (Norm,avg,Phase Corrected)')
    end
    xlabel('wavenumber (cm-1)')
    ylabel('arb')
    xlim([WindowLeft,WindowRight])

    hold off


figure(5)

    plot(wavenumberAxisCm(2,:),(angle(NormalisedSpectraAvg)))
    title('Angle (Norm,avg)')
    if PhaseCorrection == 'Y'
        title('Angle (Norm,avg,Phase Corrected)')
    end
    xlabel('wavenumber (cm-1)')
    ylabel('\phi')
    xlim([WindowLeft,WindowRight])
    ylim([-3.14,3.14])
    hold off

figure(6)
    for n = 1:Runs+1
        hold on
        plot(wavenumberAxisCm(n,:),abs(NormalisedSpectra(n,:)))
        
    end
    title('Absolute (Norm, All runs)')
    xlim([WindowLeft,WindowRight])
    title('Normalisied abs (induvidual)')
    if PhaseCorrection == 'Y'
        title('Normalisied abs (induvidual) (phase corrected)')
    end
    xlabel('wavenumber (cm-1)')
    ylabel('arb')
    legend('1','2','3','4','5','6','7','8','9','10')
    hold off


figure(7)
    for n = 1:Runs+1
        hold on
        plot(wavenumberAxisCm(n,FFTWindow(1): FFTWindow(2)),unwrap(angle(NormalisedSpectra(n,FFTWindow(1): FFTWindow(2)))))
       
    end
    title('Angle Sample (Norm, All runs)')
    xlim([WindowLeft,WindowRight])
    title('Normalisied angle (induvidual)')
    if PhaseCorrection == 'Y'
        title('Normalisied angle (induvidual) (phase corrected)')
    end
    xlabel('wavenumber (cm-1)')
    ylabel('\phi')
    legend('1','2','3','4','5','6','7','8','9','10')
    hold off


figure(10)
for n = 1:Runs+1
        hold on
        plot(wavenumberAxis(n,FFTWindow(1):FFTWindow(2)),(angle(FFTComplexNorm(n,FFTWindow(1):FFTWindow(2)).')))
        plot(wavenumberAxis(n,FFTWindow(1):FFTWindow(2)),FFTWLPFit{n}(wavenumberAxis(n,FFTWindow(1):FFTWindow(2))))

end
legend('1','1F','2','2F','3','3F','4','4F','5','5F','6','6F','7','7F','8','8F','9','9F','10','10F')
title('Normalised to Avg phase')
xlabel('wavenumber (m-1)')
ylabel('\phi')
hold off

figure(11)
 hold on
 plot(SampleXt,(dWLP+SampleWLPShift)/1e-9)
 plot(SampleXt,SampleWLPCorrectionFit/1e-9,'--')
 plot(RefXt,(dWLPRef+RefWLPShift)/1e-9)
 plot(RefXt,RefWLPCorrectionFit/1e-9,'--')
 title('WLP Drift (Fit Method)')
 xlabel('time (s)')
 ylabel('\DeltaWLP (nm)')
 legend('Ref Phase','Ref Phase Fit', 'Sample Phase','Sample Phase Fit')
 hold off

figure(12)
    plot(wavenumberAxisCm(2,:),abs(FFTComplexAvg))
    title('Absolute avg all runs (avg)')
    if PhaseCorrection == 'Y'
        title('Absolute (avg,Phase Corrected)')
    end
    xlabel('wavenumber (cm-1)')
    ylabel('arb')
    xlim([WindowLeft,WindowRight])
    hold off

figure(13)
    plot(wavenumberAxisCm(2,:),angle(FFTComplexAvg))
    title('Angle avg all runs (avg)')
    if PhaseCorrection == 'Y'
        title('Absolute (avg,Phase Corrected)')
    end
    xlabel('wavenumber (cm-1)')
    ylabel('arb')
    xlim([WindowLeft,WindowRight])
    hold off

%Amplitude Plot
figure(14)
    for n =1:Runs+1
        hold on
        plot(wavenumberAxisCm(n,:),abs(FFTComplexRef(n,:)))
        
    end
    title('Absolute (Ref)')
    xlabel('wavenumber (cm-1)')
    ylabel('arb')
    xlim([800,2600])
    hold off

%phase plot
figure(15)
    for n =1:Runs+1
        hold on
        plot(wavenumberAxisCm(n,:),angle(FFTComplexRef(n,:)))
        
    end
    title('phase (Ref)')
    xlabel('wavenumber (cm-1)')
    ylabel('\phi')
    xlim([800,2600])
    hold off

figure(16)
    hold on
    title('E_s fitted (real+imag)')
    %yyaxis left
    plot(wavenumberAxisCmWindow, real(E_s),'-b')
    plot(wavenumberAxisCmWindow, real(E_sRefA),'--b')
    legend('EsReal','realFit')
    ylabel('E_s')
    yyaxis right
    plot(wavenumberAxisCmWindow, imag(E_s),'-r')
    plot(wavenumberAxisCmWindow, imag(E_sRefA),'--r')
    xlabel('Wavenumber (cm-1)')
    legend('EsReal', 'EsImag')
    xlim([WindowLeft WindowRight])
    hold off

% figure(17)
%     hold on
%     title('E_s Ref (real+imag)')
%     %yyaxis left
%     plot(wavenumberAxisCmWindow, real(E_sRefA),'-b')
%     legend('EsReal','realFit')
%     ylabel('E_s')
%     yyaxis right
%     plot(wavenumberAxisCmWindow, imag(E_sRefA),'-r')
%     xlabel('Wavenumber (cm-1)')
%     legend('EsReal', 'EsImag')
%     xlim([WindowLeft WindowRight])
%     hold off

figure(18)
    hold on

    plot(wavenumberAxisCm(2,:),abs(NormalisedSpectraAvg))   
    plot(wavenumberAxisCm(2,:),LogRefRawAvg)
    title('Signal to noise (Norm,avg)')
    xlabel('wavenumber (cm-1)')
    ylabel('arb')
    xlim([WindowLeft-150,WindowRight+150])
    xregion(0,WindowLeft)
    xregion(WindowRight,4000)
    hold off




 if Save == 'Y'
     %save figures
     saveas(figure(1),strcat(SaveFolder,'\','Source',Source,'interferograms.fig'));
     saveas(figure(2),strcat(SaveFolder,'\','Source',Source,'Amp Raw.fig'));
     saveas(figure(3),strcat(SaveFolder,'\','Source',Source,'Angle Raw.fig'));
     saveas(figure(4),strcat(SaveFolder,'\','Source',Source,'Avg Amp.fig'));
     saveas(figure(5),strcat(SaveFolder,'\','Source',Source,'Avg Angle.fig'));
     saveas(figure(6),strcat(SaveFolder,'\','Source',Source,'Induvidual Amp.fig'));
     saveas(figure(7),strcat(SaveFolder,'\','Source',Source,'Induvidual Angle.fig'));
     saveas(figure(11),strcat(SaveFolder,'\','Source',Source,'Drift Correction.fig'));
     saveas(figure(12),strcat(SaveFolder,'\','Source',Source,'Avg Amp Raw.fig'));
     saveas(figure(13),strcat(SaveFolder,'\','Source',Source,'Avg Angle Raw.fig'));
     saveas(figure(14),strcat(SaveFolder,'\','Source',Source,'Amp Raw Ref.fig'));
     saveas(figure(15),strcat(SaveFolder,'\','Source',Source,'Angle Raw Ref.fig'));
     saveas(figure(16),strcat(SaveFolder,'\','Source',Source,'FDMDielectric.fig'))
     saveas(figure(18),strcat(SaveFolder,'\','Source',Source,'SignalToNoise.fig'))
     

    %save output data
    outputAmp(:,1) = wavenumberAxisCm(2,:);
    outputAmp(:,2) = abs(NormalisedSpectraAvg);

    outputAngle(:,1) = wavenumberAxisCm(2,:);
    outputAngle(:,2) = angle(NormalisedSpectraAvg);


    writematrix(outputAmp,strcat(SaveFolder,'\','Source',Source,'AvgAmp.txt'))
    writematrix(outputAngle,strcat(SaveFolder,'\','Source',Source,'AvgAngle.txt'))

 end

