function [Spectra,wavenumberAxisCm,SpectraRaw] = NFSpectraFunction(FileAddress,RefFileAddress,Harmonic,Source)

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

% Import the data
NFSpectraInterferogramsRaw = readtable(FileAddress, opts);

% import Ref data
NFSpectraInterferogramsRawRef = readtable(RefFileAddress, opts);


%set the source
if Source == 'A'
        WindowLeft = 800;
        WindowRight = 1250;

elseif Source == 'B'
        WindowLeft = 700;
        WindowRight = 1500;
    
elseif Source == 'C'
        WindowLeft = 1150;
        WindowRight = 1550;

elseif Source == 'D'
        WindowLeft = 1350;
        WindowRight = 1900;
        
elseif Source == 'E'
        WindowLeft = 1850;
        WindowRight = 2100;
end

%% Run twice, one for data once for Ref
for x=1:2
    %second run for Ref
    if x == 2
        NFSpectraInterferogramsRaw = NFSpectraInterferogramsRawRef;
    end

%% seperate runs into a cell array

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
ExtensionFactor = 2;

MLength = NFSpectraInterferograms{1}.M(end) - NFSpectraInterferograms{1}.M(1); % overall length of intefergoram (should be interferogram distance - 100um checked and it is) 


%need to take avg of start and end flat lines for complex column (estimated length = 400
%datapoints)
AvgLength = 400;

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

%% Offset BH Windowing
FlatTopSize = 0.3;
BHSize = 0.1;

FT = ones(round((Resolution*FlatTopSize)),1);
BH = blackmanharris(Resolution*BHSize);

Window  = [BH(1:round(size(BH)/2)); FT ; BH(round(size(BH)/2):end)];

WindowZeros = zeros(Resolution - length(Window),1);

WindowPosition = iMax - (Resolution*ExtensionFactor) - (round(length(Window)/2));

Window = [WindowZeros(1:WindowPosition,1); Window; WindowZeros(WindowPosition+1:end,1)];

%cut to size if BH goes over Right edge

 

%add zero padding
Window = [zeros((Resolution*ExtensionFactor),1); Window; zeros((Resolution*ExtensionFactor),1)];



for n = 1:Runs+1
    %avg of each 
    ZeroShift(n) = (StartAvgComplex(n)+EndAvgComplex(n))/2;

    % apply window
    NFSpectraInterferograms{n}.Complex = NFSpectraInterferograms{n}.Complex-ZeroShift(n);%bring complex data centralised on zero for windowing 
    NFSpectraInterferograms{n}.Complex = NFSpectraInterferograms{n}.Complex.*Window;%apply window
    NFSpectraInterferograms{n}.Complex = NFSpectraInterferograms{n}.Complex+ZeroShift(n);%bring complex data back to value
    
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
    Sn(n,:) = abs(NFSpectraInterferograms{n}.FFT);
    phi(n,:) = angle(NFSpectraInterferograms{n}.FFT);

    %output data file or ref file
    if x == 1
         FFTComplex(n,:) = Sn(n,:).*exp(1i*phi(n,:)); %back to complex  
    else
        FFTComplexRef(n,:) = Sn(n,:).*exp(1i*phi(n,:)); %back to complex
    end

end


%end of Mega loop
end

%average data
for n =1:size(FFTComplex,2)
    AvgSpectra(n) = mean(FFTComplex(:,n));
    AvgRef(n) = mean(FFTComplexRef(:,n));
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

Spectra = NormalisedSpectraAvg;
SpectraRaw = FFTComplexAvg;

end