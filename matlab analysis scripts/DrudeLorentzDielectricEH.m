function Output = DrudeLorentzDielectricEH(wavenumberRange,Es_inf,OmegaLv,yLv,OmegaTv,yTv)

% eqns from Vibrational modes in amorphous silicon dioxide Marta Klanjs\ek Gunde 
%% Oscilator parameters

%% Handling
freqX = wavenumberRange;

for z = 1:size(freqX,2)
    IRFreq = freqX(z);
    % Es(z) = Es_inf*((OmegaTv^2 - IRFreq^2 + 1i*IRFreq*yTv)/(OmegaLv^2 - IRFreq^2 + 1i*IRFreq*yLv));
    TH = OmegaTv^2 - IRFreq^2 - 1i * IRFreq * yTv;
    BH = OmegaLv^2 - IRFreq^2 - 1i * IRFreq * yLv;
    Es(z) = Es_inf * (TH / BH);  

end

Output = Es;


end