function [parametersOut, resnorm, residuals, exitflag, output, lambda, jacobian] = FDMMultilayerFit(Wavenumber,Complex,Tmin,Tamp,freq,Harmonic,g,R,L,W_0,W_1,Es_Ref1,Es_Ref2,Es_Ref3,Es_Ref4,dRef2,dRef3,Es_1,Es_2,Es_3,Es_4,d2,d3)

T = 1/freq;
global Es iteration

%iteration tracker

iteration = 0;


%inital guess
Es = Es_2;

initial_guess = (Es);


%options for lsqnonlin fit.
options = optimset('MaxIter',100000,'MaxFunEvals', 100000, 'FinDiffType', 'central', 'TolFUn', 1e-30, 'TolX', 1e-18);

%this is the function to minimise, the difference between the FDM and the
%measured data
function global_error = error_return(params)

    Es = real(params) + abs(imag(params))*1i; % make imaginary part always positive, without slowing down fit.
        
% Set value dielectric to be fit as Es

    %material calculation        
    FDMValue = ML_n_3layers(Tmin, Tamp, freq, R, L, g, W_0, W_1, Es_1, Es, Es_3, Es_4, d2, d3, Harmonic);
    %reference calculation
  
    RefValue = ML_n_3layers(Tmin, Tamp, freq, R, L, g, W_0, W_1, Es_Ref1, Es_Ref2, Es_Ref3, Es_Ref4, dRef2, dRef3, Harmonic);
    

    %normalise
    FDMValue = FDMValue/RefValue;

    %compare real value to calculated
    Error_vector = (Complex - FDMValue);
    global_error = Error_vector;
    
    Es;
    iteration = iteration+1;
end

%run least square fit
[parametersOut, resnorm, residuals, exitflag, output, lambda, jacobian] = lsqnonlin(@error_return, initial_guess, -inf, inf,options);

iteration;

end
