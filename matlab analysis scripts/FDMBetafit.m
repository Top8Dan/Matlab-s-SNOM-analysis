function [parametersOut, resnorm, residuals, exitflag, output, lambda, jacobian] = FDMBetafit(Wavenumber,Complex,Tmin,Tamp,freq,Harmonic,g,R,L,E_s)

T = 1/freq;
global Beta iteration

%iteration tracker

iteration = 0;

%inital guess
Beta = 1.2 + 0.05i;

initial_guess = Beta;

%options for lsqnonlin fit.
options = optimset('MaxIter',10000,'MaxFunEvals', 10000, 'FinDiffType', 'central', 'TolFUn', 1e-28, 'TolX', 1e-12);

%this is the function to minimise, the difference between the FDM and the
%measured data
function global_error = error_return(params)

    Beta = params;
        
    %material calculation        
    fiteqn = @(t) FDMHauer(Beta,R,Tmin,Tamp,freq,Harmonic,t,g,L);
    FDMValue = integral(fiteqn,-T/2,T/2,'arrayvalued', true);

    %reference calculation
  
    RefBeta = (E_s - 1)./(E_s +1);
    fiteqnRef = @(t) FDMHauer(RefBeta,R,Tmin,Tamp,freq,Harmonic,t,g,L);
    RefValue = integral(fiteqnRef, -T/2,T/2,'arrayvalued', true);
    

    %normalise
    FDMValue = FDMValue/RefValue;

    %compare real value to calculated
    Error_vector = (Complex - FDMValue);
    global_error = Error_vector;
    
    
    iteration = iteration+1;
end

%run least square fit
[parametersOut, resnorm, residuals, exitflag, output, lambda, jacobian] = lsqnonlin(@error_return, initial_guess, [-inf], [inf],options);

iteration;

end
