
function [parametersOut, resnorm, residuals, exitflag, output, lambda, jacobian] = FDMfit(XData,Amplitude,E_s,E_t,Tmin,freq,Harmonic,g)

%calculate period
T = 1/freq;
Beta = (E_s - 1)./(E_s +1);

%global variables
global Tamp Radius L iteration

%fitting tracking
iteration = 0;

%inital guesses

%these three variables are scaled e-9 in globalerror function
L = 388.608;
Radius = 60;
Tamp = 50;

initial_guess = [L Radius Tamp];

%options for lsqnonlin fit.
options = optimset('MaxIter',100000,'MaxFunEvals', 100000, 'FinDiffType', 'central', 'TolFUn', 1e-20, 'TolX', 1e-20);

%this is the function to minimise, the difference between the FDM and the
%measured data
function global_error = error_return(params)

    L = params(1)*1e-9;
    Radius = params(2)*1e-9;
    Tamp = params(3) *1e-9;
    global_error = zeros(size(XData));
    
    
        for c = 1:size(XData,1)

            Tmin = XData(c)*1e-9;


            fiteqn = @(t) FDMHauer(Beta,Radius,Tmin,Tamp,freq,Harmonic,t,g,L);
 
            FDMCurve(c) = integral(fiteqn,-T/2,T/2,'arrayvalued', true);
        end
    
    FDMCurve = abs(FDMCurve/FDMCurve(1));
    FDMCurve = (FDMCurve).';
    
    Error_vector = (Amplitude - FDMCurve);
    global_error = global_error+Error_vector;
    
    
    iteration = iteration+1
end

%use matlab least square fit
[parametersOut, resnorm, residuals, exitflag, output, lambda, jacobian] = lsqnonlin(@error_return, initial_guess, [0 10 20], [800 80 100],options);



parametersOut(1) = parametersOut(1)*1e-9;
parametersOut(2) = parametersOut(2)*1e-9;
parametersOut(3) = parametersOut(3)*1e-9;


L = parametersOut(1);
Radius = parametersOut(2);
Tamp = parametersOut(3);

end
