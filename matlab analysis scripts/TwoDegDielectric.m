function Epsilon = TwoDegDielectric(WavenumberRange,Epsilon_Inf,CarrierDensity,EffectiveMass)

Epsilon_0 = 8.854e-12; % Permittivity of free space in F/m
c = 3e10; % Speed of light in cm/s
e = 1.602e-19; % Elementary charge in coulombs
m0 = 9.109e-31; % Electron rest mass in kg

EffectiveMass = EffectiveMass*m0;

% Initialize arrays to store results
Epsilon = zeros(size(WavenumberRange));


for idx = 1:length(WavenumberRange)
    % Convert wavenumber to angular frequency
    omega = 2 * pi * c * WavenumberRange(idx);
    
    % Calculate sigma_2DEG
    sigma_2DEG = (e^2 * CarrierDensity / EffectiveMass) * (1i / omega);
    
    % Calculate the dielectric function epsilon
    Epsilon(idx) = Epsilon_Inf + (1i * sigma_2DEG) / (Epsilon_0 * omega);

end

end