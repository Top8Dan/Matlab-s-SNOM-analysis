function Sigma = TwoDegConductivity(WavenumberRange,CarrierDensity,EffectiveMass)

% Constants
e = 1.602e-19; % Elementary charge in coulombs
m0 = 9.109e-31; % Electron rest mass in kg
c = 3e10; % Speed of light in cm/s
% CarrierDensity in cm^-2
m = EffectiveMass*m0;

omega = 2 * pi * c * WavenumberRange;

% Conductivity calculation
Sigma = ((e^2 * CarrierDensity) / (m)) * (1i ./ omega);

end