function epsilon_k = DrudeFunction(WavenumberRange,Epsilon_inf,CarrierDensity,Mobility,EffectiveMass)

% Parameters - https://pubs.acs.org/doi/10.1021/acs.nanolett.8b03008
me = 9.109e-31;
n_cm3 = CarrierDensity; % Carrier density in cm^-3
mu_cm2_per_Vs = Mobility; % Mobility in cm^2/Vs
e = 1.602e-19; % Elementary charge in C
m = EffectiveMass*me; % Effective mass in kg
c = 3e10; % Speed of light in cm/s
epsilon_0 = 8.854e-12; % Permittivity of free space in F/m
epsilon_inf = Epsilon_inf; % High-frequency permittivity
k_cm = WavenumberRange; %wavenumber range in cm-1

% Convert carrier density to m^-3
n = n_cm3 * 1e6; % Carrier density in m^-3

% Calculate scattering time tau from mobility
tau = (mu_cm2_per_Vs * 1e-4) * m / e; % Convert cm^2/Vs to m^2/Vs for mobility


% Convert wavenumber from cm^-1 to m^-1
k = k_cm * 100; % Convert to m^-1

% Calculate angular frequency from wavenumber
omega = c * k; % Angular frequency in rad/s

% Calculate complex dielectric function as a function of wavenumber
numerator = n * e^2 / (epsilon_0 * m);
denominator = omega .* (omega + 1i / tau);
epsilon_k = epsilon_inf - numerator ./ denominator;

end