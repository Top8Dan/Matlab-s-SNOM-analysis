function epsilon = DiracDielectric(WavenumberRange,N_Dirac,FermiVelocity)

% Constants
e = 1.60217662e-19;     % Electron charge, C
hbar = 1.0545718e-34;   % Reduced Planck's constant, J*s
k_B = 1.38064852e-23;   % Boltzmann constant, J/K
%v_F = 5e5;              % Fermi velocity, m/s
v_F = FermiVelocity;
T = 300;                % Temperature, K
kT = 25e-3 * e;         % kT at 300 K in Joules
epsilon_0 = 8.85418782e-12; % Permittivity of free space, F/m
c = 3e10;               %speed of light cm

% Parameters
%n_Dirac = 1e12;         % Surface carrier concentration, m^-2 (example value)
EF = hbar * v_F * sqrt(4 * pi * N_Dirac); % Fermi energy, J
gamma_Dirac = 0;        % Carrier relaxation rate, assumed to be 0

% Frequency range
omega = 2 * pi * c * WavenumberRange;

% Surface conductivity calculation
sigma_Dirac = (e^2 * kT * log(2 * cosh(EF / (2 * kT))) / (2 * hbar^2 * pi)) * 1i ./ (omega + 1i * gamma_Dirac);

% Dielectric function calculation
epsilon = 1 + 1i * sigma_Dirac ./ (epsilon_0 * omega);

% % Plotting the real and imaginary parts of the dielectric function
% figure;
% subplot(2,1,1);
% plot(WavenumberRange, real(epsilon), 'LineWidth', 1.5);
% xlabel('Wavenumber (cm-1)');
% ylabel('Re(\epsilon)');
% title('Real Part of Dielectric Function \epsilon(\omega)');
% 
% subplot(2,1,2);
% plot(WavenumberRange, imag(epsilon), 'LineWidth', 1.5);
% xlabel('Wavenumber (cm-1)');
% ylabel('Im(\epsilon)');
% title('Imaginary Part of Dielectric Function \epsilon(\omega)');

% Displaying the plot
grid on;
set(gca, 'FontSize', 12);

end