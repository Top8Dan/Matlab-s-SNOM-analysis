% Parameters
me = 9.109e-31;
n_cm3 = 6e18; % Carrier density in cm^-3
mu_cm2_per_Vs = 150; % Mobility in cm^2/Vs
e = 1.602e-19; % Elementary charge in C
m = 0.475*me; % Effective mass in kg
c = 3e10; % Speed of light in cm/s
epsilon_0 = 8.854e-12; % Permittivity of free space in F/m
epsilon_inf = 10; % High-frequency permittivity

% Convert carrier density to m^-3
n = n_cm3 * 1e6; % Carrier density in m^-3

% Calculate scattering time tau from mobility
tau = (mu_cm2_per_Vs * 1e-4) * m / e; % Convert cm^2/Vs to m^2/Vs for mobility

% Define wavenumber range (in cm^-1)
k_min_cm = 1e2; % Minimum wavenumber in cm^-1
k_max_cm = 1e4; % Maximum wavenumber in cm^-1
k_cm = linspace(k_min_cm, k_max_cm, 1000); % Wavenumber array in cm^-1

% Convert wavenumber from cm^-1 to m^-1
k = k_cm * 100; % Convert to m^-1

% Calculate angular frequency from wavenumber
omega = c * k; % Angular frequency in rad/s

% Calculate complex dielectric function as a function of wavenumber
numerator = n * e^2 / (epsilon_0 * m);
denominator = omega .* (omega + 1i / tau);
epsilon_k = epsilon_inf - numerator ./ denominator;

% Separate real and imaginary parts of dielectric function
epsilon_real = real(epsilon_k);
epsilon_imag = imag(epsilon_k);

% Plot the results
figure;
subplot(2,1,1);
plot(k_cm, epsilon_real, 'b-', 'LineWidth', 2);
xlabel('Wavenumber (cm^{-1})', 'FontSize', 14);
ylabel('Real Permittivity', 'FontSize', 14);
title('Real Part of Dielectric Function vs. Wavenumber', 'FontSize', 16);
grid on;

subplot(2,1,2);
plot(k_cm, epsilon_imag, 'r-', 'LineWidth', 2);
xlabel('Wavenumber (cm^{-1})', 'FontSize', 14);
ylabel('Imaginary Permittivity', 'FontSize', 14);
title('Imaginary Part of Dielectric Function vs. Wavenumber', 'FontSize', 16);
grid on;