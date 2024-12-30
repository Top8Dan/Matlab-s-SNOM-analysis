close all
clear all

% FDM parameters
R = 35e-9; %tip radius
Tmin = 1e-9; %tip minimum distance
Tamp = 50e-9; %tip amplitude position
freq = 1; %tip frequency
T = 1/freq; %tip period
g = 0.7*exp(0.06*1i);
L = 300e-9;
eV = 1.6e-19;
Wcm = 8065;

Harmonic = 3;


% Frequency range
omega = linspace(0, 3200, 1000);

%normalisation dielectric
eps_Si = 11.7 + 0.1i;

%  Bi2Se3 paper
% lorentzian component
eps_inf = 18.2 + 0.6i;
omega_plasma = 0.75*Wcm;
omega_lorentz = 0.5e-3 *Wcm;
gamma_lorentz = 8.45*Wcm;

eps_lorentz1 = 1 + omega_plasma^2 *(0.8./(omega_lorentz^2 - omega.^2 - 1i*omega*gamma_lorentz));

eps_inf = 18.2 + 0.6i;
omega_plasma = 6.89*Wcm;
omega_lorentz = 1.39*Wcm;
gamma_lorentz = 1.56*Wcm;

eps_lorentz2 = 1 + omega_plasma^2 *(0.8./(omega_lorentz^2 - omega.^2 - 1i*omega*gamma_lorentz));

eps_lorentz = eps_lorentz1+eps_lorentz2;

%% Drude component
omega_plasma = 0.48*Wcm;
gamma_drude = 1.24e-6 *Wcm;

eps_drude = eps_inf - ((omega_plasma^2)./(omega.^2 + 1i*omega.*gamma_drude));

eps_total = eps_lorentz + eps_drude;

% lorentzian component
eps_inf = 18.2 + 0.6i;
omega_plasma = 0.75*Wcm;
omega_lorentz = 0.5e-3 *Wcm;
gamma_lorentz = 8.45*Wcm;

eps_lorentz1 = 1 + omega_plasma^2 *(0.8./(omega_lorentz^2 - omega.^2 - 1i*omega*gamma_lorentz));


% Convert to Beta
Beta = (eps_lorentz - 1)./(eps_lorentz +1);
beta_Si = (eps_Si - 1)./(eps_Si +1);

% FDM eqn

for i =1:size(omega,2)
    FDMEqn = @(t) FDMHauer(Beta(i),R,Tmin,Tamp,freq,Harmonic,t,g,L);
    FDMOutput(i)= integral(FDMEqn, -T/2,T/2,'arrayvalued', true);
end

FDMEqn_Si = @(t) FDMHauer(beta_Si,R,Tmin,Tamp,freq,Harmonic,t,g,L);
FDM_Si = integral(FDMEqn_Si, -T/2,T/2,'arrayvalued', true);

% Plotting
figure(1)
subplot(2,1,1);
plot(omega, real(eps_total), 'b', 'LineWidth', 2);
xlabel('Wavenumber (cm^{-1})');
ylabel('Real Permittivity');

subplot(2,1,2);
plot(omega, imag(eps_total), 'r', 'LineWidth', 2);

xlabel('Wavenumber (cm^{-1})');
ylabel('Imaginary Permittivity');

figure(2)
subplot(2,1,1)
plot(omega,real(Beta))
title('Beta plot')
subplot(2,1,2)
plot(omega,imag(Beta))

figure(3)
subplot(2,1,1)
plot(omega,abs(FDMOutput/FDM_Si), 'b', 'LineWidth', 2)
title('Spectra Output')
ylabel('abs')
xlabel('Wavenumber (cm^{-1})')

subplot(2,1,2)
plot(omega,angle(FDMOutput/FDM_Si), 'r', 'LineWidth', 2)
ylabel('angle')
xlabel('Wavenumber (cm^{-1})')


