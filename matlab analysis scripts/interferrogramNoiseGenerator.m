% Step 1: Define the frequency components of the spectrum
nu_center = 1000; % Central frequency (wavenumber) in cm^-1
delta_nu = 1000; % Frequency range around the center in cm^-1
N = 100; % Number of data points

% Create the frequency vector
nu = linspace(nu_center - delta_nu/2, nu_center + delta_nu/2, N);

% Step 2: Simulate the spectrum as a sum of Gaussian peaks
sigma = 50; % Standard deviation of the Gaussian peaks
S_nu = exp(-((nu - nu_center).^2) / (2 * sigma^2)); % Single Gaussian spectrum
S_nu = S_nu + exp(-((nu - (nu_center + 100)).^2) / (2 * sigma^2)); % Add another peak
S_nu = S_nu + exp(-((nu - (nu_center - 150)).^2) / (2 * sigma^2)); % Add another peak

S_nu = DrudeLorentzDielectricEH(nu,1.3,1172,90,1230,80) + DrudeLorentzDielectricEH(nu,1.3,972,90,1030,80); %SiO2;

% Step 3: Add noise to the spectrum
noise_level_spectrum = 0.0; % Adjust this value to control noise level in spectrum
S_nu_noisy = S_nu + noise_level_spectrum * randn(size(S_nu));

% Step 4: Compute the interferogram using inverse Fourier transform
I_x_noisy = ifft(S_nu_noisy); % Inverse Fourier Transform
I_x_noisy = fftshift(I_x_noisy); % Shift the zero frequency component to the center

% Step 5: Add noise to the interferogram
noise_level_interferogram = 0.02; % Adjust this value to control noise level in interferogram
I_x_noisy_real = real(I_x_noisy) + noise_level_interferogram * randn(size(I_x_noisy));

% Step 6: Generate the corresponding x values (Optical Path Difference)
x = fftshift(ifftshift((0:N-1) - floor(N/2)) / (N * (nu(2) - nu(1))));

% Step 7: Plot the noisy spectrum and noisy interferogram packet
figure;

subplot(2, 1, 1);
plot(nu, S_nu_noisy);
title('Noisy Simulated Spectrum with Multiple Peaks');
xlabel('Wavenumber (\nu) in cm^{-1}');
ylabel('Intensity');

subplot(2, 1, 2);
plot(x, I_x_noisy_real);
title('Noisy Interferogram Packet');
xlabel('Optical Path Difference (x) in cm');
ylabel('Intensity');


% Step 8: Save the spectrum and interferogram to text files

% Save the spectrum
spectrum_data = [nu.' S_nu_noisy.']; % Create a matrix with wavenumber and intensity
save('C:\Users\k68558dj\OneDrive - The University of Manchester\DataFiles\interNoiseSimulation\No Noise\simulated_spectrum.txt', 'spectrum_data', '-ascii', '-tabs');

% Save the interferogram
interferogram_data = [x.' I_x_noisy_real.']; % Create a matrix with OPD and intensity
save('C:\Users\k68558dj\OneDrive - The University of Manchester\DataFiles\interNoiseSimulation\No Noise\simulated_interferogram.txt', 'interferogram_data', '-ascii', '-tabs');


