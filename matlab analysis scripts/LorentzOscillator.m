function Output = LorentzOscillator(WavenumberRange,EpsilonInf,Amplitude,ResonanceWavenumber,Width)

% Speed of light in cm/s
c = 3e10;

%convert wavenumber to angular frequency
Omega = 2 * pi * c * WavenumberRange;
Omega_0 = 2 * pi * c * ResonanceWavenumber;
%function from : https://journals.aps.org/rmp/pdf/10.1103/RevModPhys.83.543

Output = EpsilonInf + Amplitude./(Omega_0^2 - Omega.^2 - 1i*Omega*Width);

end