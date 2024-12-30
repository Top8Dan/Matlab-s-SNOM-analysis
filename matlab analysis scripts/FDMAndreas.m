function A_eff_n = FDMAndreas(Beta, Radius, Tmin, Tamp, freq, n, t, g, L)

z = Tmin + Tamp*(1+cos(2*pi*freq*t));

f0 = (g - ((Radius + 2*z + 1.31 * Radius)/(2*L)))*...
((log((4*L)/(Radius + 4*z + 2*1.31*Radius))/(log(4*L/Radius))));

f1 = (g - ((Radius + 2*z + 0.5 * Radius)/(2*L)))*...
((log((4*L)/(Radius + 4*z + 2*0.5*Radius))/(log(4*L/Radius))));

% A_eff = freq.*(exp(-1i*n*2.0*pi*freq*t)).* (0.5*((Beta * f0)/(1-Beta*f1))+1); 
A_eff = 0.5*((Beta * f0)/(1-Beta*f1))+1; 
A_eff_n = exp(-1i*n*2*pi*freq*t).* A_eff; 


% A_eff = freq.*(exp(-1i*n*2.0*pi*freq*t)).* 0.5.* ((Beta*((g - ((Radius + 2.*(Tmin + (Tavg.*(1+cos(freq).*t))) + 0.5 .* Radius)/(2.*L))).*...
% ((log((4.*L)./(Radius + 4.*(Tmin + (Tavg.*(1+cos(freq).*t))) + Radius))./(log(4.*L./Radius)))))))...
% /(1-Beta.*((g - ((Radius + 2.*(Tmin + (Tavg.*(1+cos(freq).*t))) + 1.31 .* Radius)./(2.*L))).*...
% ((log((4.*L)./(Radius + 4.*(Tmin + (Tavg.*(1+cos(freq).*t))) + 2.*1.31.*Radius))./(log(4.*L./Radius)))))) + 1;




end