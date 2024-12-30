
% get n-th demodulated final output, final function used 
function eff_pol_n = ML_n_3layers( H_0, Full_Tapping_Amplitude, Tapping_Frequency, tip_r, length, g, W_0, W_1, epsilon1, epsilon2, epsilon3, epsilon4, d2, d3, n)
Lgr_n = 12;
beta12 = beta_ij(epsilon1, epsilon2);
beta23 = beta_ij(epsilon2, epsilon3);
beta34 = beta_ij(epsilon3, epsilon4);
eff_pol_n = Tapping_Frequency*integral(@(time) MLFDM_itg_3layers(time, Lgr_n, H_0, Full_Tapping_Amplitude, Tapping_Frequency, tip_r, length, g, W_0, W_1,beta12, beta23, beta34, d2,d3,n), -1/(2*Tapping_Frequency), 1/(2*Tapping_Frequency));
end


% % all units in SI
% 
% n= 3; % demodulation order
% % input sSNOM parameters
% R = 30e-9; % m, tip radius
% L = 200e-9; % m, tip half length 
% g = 0.7*exp(1i*0.06); 
% H_0 = 10e-9; % in-contact tip height
% TA = 35e-9; % m, tapping amplitude (half of full oscillation)
% tip_f = 230000; % Hz, tip tapping frequency
% 
% W_0 = 1.31*R*L/(L+2*R);
% W_1 = 0.5*R;
% T = 1/tip_f;
% 
% % input layer dielectric functions and thickness
% % air
% epsilon1 = 1; % air
% % 1st layer
% epsilon2 = 2.5; 
% d2 = 10e-9; % m
% % 2nd layer
% epsilon3 = 10;
% d3 = 59e-9; % m 
% % semi-infinite substrate
% epsilon4 = 11.8; % Si subsutarte 
% 
% 












function betaij = beta_ij(epsilon_i, epsilon_j)
betaij =  (epsilon_j - epsilon_i)./(epsilon_i + epsilon_j);
end 

% calculate potential and diff_potential with Gauss-Laguerre method 
function potential_at_0 = potential_z (beta12, beta23, t1, z0, Lgr_n)
% get root of n-th Laguerre polynominal
syms x
xi = roots(sym2poly(laguerreL(Lgr_n,x)));
w_i = xi./((Lgr_n+1).^2 .* laguerreL(Lgr_n+1,xi).^2) ; % weight function
f_xi = (beta12 + beta23.*exp(- xi*t1./z0))./(1+beta12.*beta23.*exp(-xi.*t1./z0)).*(1./(2*z0));
% calculate sum
potential_at_0 = sum(w_i.*f_xi); % final output
end 
function diff_potential_at_0 = diff_potential_z(beta12, beta23, t1, z0, Lgr_n)
syms x
xi = roots(sym2poly(laguerreL(Lgr_n,x)));
w_i = xi./((Lgr_n+1).^2 .* laguerreL(Lgr_n+1,xi).^2) ; % weight function
f_xi = (beta12 + beta23.*exp(- xi*t1./z0))./(1+beta12.*beta23.*exp(-xi*t1./z0)).*(xi./(4*z0.^2));
diff_potential_at_0 = sum(w_i.*f_xi);
end 

% from H_t get beta0, beta1, X0, X1
function [beta_0, beta_1, X0, X1] = get_beta_X0_X1(Lgr_n, z0, z1 , beta12, t1, beta23)
beta_0 =  potential_z (beta12, beta23, t1, z0, Lgr_n).^2 ./  diff_potential_z(beta12, beta23, t1, z0, Lgr_n);
beta_1 =  potential_z (beta12, beta23, t1, z1, Lgr_n).^2 ./  diff_potential_z(beta12, beta23, t1, z1, Lgr_n);
X0 = abs(potential_z (beta12, beta23, t1, z0, Lgr_n) ./  diff_potential_z(beta12, beta23, t1, z0, Lgr_n)) - z0;
X1 = abs(potential_z (beta12, beta23, t1, z1, Lgr_n) ./  diff_potential_z(beta12, beta23, t1, z1, Lgr_n)) - z1;
end

% ML FDM original form 
function eff_pol_to_t_int_MLFDM = MLFDM(Lgr_n, H_t, tip_r, length, g, W_0, W_1, beta12,t1, beta23)
z0 = W_0 + H_t;
z1 = W_1 + H_t;
[beta_0, beta_1, X0, X1] = get_beta_X0_X1(Lgr_n, z0, z1 , beta12, t1, beta23);
f0 = (g-((tip_r+H_t+X0))./(2*length)).*log((4.*length)./(2.*H_t+tip_r+2.*X0))./(log(4.*length./tip_r));
f1 = (g-((tip_r+H_t+X1))./(2*length)).*log((4.*length)./(2.*H_t+tip_r+2.*X1))./(log(4.*length./tip_r));
eff_pol_to_t_int_MLFDM = (2+  beta_0.*f0 ./ ( 1 - beta_1.*f1)) ;
end
% MLFDM with time dependence, intergration function
function itg_fun = MLFDM_itg(time, Lgr_n, H_0, Full_Tapping_Amplitude, Tapping_Frequency, tip_r, length, g, W_0, W_1, beta12,t1, beta23, n)
H_t = (H_0 + Full_Tapping_Amplitude./2) + (Full_Tapping_Amplitude./2) .* sin(2.*pi.*Tapping_Frequency.*time);
itg_fun = MLFDM(Lgr_n, H_t, tip_r, length, g, W_0, W_1, beta12,t1, beta23).*exp(-1i.*n.*2.*pi.*Tapping_Frequency.*time);
end 
% get n-th demodulated final output, final function used 
function eff_pol_n = ML_n(Lgr_n, H_0, Full_Tapping_Amplitude, Tapping_Frequency, tip_r, tip_length, g, W_0, W_1, epsilon1, epsilon2, epsilon3, t1, n)
beta12 = beta_ij(epsilon1, epsilon2);
beta23 = beta_ij(epsilon2, epsilon3);
eff_pol_n = Tapping_Frequency*integral(@(time) MLFDM_itg(time, Lgr_n, H_0, Full_Tapping_Amplitude, Tapping_Frequency, tip_r, tip_length, g, W_0, W_1, beta12,t1, beta23, n), -1/(2*Tapping_Frequency), 1/(2*Tapping_Frequency));
end


% ML 3 layers
% get beta0, beta1, X0, X1
function [beta_0, beta_1, X_0, X_1] = get_beta_X0_X1_3layers(Lgr_n, z0, z1 , beta12, beta23, beta34, d2,d3)
syms x
xi = roots(sym2poly(laguerreL(Lgr_n,x)));
w_i = xi./((Lgr_n+1).^2 .* laguerreL(Lgr_n+1,xi).^2) ; % weight function

beta2_34 = (beta23 + beta34.*exp(- xi*d3./z0))./(1+beta23.*beta34.*exp(-xi.*d3./z0));
f_xi_potential = (beta12 + beta2_34.*exp(- xi*d2./z0))./(1+beta12.*beta2_34.*exp(-xi.*d2./z0)).*(1./(2*z0));
potential_at_0 = sum(w_i.*f_xi_potential); % final output
f_xi_diff_potential = (beta12 + beta2_34.*exp(- xi*d2./z0))./(1+beta12.*beta2_34.*exp(-xi*d2./z0)).*(xi./(4*z0.^2));
diff_potential_at_0 = sum(w_i.*f_xi_diff_potential);
% get beta_0 and X_0
beta_0 = potential_at_0.^2./diff_potential_at_0;
X_0 =  abs(potential_at_0 ./diff_potential_at_0) -z0;
% 
beta2_34 = (beta23 + beta34.*exp(- xi*d3./z1))./(1+beta23.*beta34.*exp(-xi.*d3./z1));
f_xi_potential = (beta12 + beta2_34.*exp(- xi*d2./z1))./(1+beta12.*beta2_34.*exp(-xi.*d2./z1)).*(1./(2*z1));
potential_at_0 = sum(w_i.*f_xi_potential);
f_xi_diff_potential = (beta12 + beta2_34.*exp(- xi*d2./z1))./(1+beta12.*beta2_34.*exp(-xi*d2./z1)).*(xi./(4*z1.^2));
diff_potential_at_0 = sum(w_i.*f_xi_diff_potential);
% get beta_0 and X_0
beta_1 = potential_at_0.^2./diff_potential_at_0;
X_1 =  abs(potential_at_0 ./diff_potential_at_0) -z1;
end

% ML FDM original form 
function eff_pol_to_t_int_MLFDM = MLFDM_3layers(Lgr_n, H_t, tip_r, length, g, W_0, W_1, beta12, beta23, beta34, d2,d3)
z0 = W_0 + H_t;
z1 = W_1 + H_t;
[beta_0, beta_1, X0, X1] = get_beta_X0_X1_3layers(Lgr_n, z0, z1 , beta12, beta23, beta34, d2,d3);
f0 = (g-((tip_r+H_t+X0))./(2*length)).*log((4.*length)./(2.*H_t+tip_r+2.*X0))./(log(4.*length./tip_r));
f1 = (g-((tip_r+H_t+X1))./(2*length)).*log((4.*length)./(2.*H_t+tip_r+2.*X1))./(log(4.*length./tip_r));
eff_pol_to_t_int_MLFDM = (2+  beta_0.*f0 ./ ( 1 - beta_1.*f1)) ;
end
% MLFDM with time dependence, intergration function
function itg_fun = MLFDM_itg_3layers(time, Lgr_n, H_0, Full_Tapping_Amplitude, Tapping_Frequency, tip_r, length, g, W_0, W_1,beta12, beta23, beta34, d2,d3,n)
H_t = (H_0 + Full_Tapping_Amplitude./2) + (Full_Tapping_Amplitude./2) .* sin(2.*pi.*Tapping_Frequency.*time);
itg_fun = MLFDM_3layers(Lgr_n, H_t, tip_r, length, g, W_0, W_1, beta12, beta23, beta34, d2,d3).*exp(-1i.*n.*2.*pi.*Tapping_Frequency.*time);
end 

