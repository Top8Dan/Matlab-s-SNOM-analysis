function [Fun] = Dielectric_Au(w)

% Ordal et al., Applied Optics, Vol. 22, No. 7, 1983

 epsilon_inf = 1.0;
 wp = 7.25e6; % plasma frequency
 wt = 2.16e4; % damping frequency
 Fun = epsilon_inf - ((wp^2)./(w.^2 + 1i*w.*wt));

%  Fun = (epsilonAu-1)./(epsilonAu+1);

end