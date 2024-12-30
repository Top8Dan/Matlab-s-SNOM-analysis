clear all
close all
%Lorentzian curve
Start = 0;
End = 3000;

x0 = 1200;
gamma = 300;

x = linspace(Start,End,3200);

    L= (gamma/(2*pi)) ./ ((x - x0).^2 + (gamma/2)^2);
    LImag = -(gamma/(2*pi)) * (gamma/2) ./ ((x - x0).^2 + (gamma/2)^2);

    LComplex = complex(L,LImag);

hold on
plot(x,real(L))
plot(x,(LImag))
hold off

figure(2)
plot(x,abs(LComplex))

figure(3)
plot(x,angle(LComplex))
