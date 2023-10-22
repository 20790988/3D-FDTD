%In microstrip

epsilon_0 = 8.8542e-12;
mu_0 = 1.2566e-6;


h = 0.3e-6;
w = 4e-6;
e_r = 4.6;

c = 1/sqrt(epsilon_0*mu_0);

fT1 = num2eng((c/(2*pi*h))*sqrt(2/(e_r-1))*atan(e_r))

fCT = num2eng(c/(sqrt(e_r)*(2*w+0.8*h)))

fT2 = num2eng(c/(4*h*sqrt(e_r-1)))

fT4 = num2eng(c/(2*h*sqrt(e_r)))

