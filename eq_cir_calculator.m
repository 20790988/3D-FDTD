epsilon_0 = 8.8542e-12;
mu_0 = 1.2566e-6;
h = 0.3e-6;
W = 4e-6;
e_r=7.29;

%assume W/h >1

K = (h/W)*(W/h+2.42-0.44*(h/W)+(1-h/W).^6);

C = W*K*((9.5*e_r+1.25)*W/h+5.2*e_r+7)
%pF

L = 100*h*(4*sqrt(W/h)-4.21)
%nH