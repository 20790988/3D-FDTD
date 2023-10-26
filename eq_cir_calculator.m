epsilon_0 = 8.8542e-12;
mu_0 = 1.2566e-6;

%microstrip

% h = 0.3e-6;
% W = 4e-6;
% e_r=7.29;
% % e_r=4.6;
% 
% %assume W/h >1
% 
% K = (h/W)*(W/h+2.42-0.44*(h/W)+(1-h/W).^6);
% 
% C = W*K*((9.5*e_r+1.25)*W/h+5.2*e_r+7)*1e-12
% %pF
% 
% L = 100*h*(4*sqrt(W/h)-4.21)*1e-9
% %nH

% stripline

w = 4.4e-6;
b = 0.6e-6;
t = 0.2e-6;
e_r = 4.6;
% e_r = 9.76;

f = 1e0:1e9:500e9;

c = 1/sqrt(epsilon_0*mu_0);

% right angle

D = w+2*b/pi*log(2);

gam_g = c./(sqrt(e_r).*f);

X_a = D*(1.756+ 4*((D./gam_g).^2))./gam_g;

X_b = D*(0.0725-0.159*((gam_g./D).^2))./gam_g;

L_1 = X_a./(2*pi*f);
C_1 = -1./(X_b*2*pi.*f);

% arb angle

D = w + 2*b/pi*log(2)+t/pi*(1-log(2*t/b));

x = 0.5*(1+90/180);

psi = 0.5223*log(x)+0.394;

X_a = 2*D*(psi+1.9635+1./x)./gam_g;

X_b = (-gam_g*cot(pi/2/2))/(2*pi*D);


L_2 = X_a./(2*pi*f);
C_2 = -1./(X_b*2*pi.*f);


figure(1)
hold on
plot(f/1e9,L_1./1e-12,'m','LineWidth',2);
xlim([0 500]);
xlabel('freq (GHz)');
ylabel('L (nH)');
plot(f/1e9,L_2./1e-12,'k','LineWidth',2);
grid on
hold off
legend('Bahl & Garg 1978','Pramanick & Bhartia 2016','Location','southeast')
% ylim([0 0.33]);


figure(2);
hold on
plot(f/1e9,C_1/1e-15,'m','LineWidth',2);
xlim([0 500]);

xlabel('freq (GHz)');
ylabel('C (fF)');
plot(f/1e9,C_2/1e-15,'k','LineWidth',2);
hold off
% ylim([0 3.2]);
grid on
% legend('FDTD','Formula from Javadzadeh et al','Location','southeast')
