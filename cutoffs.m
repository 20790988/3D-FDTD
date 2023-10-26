%In microstrip

epsilon_0 = 8.8542e-12;
mu_0 = 1.2566e-6;
c = 1/sqrt(epsilon_0*mu_0);


% h = 0.3e-6;
% w = 4e-6;
% e_r = 4.6;
% 
% 
% fT1 = num2eng((c/(2*pi*h))*sqrt(2/(e_r-1))*atan(e_r))
% 
% fCT = num2eng(c/(sqrt(e_r)*(2*w+0.8*h)))
% 
% fT2 = num2eng(c/(4*h*sqrt(e_r-1)))
% 
% fT4 = num2eng(c/(2*h*sqrt(e_r)))

%in stripline

w = 4.4e-6;
b = 0.6e-6;
e_r = 4.6;


%convert w, b to cm 
w_ = w*100;
b_ = b*100;

f_T= num2eng((15/(b_*sqrt(e_r)*(w_/b_+pi/4)))*1e9)


%b/lamda_c = 4d/b
x_1 = [0.0 0.2 0.3 0.35 0.4 0.45 0.5];
f_1 = [0.882 0.917 0.968 1.016 1.070 1.180 1.586];

x_2 = 0.001:0.0001:0.5;
f_2 = 1./(sqrt(e_r)*x_2)-2*w/b;

figure(1);
hold on
plot(x_1,f_1,'-o','LineWidth',2);
plot(x_2,f_2,'LineWidth',2);

hold off
grid on

xlabel('$\frac{b}{\lambda_c}$','Interpreter','latex','fontsize',16)
ylabel('$\frac{4d}{b}$','Interpreter','latex','fontsize',16)
legend('Values from Matthaei et al','Curve for RSFQ line')
ylim([-15 5])
x=0.03;

lambda_c = (1/x)*b;

f_c = num2eng(c/(lambda_c))







