clear
close all
load('sparams.mat');

% z11 = 0;
% z12 = 0;
% z21 = z12;
% z22 = 0;
% z0 = 0;
% 
% del_z = (z11 +z0).*(z22+ z0) -z12.*z21;
% 
% s11 = ((z11-z0).*(z22+z0)-z12*z21)./del_z;
% s12 = (2*z12.*z0)./del_z;
% s21 = s12;
% s22 = ((z11+z0).*(z22-z0)-z12*z21)./del_z;


%or


s11 = sn1(1,:);
s12 = sn1(2,:);
s21 = s12;
s22 = s11;
z0 = 14;

% f = 0;

del_s = (1-s11).*(1-s22)-s12.*s21;

z11 = z0.*((1+s11).*(1-s22)+s12.*s21)./del_s;

z12 = z0.*(2*s12)./del_s;
z21 = z12;

z22 = z0.*((1-s11).*(1+s22)+s12.*s21)./del_s;

L = (z11-z12)./(j*2*pi*f);

C = 1./(j*2*pi*f.*z12);

figure(1);
hold on
plot(f/1e9,real(L)./1e-12);
xlim([0 500]);
hold off
xlabel('freq (GHz)');
ylabel('L (pH)');
ylim([0 0.15]);

figure(2);
hold on
plot(f/1e9,real(C)/1e-15);
xlim([0 500]);
hold off
xlabel('freq (GHz)');
ylabel('C (fF)');
ylim([0 2.5]);



