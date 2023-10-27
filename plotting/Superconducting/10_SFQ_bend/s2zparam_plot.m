clear
close all
load("LC_10.mat");
L_(:,1) = L;
C_(:,1) = C;

load("LC_11.mat");
L_(:,2) = L;
C_(:,2) = C;

load("LC_12.mat");
L_(:,3) = L;
C_(:,3) = C;

load("LC_13.mat");
L_(:,4) = L;
C_(:,4) = C;

% newcolors = [0 0 255
%              255 0 0
%              255 0 255
%              0 255 255]./255;
         
colororder(acc_colormap);

figure(1);
hold on
plot(f/1e9,real(L_)./1e-12,'LineWidth',2);
xlim([0 500]);
hold off
xlabel('freq (GHz)');
ylabel('L (pH)');
% yline(0.0383316e-9/1e-12,'b--','LineWidth',2);  
% yline(0.00932666e-9/1e-12,'r--','LineWidth',2);
grid on
legend('Manhattan','ch = 2.55','ch = 4.35','ch = 4.95','Location','southeast')
ylim([0 0.35]);


figure(2);
colororder(acc_colormap);
hold on
plot(f/1e9,real(C_)/1e-15,'LineWidth',2);
xlim([0 500]);
hold off
xlabel('freq (GHz)');
ylabel('C (fF)');
% yline(34.0121e-15./1e-15,'b--','LineWidth',2)
% yline(33.4045e-15./1e-15,'r--','LineWidth',2)
ylim([0 3.5]);
grid on
% legend('FDTD','Formula from Javadzadeh et al','Location','southeast')

function map = acc_colormap()
    map = [204 121 167
    213 94 0;
    0 114 178;
    240 228 66;
    0 158 115;
    86 180 233;
    230 159 0;
    0 0 0;]./255;
end


