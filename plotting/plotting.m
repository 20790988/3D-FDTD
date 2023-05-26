clear
close all

load("colorblind_colormap.mat");

epsilon_0 = 8.8542e-12;
mu_0 = 1.2566e-6;
epsilon_r = 9.6;

d = 0.6e-3;
w = 0.6e-3;

e_eff_0 = (epsilon_r+1)/2+(epsilon_r-1)/(2*sqrt(1+12*d/w));

Z_0 = 0;

if w/d <= 1
    Z_0 = (60/sqrt(e_eff_0))*log((8*d)/(w)+(w)/(4*d));
else
    Z_0 = (120*pi)/(sqrt(e_eff_0)*(w/d+1.393+0.667*log(w/d+1.444)));
end

e_eff = 7;

monitor = load("monitor_bigger.mat");
monitor_values = monitor.monitor_values;

delta_t = monitor.delta_t;
delta_x = monitor.delta_z;

Ez = monitor_values{1};
N = 1600;

num_monitors = length(monitor_values);

for i = 1:1:num_monitors
    Ez = monitor_values{i};
    voltage_temp = sum(Ez,2);
    voltage_temp = squeeze(voltage_temp)*6*delta_x;
    voltage(i,:) = -voltage_temp(1:N);
end



t = (0:N-1)*delta_t;
% plot(t,voltage);

reference = gauspuls(t-35e-12,40e9,1).*1.65;




N_fft = N*8;
voltage(1,1:600) = 0;
F_k = fft(voltage,N_fft,2)/N_fft;
F_ref = fft(reference,N_fft)/N_fft;

% mag_phase_plot((0:N_fft-1)/N_fft/delta_t/1e9,F_k(index(1),:),'f (GHz)',2,'b');

x = (0:N_fft-1)/N_fft/delta_t/1e9;
f = (0:N_fft-1)/N_fft/delta_t;

G_f = (0.6+0.009*Z_0).*(((f/1e9)./(Z_0/(8*pi*(d*100)))).^2);
e_eff_f = epsilon_r-(epsilon_r-e_eff_0)./(1+G_f);

f_x_ref = F_ref;
f_x_1 = F_k(1,:);
f_x_2 = F_k(7,:);
f_x_3 = F_k(13,:);

% time_plot(t/1e-9,[voltage(1:6,:)],1,{'r','g','b','k'})

time_plot(t/1e-12,[reference; voltage(1,:); voltage(15,:); voltage(13,:)],2,{'k','r','b--','m-.'})
legend('Ref','Port 1','Port 2','Port 3');

xlabel('time (ps)')
ylabel('voltage (V)')

%phase correction
f_x_ref = f_x_ref.*exp(j*2*pi*f.*sqrt(mu_0*epsilon_0*e_eff_f)*(-6.6e-3));
f_x_1 = f_x_1.*exp(-j*2*pi*f.*sqrt(mu_0*epsilon_0*e_eff_f)*(-8.6e-3));
f_x_2 = f_x_2.*exp(-j*2*pi*f.*sqrt(mu_0*epsilon_0*e_eff_f)*(-8.9e-3));
f_x_3 = f_x_3.*exp(-j*2*pi*f.*sqrt(mu_0*epsilon_0*e_eff_f)*(-8.9e-3));

% mag_phase_plot(x,[f_x_ref;f_x_1;f_x_2;f_x_3],3,{'k','r','b--','m-.'})
% xlim([0 60])
% xlabel('freq (GHz)')

s11 = f_x_1./f_x_ref;
s21 = f_x_2./f_x_ref;
s31 = f_x_3./f_x_ref;

mag_phase_plot(x,[s11;s21;s31],4,{'r','b--','m-.'})
xlim([0 60])
xlabel('freq (GHz)')
subplot(2,1,1)
legend('S11','S21','S31')
yticks([0:0.25:1.25])

c = 1/sqrt(epsilon_0*mu_0);
%TM_0 surface wave mode critical freq
fT1 = (c/(2*pi*d))*sqrt(2/(epsilon_r-1))*atan(epsilon_r)/1e9

%TE_1 surface wave
fT2 = c/(4*d*sqrt(epsilon_r-1))/1e9

subplot(2,1,2);
xline(fT1,'k','LineWidth',2);
xline(fT2,'k','LineWidth',2);

 
function mag_phase_plot(x,f_x,fig_no,linestyles)
    
    figure(fig_no);

    s = size(f_x);
    s = s(1);
    linestyles{s+1} = 0;

    mag = abs(f_x);
    pl1 = subplot(2,1,1);

    hold on
    for i = 1:s
        if ~isempty(linestyles{i})
            plot(x,mag(i,:),linestyles{i});
        else
            plot(x,mag(i,:));
        end
    end
    hold off
    
    ylabel('Magnitude')
    grid on

    pl2 = subplot(2,1,2);
    pha = angle(f_x)*180/pi;
    
    hold on
    for i = 1:s
        if ~isempty(linestyles{i})
            plot(x,pha(i,:),linestyles{i});
        else
            plot(x,pha(i,:));
        end
    end
    hold off
    ylabel('Phase (deg)')
    
    yticks([-180 -90 0 90 180])

    grid on
    
    linkaxes([pl1, pl2],'x');
end

function time_plot(x,f_x,fig_no,linestyles)
    
    figure(fig_no);
    
    s = size(f_x);
    s = s(1);
    linestyles{s+1} = 0;
    hold on
    for i = 1:s
        if ~isempty(linestyles{i})
            plot(x,f_x(i,:),linestyles{i});
        else
            plot(x,f_x(i,:));
        end
    end
    hold off
    grid on
end