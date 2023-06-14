clear
close all

% load("colorblind_colormap.mat");

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

monitor = load("monitor.mat");
monitor_values = monitor.monitor_values;

delta_t = monitor.delta_t;
delta_x = monitor.delta_z;

Ez = monitor_values{1};
N = 5000;

num_monitors = length(monitor_values);

for i = 1:1:num_monitors
    Ez = monitor_values{i};
    num_cells = size(Ez,2);
    voltage_temp = sum(Ez,2);
    voltage_temp = squeeze(voltage_temp)*num_cells*delta_x;
    voltage(i,:) = -voltage_temp(1:N);
end



t = (0:N-1)*delta_t;
plot(t,voltage);

reference = gauspuls(t-35e-12,40e9,1).*2e-3;

N_fft = N*8;
% voltage(1,1:600) = 0;
F_k = fft(voltage,N_fft,2)/N_fft;
F_ref = fft(reference,N_fft)/N_fft;

% mag_phase_plot((0:N_fft-1)/N_fft/delta_t/1e9,F_k(index(1),:),'f (GHz)',2,'b');

x = (0:N_fft-1)/N_fft/delta_t/1e9;
f = (0:N_fft-1)/N_fft/delta_t;

G_f = (0.6+0.009*Z_0).*(((f/1e9)./(Z_0/(8*pi*(d*100)))).^2);
e_eff_f = epsilon_r-(epsilon_r-e_eff_0)./(1+G_f);

e_eff_f = 1;

f_x_ref = F_ref;
f_x_1 = F_k(1,:);
f_x_2 = F_k(7,:);

% time_plot(t/1e-9,[voltage(1:6,:)],1,{'r','g','b','k'})

time_plot(t/1e-12,[reference; voltage(1,:); voltage(7,:)],2,{'k','r--','b-.'})
legend('Ref','Port 1','Port 2');

xlabel('time (ps)')
ylabel('voltage (V)')

%phase correction
% f_x_ref = f_x_ref.*exp(j*2*pi*f.*sqrt(mu_0*epsilon_0*e_eff_f)*(-6.6e-3));
f_x_1 = f_x_1.*exp(-j*2*pi*f.*sqrt(mu_0*epsilon_0*e_eff_f)*(-2e-3));
f_x_2 = f_x_2.*exp(-j*2*pi*f.*sqrt(mu_0*epsilon_0*e_eff_f)*(-4e-3));

% mag_phase_plot(x,[f_x_ref;f_x_1;f_x_2;f_x_3],3,{'k','r','b--','m-.'})
% xlim([0 60])
% xlabel('freq (GHz)')

s11 = f_x_1./f_x_ref;
s21 = f_x_2./f_x_ref;

mag_phase_plot(x,[s11;s21],4,{'r--','b-.'})
xlim([0 60])
xlabel('freq (GHz)')
subplot(2,1,1)
legend('S11','S21')
yticks([0:0.25:1.25])

c = 1/sqrt(epsilon_0*mu_0);
%TM_0 surface wave mode critical freq
fT1 = (c/(2*pi*d))*sqrt(2/(epsilon_r-1))*atan(epsilon_r)/1e9;

%TE_1 surface wave
fT2 = c/(4*d*sqrt(epsilon_r-1))/1e9;

subplot(2,1,2);
xline(fT1,'k','LineWidth',1.5);
xline(fT2,'k','LineWidth',1.5);

 
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