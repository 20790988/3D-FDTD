clear
close all

load("colorblind_colormap.mat");

epsilon_0 = 8.8542e-12;
mu_0 = 1.2566e-6;
e_eff = 6.5;

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
    voltage(i,:) = voltage_temp(1:N);
end



t = (0:N-1)*delta_t;
% plot(t,voltage);

reference = gauspuls(t-35e-12,40e9,1).*1.65;




N_fft = N*8;
F_k = fft(voltage,N_fft,2)/N_fft;
F_ref = fft(reference,N_fft)/N_fft;

% mag_phase_plot((0:N_fft-1)/N_fft/delta_t/1e9,F_k(index(1),:),'f (GHz)',2,'b');

x = (0:N_fft-1)/N_fft/delta_t/1e9;
f = (0:N_fft-1)/N_fft/delta_t;

voltage(1,1:600) = 0;

f_x_ref = F_ref;
f_x_1 = F_k(1,:);
f_x_2 = F_k(7,:);
f_x_3 = F_k(13,:);

time_plot(t/1e-9,[reference; voltage(1,:); voltage(15,:); voltage(13,:)],2,{'k','r','b','m'})

xlabel('time (ns)')
ylabel('voltage (V)')


s11 = f_x_1./f_x_ref;
s21 = f_x_2./f_x_ref;
s31 = f_x_3./f_x_ref;

mag_phase_plot(x,[s11;s21;s31],3,{'r','b','m'})
xlim([0 60])
subplot(2,1,1)
legend('S11','S21','S31')

%     %phase correction
%     f_x_4 = f_x_4.*exp(j*2*pi*f*sqrt(mu_0*epsilon_0*e_eff)*(-3.36e-3));
% %     f_x_4(3,840:end) = 0;
%     f_x_1 = f_x_1.*exp(-j*2*pi*f*sqrt(mu_0*epsilon_0*e_eff)*(-6.3e-3));
% 
%     f_x_ref = fft(reference,N_fft)/N_fft;
%     f_x_ref = f_x_ref.*exp(j*2*pi*f*sqrt(mu_0*epsilon_0*e_eff)*(-7.2e-3));
% 
%     s11 = f_x_1./f_x_ref;
% %     s11 = s11.*exp(-2*j*2*pi*f*sqrt(mu_0*epsilon_0*e_eff)*(-5.22e-3));
 
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