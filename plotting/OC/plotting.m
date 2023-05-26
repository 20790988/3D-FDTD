clear
close all

epsilon_0 = 8.8542e-12;
mu_0 = 1.2566e-6;
e_eff = 6.5;
%  e_eff = 7;

monitor = load("monitor.mat");
monitor_values = monitor.monitor_values;

delta_t = monitor.delta_t;
delta_x = monitor.delta_z;

Ez = monitor_values{1};
N = 2200;

num_monitors = length(monitor_values);

for i = 1:1:num_monitors
    Ez = monitor_values{i};
    voltage_temp = sum(Ez,2);
    voltage_temp = squeeze(voltage_temp)*6*delta_x;
    voltage(i,:) = -voltage_temp(1:N);
end

t = (0:N-1)*delta_t;
% plot(t,voltage);

reference = gauspuls(t-35e-12,40e9,1).*2.75;


N_fft = N*8;
voltage(1,1:945) = 0;
F_k = fft(voltage,N_fft,2)/N_fft;
F_ref = fft(reference,N_fft)/N_fft;

% mag_phase_plot((0:N_fft-1)/N_fft/delta_t/1e9,F_k(index(1),:),'f (GHz)',2,'b');

x = (0:N_fft-1)/N_fft/delta_t/1e9;
f = (0:N_fft-1)/N_fft/delta_t;

f_x_ref = F_ref;
f_x_1 = F_k(1,:);
f_x_4 = F_k(4,:);

% time_plot(t/1e-9,[voltage(4,:)],1,{'r','g','b','k'})

time_plot(t/1e-12,[reference; voltage(1,:); voltage(4,:)],2,{'k','r','b--'})
legend('x=3.00 mm (Ref)','x=0.96 mm (Port 1)','x=3.84 mm');

xlabel('time (ps)')
ylabel('voltage (V)')

%phase correction
f_x_ref = f_x_ref.*exp(j*2*pi*f*sqrt(mu_0*epsilon_0*e_eff)*(-6.18e-3));
f_x_1 = f_x_1.*exp(-j*2*pi*f*sqrt(mu_0*epsilon_0*e_eff)*(-8.22e-3));

% mag_phase_plot(x,[f_x_ref;f_x_1;f_x_2;f_x_3],3,{'k','r','b--','m-.'})
% xlim([0 60])
% xlabel('freq (GHz)')

s11 = f_x_1./f_x_ref;

mag_phase_plot(x,[s11],4,{'r','b--','m-.'})
xlim([0 60])
xlabel('freq (GHz)')
subplot(2,1,1)
legend('S11')
yticks([0:0.25:1.25])

  



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