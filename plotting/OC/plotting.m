clear
close all

epsilon_0 = 8.8542e-12;
mu_0 = 1.2566e-6;
e_eff = 9.6;
%  e_eff = 7;

monitor = load("monitor_bwl.mat");
monitor = monitor.monitor_values;

delta_t = 115.55e-15;
delta_x = 60e-6;

Ez = monitor{1};
N = 2700;

num_monitors = length(monitor);

for i = 1:1:num_monitors
    Ez = monitor{i};
    voltage_temp = sum(Ez,2);
    voltage_temp = squeeze(voltage_temp)*6*delta_x;
    voltage(i,:) = voltage_temp(1:N);
end

t = (0:N-1)*delta_t;
% plot(t,voltage);

voltage(3,840:end) = 0;

reference = gauspuls(t-35e-12,40e9,1).*2.7;

figure(1);
hold on
temp_t = t/1e-9;
plot(temp_t,-voltage(1,:),'b');
plot(temp_t,-voltage(3,:),'r--');
plot(temp_t,reference,'k-.');

hold off
xlabel('Time (ns)')
ylabel('Voltage (V)')
legend('x = -1.02 mm','x = 0.9 mm','x = 0 mm (reference)')
grid on

N_fft = N*8;
F_k = fft(-voltage,N_fft,2);

% mag_phase_plot((0:N_fft-1)/N_fft/delta_t/1e9,F_k(index(1),:),'f (GHz)',2,'b');


x = (0:N_fft-1)/N_fft/delta_t/1e9;
f = (0:N_fft-1)/N_fft/delta_t;

    f_x_1 = F_k(1,:)/N_fft;
    f_x_4 = F_k(3,:)/N_fft;

    %phase correction
    f_x_4 = f_x_4.*exp(j*2*pi*f*sqrt(mu_0*epsilon_0*e_eff)*(-3.36e-3));
%     f_x_4(3,840:end) = 0;
    f_x_1 = f_x_1.*exp(-j*2*pi*f*sqrt(mu_0*epsilon_0*e_eff)*(-6.3e-3));

    f_x_ref = fft(reference,N_fft)/N_fft;
    f_x_ref = f_x_ref.*exp(j*2*pi*f*sqrt(mu_0*epsilon_0*e_eff)*(-7.2e-3));


    figure(2);
    pl1 = subplot(2,1,1);
    hold on
    mag = abs(f_x_1);
    plot(x,mag','b');

    mag = abs(f_x_4);
    plot(x,mag','r--');

    mag = abs(f_x_ref);
    plot(x,mag','k-.');
   
    hold off
    ylabel('Magnitude')
    grid on

%     yticks(0:0.2:1.2);

    pl2 = subplot(2,1,2);
    hold on
    pha = angle(f_x_1)*180/pi;
    plot(x,pha','b');
    pha = angle(f_x_4)*180/pi;
    plot(x,pha','r--');
    pha = angle(f_x_ref)*180/pi;
    plot(x,pha','k-.');

    ylabel('Phase (deg)')
    hold off
%     yticks([-pi -pi/2 0 pi/2 pi])
%     yticklabels({'\pi', '-\pi/2', '0', '\pi/2', '\pi'})

    yticks([-180 -90 0 90 180])
%     yticklabels({'\pi', '-\pi/2', '0', '\pi/2', '\pi'})
    xlabel('Frequency (GHz)');

    grid on
    
    linkaxes([pl1, pl2],'x');

    xlim([0,60]);

    s11 = f_x_1./f_x_ref;
%     s11 = s11.*exp(-2*j*2*pi*f*sqrt(mu_0*epsilon_0*e_eff)*(-5.22e-3));

    figure(3);
    pl1 = subplot(2,1,1);
    hold on
    mag = real(s11);
    plot(x,mag','b');
    mag = imag(s11);
    plot(x,mag','r');

    hold off
    ylabel('Magnitude')
    grid on

%     ylim([0, 1.2])
    legend('S_{11}');
    yticks(0:0.2:1.2)

    pl2 = subplot(2,1,2);
    hold on
%     pha = angle(s11)*180/pi;
    plot(x,pha','b');
    
    ylabel('Phase (deg)')
    hold off
    
%     yticklabels({'\pi', '-\pi/2', '0', '\pi/2', '\pi'})
    yticks([-180 -90 0 90 180])
    xlabel('Frequency (GHz)');

    grid on
    
    linkaxes([pl1, pl2],'x');
   

xlim([0,60]);
% subplot(2,1,1);
% ylim([0 0.1])
hold off




function mag_phase_plot(x,f_x,x_text,fig_no)
    figure(fig_no);
    pl1 = subplot(2,1,1);
    mag = abs(f_x);
    plot(x,mag');
    ylabel('Magnitude')
    grid on

    pl2 = subplot(2,1,2);
    pha = angle(f_x);
    
    mag_min = 0.001*max(mag,[],"all");

%     pha(:) = 0;
    plot(x,pha');
    ylabel('Phase (rad)')
    
    temp_lim = ylim;
    if(max(abs(temp_lim))<pi)
        ylim([-pi pi])
    end
    yticks([-pi -pi/2 0 pi/2 pi])
    yticklabels({'\pi', '-\pi/2', '0', '\pi/2', '\pi'})
    xlabel(x_text);

    grid on
    
    linkaxes([pl1, pl2],'x');
end