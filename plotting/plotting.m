clear

monitor = load("monitor.mat");
monitor = monitor.monitor_values;

delta_t = 192.58e-15;
delta_x = 100e-6;

Ez = monitor{1};
N = 1038;

num_monitors = length(monitor);

for i = 1:1:num_monitors
    Ez = monitor{i};
    voltage_temp = sum(Ez,2);
    voltage_temp = squeeze(voltage_temp)*6*delta_x;
    voltage(i,:) = voltage_temp(1:N);
end



t = (0:N-1)*delta_t;

figure(1);
hold on
temp_t = t/1e-9;
plot(temp_t,-voltage(1,:),'b');
plot(temp_t,-voltage(9,:),'r--');
plot(temp_t,-voltage(3,:),'k-.');

hold off
xlabel('Time (ns)')
ylabel('Voltage (V)')
legend('Port 1 [1 mm]','Port 2 [9 mm]','Reference [3 mm]')
grid on

N_fft = N*4;
F_k = fft(-voltage,N_fft,2);

% mag_phase_plot((0:N_fft-1)/N_fft/delta_t/1e9,F_k(index(1),:),'f (GHz)',2,'b');


x = (0:N_fft-1)/N_fft/delta_t/1e9;
    f_x_1 = F_k(1,:)/N_fft;
    f_x_9 = F_k(9,:)/N_fft;
    f_x_3 = F_k(3,:)/N_fft;

    figure(2);
    pl1 = subplot(2,1,1);
    hold on
    mag = abs(f_x_1);
    plot(x,mag','b');

    mag = abs(f_x_9);
    plot(x,mag','r--');
    
    mag = abs(f_x_3);
    plot(x,mag','k-.');

    hold off
    ylabel('Magnitude')
    grid on

%     yticks(0:0.2:1.2);

    pl2 = subplot(2,1,2);
    hold on
    pha = angle(f_x_1);
    plot(x,pha','b');
    pha = angle(f_x_9);
    plot(x,pha','r--');
    pha = angle(f_x_3);
    plot(x,pha','k-.');

    ylabel('Phase (rad)')
    hold off
    yticks([-pi -pi/2 0 pi/2 pi])
    yticklabels({'\pi', '-\pi/2', '0', '\pi/2', '\pi'})
    xlabel('Frequency (GHz)');

    grid on
    
    linkaxes([pl1, pl2],'x');




xlim([0,60]);


F_k = fft(voltage,N_fft,2);
    
    x = (0:N_fft-1)/N_fft/delta_t/1e9;
    f_x_1 = F_k(1,:)./F_k(3,:);
    f_x_2 = F_k(9,:)./F_k(3,:);

    figure(3);
    pl1 = subplot(2,1,1);
    hold on
    mag = abs(f_x_1);
    plot(x,mag','b*');

    mag = abs(f_x_2);
    plot(x,mag','r+');
    hold off
    ylabel('Magnitude')
    grid on

    yticks(0:0.2:1.2);
    legend('S_{11}','S_{21}');

    pl2 = subplot(2,1,2);
    hold on
    pha = angle(f_x_1);
    plot(x,pha','b*');
    pha = angle(f_x_2);

    plot(x,pha','r+');
    ylabel('Phase (rad)')
    hold off
    yticks([-pi -pi/2 0 pi/2 pi])
    yticklabels({'\pi', '-\pi/2', '0', '\pi/2', '\pi'})
    xlabel('Frequency (GHz)');

    grid on
    
    linkaxes([pl1, pl2],'x');
   

xlim([0,60]);
% subplot(2,1,1);
% ylim([0 0.1])
hold off


t = (0:N).*delta_t;



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