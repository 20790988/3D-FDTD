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




N_fft = N*4;


% mag_phase_plot((0:N_fft-1)/N_fft/delta_t/1e9,F_k(index(1),:),'f (GHz)',2,'b');





t = (0:N_fft-1).*delta_t;
t_temp = t-60e-12;

source_signal = gauspuls(t_temp,50e9,1);
x = (0:N_fft-1)/N_fft/delta_t/1e9;

figure(1)
plot(t,source_signal);

mag_phase_plot(x,fft(source_signal)/N_fft,'',2);


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