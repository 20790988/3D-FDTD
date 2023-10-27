clear
close all

%====================SETTINGS=====================%

result_filename = "monitor_14.mat";

%Line properties
    %distance between plates and width of line in meters
%     d = 0.3e-6;
%     w = 4e-6;

    epsilon_r = 4.6;

    is_microstrip = false;

%Sampling

    % max number of samples used
    N_max = Inf;

    %upsampling factor for fft. N*k samples used
    k_fft = 1000;
    N_fft_max = 8000;

    %how many samples of port 1 voltage are set to zero
    N_zero = 0;

%Indices
    port_index = [1,29];

    number_of_ports = length(port_index);

    reference_index = 0;
        %for index 0 the source values are used

    %phase correction distance measured away from scattering object
    % in meters
    phase_distance_ref = -10.8e-6;
    phase_distance = [-11.8e-6 -11.8e-6];


% Theoretical values
%     line_length = 30e-6;

%=========================================%

load('e_eff.mat');
load_eff_f = e_eff_from_sim;
load_f = f;

load('sparams_equiv_circuit.mat');

s11_eq_cir = s11;
s21_eq_cir = s21;
f_eq_cir = f;


epsilon_0 = 8.8542e-12;
mu_0 = 1.2566e-6;

if is_microstrip
    e_eff_0 = (epsilon_r+1)/2+(epsilon_r-1)/(2*sqrt(1+12*d/w));

    if w/d <= 1
        Z_0 = (60/sqrt(e_eff_0))*log((8*d)/(w)+(w)/(4*d));
    else
        Z_0 = (120*pi)/(sqrt(e_eff_0)*(w/d+1.393+0.667*log(w/d+1.444)));
    end
end

monitor = load(result_filename);
monitor_values = monitor.monitor_values;

delta_t = monitor.delta_t;
delta_x = monitor.delta_z;

Ez = cell2mat(monitor_values{1}(3));

N = length(Ez);

if N_max < N
    N = N_max;
end

num_monitors = length(monitor_values);

for i = 1:1:num_monitors
    Ez = cell2mat(monitor_values{i}(3));
    num_cells = size(Ez,2);
    voltage_temp = sum(Ez,2);
    voltage_temp = squeeze(voltage_temp)*num_cells*delta_x;
    voltage(i,:) = voltage_temp(1:N);
end


t = (0:N-1)*delta_t;
figure(1)
plot(t./1e-12,voltage);
grid on
xlabel('Time (ps)')
ylabel('Voltage')


if reference_index == 0 
    reference = monitor.source_from_bootstrap;
    reference = permute(reference,[2 3 1]);
    middle = floor(size(reference,1)/2);
    reference = reference(middle,:,:);
    num_cells = size(reference,2);
    voltage_temp = sum(reference,2);
    voltage_temp = squeeze(voltage_temp)*num_cells*delta_x;
    voltage_temp(N) = 0;

    reference = -voltage_temp(1:N)';
else
    reference = voltage(reference_index,:);
end

N_fft = N*k_fft;

% reference(N_zero:end) = 0;
% voltage(port_1_index,1:N_zero) = 0;


F_k = fft(voltage,N_fft,2)/N_fft;
F_ref = fft(reference,N_fft)/N_fft;

x = (0:N_fft-1)/N_fft/delta_t/1e9;
f = (0:N_fft-1)/N_fft/delta_t;

if N_fft_max< length(F_k)
    F_k = F_k(:,1:N_fft_max);
    F_ref = F_ref(:,1:N_fft_max);

    x = x(1:N_fft_max);
    f = f(1:N_fft_max);
end



if is_microstrip
    G_f = (0.6+0.009*Z_0).*(((f/1e9)./(Z_0/(8*pi*(d*100)))).^2);
    e_eff_f = epsilon_r-(epsilon_r-e_eff_0)./(1+G_f);
else
    e_eff_f = epsilon_r;
end

f_x_ref = F_ref;
f_x = F_k(port_index,:);

time_plot(t/1e-12,[reference; ...
    voltage(port_index(1),:); ...
    voltage(port_index(2),:)],2,{'k','r--','b-.'})

legend('Ref','Port 1','Port 2');

xlabel('time (ps)')
ylabel('voltage (V)')

e_eff_f_from_file = interp1(load_f,load_eff_f,f,'linear');

%phase correction
f_x_ref = f_x_ref.*exp(j*2*pi*f.*sqrt(mu_0*epsilon_0*e_eff_f_from_file)*(phase_distance_ref));
f_x = f_x.*exp(-j*2*pi*f.*sqrt(mu_0*epsilon_0*e_eff_f_from_file).*(phase_distance'));

sn1 = f_x./f_x_ref;


%theoretical PEC

% s11_th = zeros(1,length(sn1));
% beta = 2*pi*f.*sqrt(mu_0*epsilon_0.*e_eff_f);
% s21_th = exp(-1i*beta*line_length);

%inductex mag

% freq_inductex=[1E10,1.325E11,2.55E11,3.775E11,5E11];
% x_inductex = freq_inductex./1e9;
% 
% s21_inductex_mag=[-7.33162E-6,-0.00127104,-0.00455141,-0.00943939,-0.0153382];
% s11_inductex_mag=[-57.7176,-35.3275,-29.7882,-26.6209,-24.5127];

%toepfer

% S11 = readmatrix('S11_toepfer_extract.csv');
% S21 = readmatrix('S21_toepfer_extract.csv');

figure(9)
hold on


plot(x,20*log10(abs(sn1(1,:))),'r','LineWidth',2);
plot(x,20*log10(abs(sn1(2,:))),'r--','LineWidth',2);
% 
% plot(x_inductex,s11_inductex_mag,'k');
% plot(x_inductex,s21_inductex_mag,'k--');

% plot(S11(:,1),S11(:,2),'m');
% plot(S21(:,1),S21(:,2),'m--');

plot(f_eq_cir/1e9,20*log10(abs(s11_eq_cir)),'b');
plot(f_eq_cir/1e9,20*log10(abs(s21_eq_cir)),'b--');

hold off
ylabel('Magnitude')
grid on
xlim([0 500]);
xlabel('freq (GHz)')
ylabel('Magnitude (dB)')
legend('FDTD s11','FDTD s21','Eqv cir s11','Eqv cir s21');

ylim([-70 10])


% figure(5);
% hold on
% plot(x,20*log10(abs(sn1(2,:))),'b--',LineWidth=2);
% 
% % plot(x,20*log10(abs(s21_th)),'b-.',LineWidth=2)
% 
% plot(x_inductex,s21_inductex_mag,'r-.',LineWidth=2)
% 
% 
% plot(S21(:,1),S21(:,2),'k:',LineWidth=2);
% hold off
% ylabel('Magnitude')
% grid on
% xlim([0 500]);
% xlabel('freq (GHz)')
% ylabel('Magnitude (dB)')
% legend('FDTD s21','Inductex s21','Toepfer S21');
% ylim([-1 1]);
% 
% figure(6);
% hold on
% plot(x,20*log10(abs(sn1(1,:))),'b--',LineWidth=2)
% plot(x_inductex,s11_inductex_mag,'r-.',LineWidth=2)
% plot(S11(:,1),S11(:,2),'k:',LineWidth=2);
% 
% hold off
% ylabel('Magnitude')
% grid on
% xlim([0 500]);
% xlabel('freq (GHz)')
% ylabel('Magnitude (dB)')
% legend('FDTD s11','Inductex s11','Toepfer S11');
% ylim([-70 0])

% figure(7);
% hold on
% plot(x,-2*pi.*f.*(4)*1e-6./angle(sn1(2,:)),'b--',LineWidth=2);
% hold off
% ylabel('Phase Velocity (m/s)')
% grid on
% xlim([0 500]);
% xlabel('freq (GHz)')
% legend('FDTD');
% ylim([0 3e8])


figure(8);
hold on
plot(x,angle(sn1(1,:)).*180/pi,'r');
plot(x,angle(sn1(2,:)).*180/pi,'r--');

plot(f_eq_cir/1e9,angle(s11_eq_cir).*180/pi,'b');
plot(f_eq_cir/1e9,angle(s21_eq_cir).*180/pi,'b--');
hold off
ylabel('Phase angle (degrees)')
grid on
xlim([0 500]);
xlabel('freq (GHz)')
legend('FDTD s11','FDTD s21','Eqv cir s11','Eqv cir s21');


save('sparams.mat','f','sn1');


% c = 1/sqrt(epsilon_0*mu_0);
% %TM_0 surface wave mode critical freq
% fT1 = (c/(2*pi*d))*sqrt(2/(epsilon_r-1))*atan(epsilon_r)/1e9;
% 
% %TE_1 surface wave
% fT2 = c/(4*d*sqrt(epsilon_r-1))/1e9;


 
function mag_phase_plot(x,f_x,fig_no,linestyles,dB_flag)
    
    figure(fig_no);

    s = size(f_x);
    s = s(1);
    linestyles{s+1} = 0;
    
    mag = abs(f_x);
    if exist('dB_flag','var') && dB_flag
        mag = 20*log10(mag);
    end
    
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
    
%     yticks([-180 -90 0 90 180])

    grid on
    
    linkaxes([pl1, pl2],'x');
end

function real_imag_plot(x,f_x,fig_no,linestyles)
    
    figure(fig_no);

    s = size(f_x);
    s = s(1);
    linestyles{s+1} = 0;

    re = real(f_x);
    pl1 = subplot(2,1,1);

    hold on
    for i = 1:s
        if ~isempty(linestyles{i})
            plot(x,re(i,:),linestyles{i});
        else
            plot(x,re(i,:));
        end
    end
    hold off
    
    ylabel('Real')
    grid on

    pl2 = subplot(2,1,2);
    im = imag(f_x);
    
    hold on
    for i = 1:s
        if ~isempty(linestyles{i})
            plot(x,im(i,:),linestyles{i});
        else
            plot(x,im(i,:));
        end
    end
    hold off
    ylabel('Imag')
    
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



function pulse = gaus_derv(t,mu,sigma)
    pulse = (t-mu).*exp(-((t-mu).^2)/(2*sigma^2));
    pulse =  pulse./max(pulse,[],'all');
end