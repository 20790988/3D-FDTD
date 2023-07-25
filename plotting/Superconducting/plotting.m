clear
close all

%====================SETTINGS=====================%

result_filename = "monitor_08_fix_3_long.mat";

%Line properties
    %distance between plates and width of line in meters
    d = 0.3e-6;
    w = 4e-6;

    epsilon_r = 4;

    is_microstrip = false;

%Sampling

    % max number of samples used
    N_max = Inf;

    %upsampling factor for fft. N*k samples used
    k_fft = 1000;
    N_fft_max = 1000;

    %how many samples of port 1 voltage are set to zero
    N_zero = 0;

%Indices
    port_1_index = 1;
    port_2_index = 13;

    reference_index = 2;
        %for index 0 the source values are used

    %phase correction distance measured away from scattering object
    % in meters
    phase_distance_ref = 0;
    phase_distance_1 = 0;
    phase_distance_2 = 0;


% Theoretical values
    line_length = 11e-6;

%=========================================%

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
    reference = monitor.source_val_E(1:N).*0.0011;
else
    reference = voltage(reference_index,:);
end

N_fft = N*k_fft;

voltage(1,1:N_zero) = 0;

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
f_x_1 = F_k(port_1_index,:);
f_x_2 = F_k(port_2_index,:);

time_plot(t/1e-12,[reference; ...
    voltage(port_1_index,:); ...
    voltage(port_2_index,:)],2,{'k','r--','b-.'})

legend('Ref','Port 1','Port 2');

xlabel('time (ps)')
ylabel('voltage (V)')

%phase correction
f_x_ref = f_x_ref.*exp(j*2*pi*f.*sqrt(mu_0*epsilon_0*e_eff_f)*(phase_distance_ref));
f_x_1 = f_x_1.*exp(-j*2*pi*f.*sqrt(mu_0*epsilon_0*e_eff_f)*(phase_distance_1));
f_x_2 = f_x_2.*exp(-j*2*pi*f.*sqrt(mu_0*epsilon_0*e_eff_f)*(phase_distance_2));

s11 = f_x_1./f_x_ref;
s21 = f_x_2./f_x_ref;

%theoretical
s11_th = zeros(1,length(s11));
beta = 2*pi*f.*sqrt(mu_0*epsilon_0.*e_eff_f);
s21_th = exp(-1i*beta*line_length);

mag_phase_plot(x,[s11;s21;s11_th;s21_th],4,{'r--','b-.','k--','k-.'});
xlim([0 100]);
xlabel('freq (GHz)')
subplot(2,1,1)
legend('S11','S21')
yticks([0:0.25:1.25])
ylim([0,1.25]);

real_imag_plot(x,[s11;s21],5,{'r--','b-.','k--','k-.'});
xlim([0 100]);
xlabel('freq (GHz)')
subplot(2,1,2)
ylim([-1.25,1.25]);

subplot(2,1,1)
ylim([-1.25,1.25]);
legend('S11','S21')

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