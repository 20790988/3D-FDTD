clear

% fstart = 3e9;
% fstop = 5e9;
% fstep = 0.1e9;
% 
% fs = 100*fstop;
% 
% N = fs/fstep;
% ts = 1/fs;
% T = (N-1)*ts;
% 
% n = 0:N-1;
% 
% fn = cos(n*2*pi*2e9/fs);
% 
% % time_plot(n*ts/1e-9,fn,'t (ns)',1);
% 
% % F_k = load('tempmat.mat').F_k;
% F_k = zeros(1,N);
% 
% % F_k(fstart/fs*N+1:fstop/fs*N+1) = 1;
% 
% F_k = normpdf(n/N*fs,(fstop-fstart)/2+fstart,(fstop-fstart)/2);
% F_k = F_k/max(F_k);
% 
% %mirroring step
% F_k(N:-1:N/2+2) = F_k(2:N/2);
% 
% mag_phase_plot(n/N*fs/1e9,F_k,'freq GHz',2);
% 
% time_plot(n*ts/1e-9,ifft(F_k),'t (ns)',3,true);
% 
% % mag_phase_plot(n/N*fs/1e9,fft(ifft(F_k)),'freq GHz',4);


delta_t = 3.8516e-12;
T_max = 519*delta_t;
t = 0:delta_t:(T_max-delta_t);
t = t-0.4e-9;
t0 = 350*delta_t;
T = 40*delta_t;

% f_t = exp(-(t-t0).^2./(T^2));
f_t = gauspuls(t,4e9,1);

time_plot(t/1e-9,f_t,'t (ns)',1);

mag_phase_plot((0:518)/518/delta_t/1e9,fft(f_t),'f (GHz)',3)

function mag_phase_plot(x,f_x,x_text,fig_no)
    figure(fig_no);
    pl1 = subplot(2,1,1);
    mag = abs(f_x);
    stem(x,mag);
    ylabel('Magnitude')
    grid on

    pl2 = subplot(2,1,2);
    pha = angle(f_x);
    pha(mag<1e-8) = 0;
    stem(x,pha);
    ylabel('Phase (rad)')
    
    temp_lim = ylim;
    if(max(abs(temp_lim))<pi)
        ylim([-pi pi])
    end
    xlabel(x_text);
    grid on
    
    linkaxes([pl1, pl2],'x');
end

function time_plot(x,f_x,x_text,fig_no,flip)
    figure(fig_no);
    N = length(x);
    
    if exist('flip','var') && flip
        f_x = [f_x(N/2+1:N) f_x(1:N/2)];
        x = [-x(N/2+1:-1:2) x(1:N/2)];
    end
    
    plot(x,f_x);
    ylabel('Amplitude')
    xlabel(x_text);
    grid on
end
