clear
close all

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

f = 0:1e9:1e15;

t_n = 4.3e-15*(4.2/9.3)^(-4);

scale_1 = 1./((j*2*pi*f*t_n)+1);

t_n = 61e-15;
scale_2 = 1./((j*2*pi*f*t_n)+1);

mag_phase_plot(f,[scale_1;scale_2],1,{'r','b'});


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