clear

fstart = 2e9;
fstop = 5e9;
fstep = 0.1e9;

fs = 4*fstop;

N = fs/fstep;
ts = 1/fs;
T = (N-1)*ts;

n = 0:N-1;

fn = cos(n*2*pi*2e9/fs);

% time_plot(n*ts/1e-9,fn,'t (ns)',1);

% F_k = load('tempmat.mat').F_k;
F_k = zeros(1,N);

F_k(fstart/fs*N+1:fstop/fs*N+1) = 1;

% F_k=normpdf(n/N*fs,(fstop-fstart)/2+fstart,(fstop-fstart)/2);

%mirroring step
F_k(N:-1:N/2+2) = F_k(2:N/2);

mag_phase_plot(n/N*fs/1e9,F_k,'freq GHz',2);

time_plot(n*ts/1e-9,ifft(F_k),'t (ns)',3,true);

mag_phase_plot(n/N*fs/1e9,fft(ifft(F_k)),'freq GHz',4);

function mag_phase_plot(x,f_x,x_text,fig_no)
figure(fig_no);
pl1 = subplot(2,1,1);
mag = abs(f_x);
plot(x,mag);
ylabel('Magnitude')

pl2 = subplot(2,1,2);
pha = angle(f_x);
pha(mag<1e-8) = 0;
plot(x,pha);
ylabel('Phase (rad)')

temp_lim = ylim;
if(max(abs(temp_lim))<pi)
    ylim([-pi pi])
end
xlabel(x_text);

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
end
