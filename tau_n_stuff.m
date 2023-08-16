clear
close all

j = 1i;

lambda_0 = 85e-9;
T_c = 9.25;
T_op = 4.2;

tau_n = 61e-15;
mu_0 = 1.2566e-6;

cond_DC = 6.7e6;

f = 1e9:1e9:500e9;
w = 2*pi*f;

lambda = lambda_0/(sqrt(1-(T_op/T_c)^4));

cond_n_0 = cond_DC*(T_op/T_c).^4;

cond_n = cond_DC*(T_op/T_c).^4./(j*w*tau_n+1);

cond_s = 1./(lambda.^2*mu_0*j*w);


real_imag_plot(f./1e9,cond_n,1,{'r'});
real_imag_plot(f./1e9,cond_s,2,{'b'});

% mag_phase_plot(f./1e9,cond_n,1,{'r'});
% mag_phase_plot(f./1e9,cond_s,2,{'b'});


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