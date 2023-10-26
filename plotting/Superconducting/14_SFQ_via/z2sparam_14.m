clear
close all

f = 1e9:1e9:500e9;

L = 0.355e-12;
C = 11.77e-15;

z11 = j*pi*2*f*L+1./(j*2*pi*f*C);
z12 = 1./(j*2*pi*f*C);
z21 = z12;
z22 = z11;
z0 = 5;

del_z = (z11 +z0).*(z22+ z0) -z12.*z21;

s11 = ((z11-z0).*(z22+z0)-z12.*z21)./del_z;
s12 = (2*z12.*z0)./del_z;
s21 = s12;
s22 = ((z11+z0).*(z22-z0)-z12.*z21)./del_z;

mag_phase_plot(f/1e9,s11,1,{'k'},true);
mag_phase_plot(f/1e9,s21,2,{'k'},true);

save('sparams_equiv_circuit.mat','s11','s21','f');

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



