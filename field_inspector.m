clear -except monitor
close all

%====================SETTINGS=====================%
db_axis = true;
%=========================================%
result_filename = "field_cap.mat";

%Line properties
    %distance between plates and width of line in meters
    

monitor = load(result_filename);
monitor_values = monitor.bootstrap_values{1};

delta_t = monitor.delta_t;
delta_x = monitor.delta_z;

Ex = monitor_values{1};
Ey = monitor_values{2};
Ez = monitor_values{3};

Hx = monitor_values{4};
Hy = monitor_values{5};
Hz = monitor_values{6};

N = length(Ez);

num_monitors = length(monitor_values);

figure(2)
s = size(Ez);
plot(squeeze(Ez(floor(s(1)/2),10,:)));

figure(1)

Hz = sqrt(Hx.^2+Hy.^2+Hz.^2);
Ez = sqrt(Ex.^2+Ey.^2+Ez.^2);

E_max = max(abs(Ez),[],"all");
H_max = max(abs(Hz),[],"all");


Ez = permute(Ez,[2 1 3]);
Hz = permute(Hz,[2 1 3]);
for ii = 1:N
    if sum(Ez(:,:,ii),"all") == 0
        continue;
    end
    figure(1)
    colormap jet
    surf(Ez(:,:,ii),LineStyle="none");
    title(ii);
    axis equal
    bar = colorbar();
    ylabel(bar,'|E_{tot}| (V/m)');
    view([0 0 1]);
    if db_axis
        set(gca,'ColorScale','log')
        clim([E_max*1e-3 E_max]);
    else
        clim([-E_max E_max]);
    end

    figure(2)
    colormap jet
    surf(Hz(:,:,ii),LineStyle="none");
    title(ii);
    axis equal
    bar = colorbar();
    ylabel(bar,'|H_{tot}| (A/m)');
    view([0 0 1]);
    if db_axis
        set(gca,'ColorScale','log')
        clim([H_max*1e-3 H_max]);
    else
        clim([-H_max H_max]);
    end

    pause(0.001)
end