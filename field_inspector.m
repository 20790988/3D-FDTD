clear -except monitor
close all

%====================SETTINGS=====================%
db_axis = false;
%=========================================%
result_filename = "field_cap.mat";

%Line properties
    %distance between plates and width of line in meters
    

monitor = load(result_filename);
monitor_values = monitor.monitor_values{1};

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

figure(1)

E_max = max(Ez,[],"all");


for ii = 1:N
    if sum(Ez(:,:,ii),"all") == 0
        continue;
    end
    colormap jet
    surf(Ez(:,:,ii),LineStyle="none");
    title(ii);
    axis equal
    bar = colorbar();
    ylabel(bar,'|A_{tot}| (V/m)');
    view([0 0 1]);
    if db_axis
        set(gca,'ColorScale','log')
        clim([1e-6 E_max]);
    else
        clim([-E_max E_max]);
    end
    pause(0.001)
end