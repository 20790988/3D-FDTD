clear -except monitor
close all

%====================SETTINGS=====================%
db_axis = true;
%=========================================%
result_filename = "field_cap_2.mat";

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

% figure(2)
% s = size(Ez);
% plot(squeeze(Ez(floor(s(1)/2),floor(s(2)/2),:)));

figure(1)

Ez = sqrt(Hx.^2+Hy.^2+Hz.^2);

E_max = max(abs(Ez),[],"all");

Ez = permute(Ez,[2 1 3]);
for ii = 700:N
    if sum(Ez(:,:,ii),"all") == 0
        continue;
    end
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

    pause(0.001)
end