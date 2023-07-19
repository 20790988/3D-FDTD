clear
close all

%====================SETTINGS=====================%
db_axis = true;
%=========================================%
result_filename = "full_field.mat";

%Line properties
    %distance between plates and width of line in meters
    

monitor = load(result_filename);
monitor_values = monitor.monitor_values;

delta_t = monitor.delta_t;
delta_x = monitor.delta_z;

Ez = monitor_values{1};

N = length(Ez);

num_monitors = length(monitor_values);

figure(1)

E_max = max(Ez,[],"all");


for ii = 1:N
    colormap jet
    surf(abs(Ez(:,:,ii)),LineStyle="none");
    title(ii);
    axis equal
    bar = colorbar();
    ylabel(bar,'|A_{tot}| (V/m)');
    view([0 0 1]);
    if db_axis
        set(gca,'ColorScale','log')
    end
    clim([1e-6 E_max]);
    pause(0.005)
end

t = (0:N-1)*delta_t;
figure(1)
plot(t./1e-12,voltage);
grid on
xlabel('Time (ps)')
ylabel('Voltage')