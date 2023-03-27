function [source_x, source_y, source_z, source_signal] = import_source(N_x,N_y,N_z,N_t_max,delta_t)
    source_x = 25;
    source_y = N_y/2;
    source_z = N_z/2;
    
    t = 0:delta_t:N_t_max*delta_t;
    freq = 10e9;
    source_signal = sin(2*pi*freq*t);
end