function [param, grid, source, monitor] = model_waveguide()
    
    param = struct('material',0);
    source = struct('coord',0);
    source.t_max = 0;
    source.value = {0,0,@source_func};
    
    global grid_pause_on_unaligned grid_error_tolerance grid_max_error
    global floor_tolerance
    floor_tolerance = 1e-9;
    grid_max_error = 0;
%====================SIMULATION PARAMETERS====================%

    % Material specification
    sigma = [0 0];
    sigma_m = [0 0];    
    epsilon_r = [1 9.6];
    mu_r = [1 1];
    
    PEC = 0;
    FREE_SPACE = 1;
    DIELECTRIC = 2;

    % Material at Border
    param.border_material_index = 1;
    
    % Cell size in units 
    unit = 1e-3;
    
    delta_x = 0.1;
    delta_y = delta_x;
    delta_z = delta_x;
    
    %Grid size and variables
    w = 0.6;
    th = 0.2;
    port_offset = 1;

    M_x = 20;
    M_y = 6;
    M_z = 6;
    
    % Grid alignment behaviour
    grid_pause_on_unaligned = true;
    grid_error_tolerance = 0.5;

    % Simulation length in seconds
    param.M_t_max = 400e-12;
    
%============================================================%

    param.material = [sigma;sigma_m;epsilon_r;mu_r];
    delta = {delta_x,delta_y,delta_z};
    param.delta = {delta_x*unit,delta_y*unit,delta_z*unit};

    N = m_to_n(M_x,M_y,M_z,delta);   
    param.N = {N{1},N{2},N{3}};

    grid = ones(N{1},N{2},N{3});

%====================MODEL SETUP====================%
    grid(:,:,:) = FREE_SPACE;


    grid = add_cuboid(grid,delta,0,M_x, ...
        M_y/2-w/2,M_y/2+w/2, ...
        M_z/2-w/2-th,M_z/2-w/2, ...
        PEC);


    grid = add_cuboid(grid,delta,0,M_x, ...
        M_y/2-w/2,M_y/2+w/2, ...
        M_z/2+w/2+th+0.1,M_z/2+w/2+0.1, ...
        PEC);

%     grid = add_cuboid(grid,delta,0,M_x, ...
%         M_y/2-w/2,M_y/2+w/2, ...
%         2.8,3.3, ...
%         DIELECTRIC);
  
%====================SOURCE PROPERTIES====================%
    %if t_max = 0 (default), source will continue as long as simulation
   
    source.coord = m_to_n(port_offset, 2.8:delta_y:3.3, 2.8:delta_z:3.3, delta);

%     source.t_max = 0.6e-9;

%====================MONITOR SETUP====================%
    monitor(1).name = 'port_1';
    monitor(1).coords = m_to_n(port_offset-0.1, M_y/2, 2.8:delta_z:3.3, delta);

    monitor(2).name = 'port_2';
    monitor(2).coords = m_to_n(M_x-port_offset, M_y/2, 2.8:delta_z:3.3, delta);

    monitor(3).name = 'port_ref';
    monitor(3).coords = m_to_n(port_offset+0.1, M_y/2, 2.8:delta_z:3.3, delta);

%==================================================%

    if grid_max_error < grid_error_tolerance
        fprintf('Model setup sucsessful\n');  
        fprintf('Max gridding error:%.3e of cell\n',grid_max_error);
    else
        fprintf('Gridding error larger than tolerance:%.3e of cell\n',grid_max_error);
        if grid_pause_on_unaligned
            param = 0;
            grid = 0;
            source = 0;
            fprintf('Model setup unsucsessful\n');
        end
    end
end

%====================SOURCE SIGNAL====================%
function source_signal = source_func(t)
%     t0 = 0.3e-9;
%     T = 7.7032e-11;
%     source_signal = exp(-(t-t0).^2./(T^2));

    t = t-60e-12;
    a = 1/(0.6e-3);
    source_signal = gauspuls(t,25e9,1)*a;

%     f = 5e9;pack
%     source_signal = sin(2*pi*f*t);
end
%==================================================%

function grid = add_cuboid(grid,delta,xmin,xmax,ymin,ymax,zmin,zmax,material)
    
    if xmin>xmax
        [xmin,xmax] = swop(xmin,xmax);
    end

    if ymin>ymax
        [ymin,ymax] = swop(ymin,ymax);
    end

    if zmin>zmax
        [zmin,zmax] = swop(zmin,zmax);
    end

    min = m_to_n(xmin,ymin,zmin,delta);
    max = m_to_n(xmax,ymax,zmax,delta);

    grid(min{1}:max{1},min{2}:max{2},min{3}:max{3}) = material;

end

function grid = add_tube(grid,delta,direction,xmin,xmax,ymin,ymax,zmin,zmax,material)
   
    min = m_to_n(xmin,ymin,zmin,delta);
    max = m_to_n(xmax,ymax,zmax,delta);
    
    grid(min{1},min{2}:max{2},min{3}:max{3}) = material;
    grid(max{1},min{2}:max{2},min{3}:max{3}) = material;
    grid(min{1}:max{1},min{2},min{3}:max{3}) = material;
    grid(min{1}:max{1},max{2},min{3}:max{3}) = material;

end


function point = m_to_n(x,y,z,delta)
  
    ii = x/delta{1}+1;
    jj = y/delta{2}+1;
    kk = z/delta{3}+1;
    if (~is_int(ii) || ~is_int(jj) || ~is_int(kk))
        fprintf('Gridding error at dimension\n%.3e %.3e %.3e\n',x,y,z);
    end

    point = {floor_(ii),floor_(jj),floor_(kk)};
end

function result = is_int(num)
    global grid_error_tolerance grid_max_error
    fl = floor_(num);
    diff = max(abs(fl-num));
    result = true;
    if (diff>grid_error_tolerance)
        result = false;
    end

    if diff> grid_max_error
        grid_max_error = diff;
    end
end

function [b,a] = swop(a,b)
    return
end

function result = floor_(num)
global floor_tolerance

    result = floor(num+floor_tolerance);
end
