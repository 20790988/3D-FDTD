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
    
    delta_x = 0.06;
    delta_y = delta_x;
    delta_z = delta_x;
    
    %Grid size and variables

    M_x = 13.38;
    M_y = 6.96;
    M_z = 2.4;

    % Grid alignment behaviour
    grid_pause_on_unaligned = true;
    grid_error_tolerance = 0.5;

    % Simulation length in seconds
    param.M_t_max = 0.8e-9;
    
%============================================================%

    param.material = [sigma;sigma_m;epsilon_r;mu_r];
    delta = {delta_x,delta_y,delta_z};
    param.delta = {delta_x*unit,delta_y*unit,delta_z*unit};

    N = m_to_n(M_x,M_y,M_z,delta,{0,0,0});   
    param.N = {N{1},N{2},N{3}};

    grid = ones(N{1},N{2},N{3});
    
%====================MODEL SETUP====================%
    grid(:,:,:) = FREE_SPACE;
    
    origin = {1.98,M_y/2,0.12};

    %conductor
     grid = add_cuboid(grid,delta,-1.98,7.2, ...
        -0.3,0.3, ...
        0.6,0.72, ...
        PEC, ...
        origin);

    %ground plane
    grid = add_cuboid(grid,delta,0,M_x, ...
        0,M_y, ...
        0,0.12, ...
        PEC, ...
        {0,0,0});

grid = add_cuboid(grid,delta,-1.98,M_x-1.98, ...
        -M_y/2,M_y/2, ...
        0.06,0.54, ...
        DIELECTRIC, ...
        origin);
  
%====================SOURCE PROPERTIES====================%
    %if t_max = 0 (default), source will continue as long as simulation
   
    source.coord = m_to_n(0, -0.3:delta_y:0.3, 0.06:delta_z:0.54, delta, origin);

%     source.t_max = 0.6e-9;

%====================MONITOR SETUP====================%
    for ii = 1:7
        str = sprintf('port_%d',ii);
        monitor(ii).name = str;
        monitor(ii).coords = m_to_n(ii*0.96, M_y/2, 0.12:delta_z:0.72, delta,{0,0,0});
    end
   
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
function [source_signal_E,source_signal_H] = source_func(t,delta_t)


    epsilon_0 = 8.8542e-12;
    mu_0 = 1.2566e-6;
    e_eff = 6.493;
    eta = sqrt(mu_0./(epsilon_0*e_eff));
    a = 1/(0.6e-3);

    t0 = 36e-12;
    T = 4.115e-12;
    source_signal_E = -exp(-(t-t0).^2./(T^2))*a;
    source_signal_H = exp(-(t-t0).^2./(T^2))*a/eta;

%     t = t-35e-12;
%     source_signal_E = -gauspuls(t,40e9,1)*a;
%     
%     t = t-(0.366*delta_t);
%     source_signal_H = gauspuls(t,40e9,1)*a/eta;
end
%==================================================%

function grid = add_cuboid(grid,delta,xmin,xmax,ymin,ymax,zmin,zmax,material,relative_to)
    
    if xmin>xmax
        [xmin,xmax] = swop(xmin,xmax);
    end

    if ymin>ymax
        [ymin,ymax] = swop(ymin,ymax);
    end

    if zmin>zmax
        [zmin,zmax] = swop(zmin,zmax);
    end
      
    min = m_to_n(xmin,ymin,zmin,delta,relative_to);
    max = m_to_n(xmax,ymax,zmax,delta,relative_to);

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


function point = m_to_n(x,y,z,delta,relative_to)
      
    ii = (x+relative_to{1})/delta{1}+1;
    jj = (y+relative_to{2})/delta{2}+1;
    kk = (z+relative_to{3})/delta{3}+1;
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
