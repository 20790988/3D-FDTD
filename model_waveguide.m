function [param, grid, source] = model_waveguide()
    
    param = struct('material',0);
    source = struct('coord',0);
    source.value = {@source_func,@source_func,@source_func};
    
    global grid_pause_on_unaligned grid_error_tolerance grid_max_error
    grid_max_error = 0;
%====================SIMULATION PARAMETERS====================%

    % Material specification
    sigma = [0 0];
    sigma_m = [0 0];    
    epsilon_r = [1 2];
    mu_r = [1 1];
    
    PEC = 0;
    FREE_SPACE = 1;
    DIELECTRIC = 2;

    % Material at Border
    param.border_material_index = 1;
    
    % Grid and cell size in units 
    unit = 1e-3;

    M_x = 51;
    M_y = 21;
    M_z = 500;

    delta_x = 1;
    delta_y = delta_x;
    delta_z = delta_x;
      
    % Grid alignment behaviour
    grid_pause_on_unaligned = false;
    grid_error_tolerance = 1;

    % Simulation length in seconds
    param.M_t_max = 100e-9;
    
%============================================================%

    param.material = [sigma;sigma_m;epsilon_r;mu_r];
    delta = {delta_x,delta_y,delta_z};
    param.delta = {delta_x*unit,delta_y*unit,delta_z*unit};

    N = m_to_n(M_x,M_y,M_z,delta);   
    param.N = {N{1},N{2},N{3}};

    grid = ones(N{1},N{2},N{3});

%====================MODEL SETUP====================%
    
    grid(:,:,:) = FREE_SPACE;

    grid = add_tube_z(grid,delta,3,48,3,18,0,M_z,PEC);
  
%====================SOURCE POSITION====================%

    source.coord = m_to_n([5:delta_x:46],[5:delta_y:16],250,delta);

%==================================================%
    if grid_max_error < grid_error_tolerance
        fprintf('Model setup sucsessful\n');    
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
    freq = 5e9;
    source_signal = sin(2*pi*freq*t);
end
%==================================================%

function grid = add_cuboid(grid,delta,xmin,xmax,ymin,ymax,zmin,zmax,material)
   
    min = m_to_n(xmin,ymin,zmin,delta);
    max = m_to_n(xmax,ymax,zmax,delta);

    grid(min{1}:max{1},min{2}:max{2},min{3}:max{3}) = material;

end

function grid = add_tube_z(grid,delta,xmin,xmax,ymin,ymax,zmin,zmax,material)
   
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

    point = {floor(ii),floor(jj),floor(kk)};
end

function result = is_int(num)
    global grid_error_tolerance grid_max_error
    fl = floor(num);
    diff = abs(fl-num);
    result = true;
    if (diff>grid_error_tolerance)
        result = false;
    end

    if diff> grid_max_error
        grid_max_error = diff;
    end
end

