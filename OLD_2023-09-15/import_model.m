function [param, grid, source] = import_model()
    
    param = struct('material',0);
    source = struct('coord',0);
    source.value = {0,0,@source_func};
    
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

    M_x = 300;
    M_y = 300;
    M_z = 300;

    delta_x = 4e-3;
    delta_y = delta_x;
    delta_z = delta_x;
      
    % Grid alignment behaviour
    grid_pause_on_unaligned = false;
    grid_error_tolerance = 1;

    % Simulation length in seconds
    param.M_t_max = 2e-9;
    
%============================================================%

    param.material = [sigma;sigma_m;epsilon_r;mu_r];
    param.delta = {delta_x,delta_y,delta_z};

    N = m_to_n(M_x,M_y,M_z,unit,param.delta);   
    param.N = N;
    grid = ones(N{1},N{2},N{3});

%====================MODEL SETUP====================%
    
    grid(:,:,:) = FREE_SPACE;
    grid = add_cuboid(grid,param.delta,unit,180,240,75,225,120,180,PEC);
  
%====================SOURCE POSITION====================%

    source.coord = m_to_n(150,[60:delta_x:240],150,unit,param.delta);

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
    freq = 10e9;
    source_signal = sin(2*pi*freq*t);
end
%==================================================%

function grid = add_cuboid(grid,delta,unit,xmin,xmax,ymin,ymax,zmin,zmax,material)
   
    min = m_to_n(xmin,ymin,zmin,unit,delta);
    max = m_to_n(xmax,ymax,zmax,unit,delta);

    grid(min{1}:max{1},min{2}:max{2},min{3}:max{3}) = material;

end

function point = m_to_n(x,y,z,unit,delta)
  
    ii = x*unit/delta{1}+1;
    jj = y*unit/delta{2}+1;
    kk = z*unit/delta{3}+1;
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

