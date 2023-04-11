function [param, grid, source] = model_waveguide()
    
    param = struct('material',0);
    source = struct('coord',0);
    source.t_max = 0;
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
    M_y = 50+40;
    M_z = 20+40;

    delta_x = 1;
    delta_y = delta_x;
    delta_z = delta_x;
      
    % Grid alignment behaviour
    grid_pause_on_unaligned = false;
    grid_error_tolerance = 1;

    % Simulation length in seconds
    param.M_t_max = 3e-9;
    
%============================================================%

    param.material = [sigma;sigma_m;epsilon_r;mu_r];
    delta = {delta_x,delta_y,delta_z};
    param.delta = {delta_x*unit,delta_y*unit,delta_z*unit};

    N = m_to_n(M_x,M_y,M_z,delta);   
    param.N = {N{1},N{2},N{3}};

    grid = ones(N{1},N{2},N{3});

%====================MODEL SETUP====================%
    
    grid(:,:,:) = FREE_SPACE;
    ofs = 20;
    grid = add_cuboid(grid,delta,ofs,M_x-ofs,0+ofs,M_y-ofs,0+ofs,0+ofs,PEC);
    grid = add_cuboid(grid,delta,ofs,M_x-ofs,0+ofs,M_y-ofs,M_z-ofs,M_z-ofs,PEC);
  
%====================SOURCE PROPERTIES====================%
    %if t_max = 0 (default), source will continue as long as simulation

    source.coord = m_to_n(M_x/2, ofs+1:delta_y:M_y-ofs-1, ofs+1:delta_z:M_z-ofs-1, delta);
%     source.t_max = 7.7032e-10;

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
    t0 = 6.7403e-10;
    T = 7.7032e-11;
    source_signal = exp(-(t-t0).^2./(T^2));
end
%==================================================%

function grid = add_cuboid(grid,delta,xmin,xmax,ymin,ymax,zmin,zmax,material)
   
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

