function [param, grid, source] = import_model()
    
    param = struct('material',0);
    source = struct('coord',0);
    source.value = {0,0,@source_func};
    
% Material specification
    sigma = [0 60e6 0];
    sigma_m = [0 0 0];    
    epsilon_r = [1 1 2];
    mu_r = [1 1 1];
    
    param.material = [sigma;sigma_m;epsilon_r;mu_r];

% Material at Border
    param.border_material_index = 1;

% Grid and cell size

    M_x = 300;
    M_y = 300;
    M_z = 300;

    unit = 1e-3;

    delta_x = 3e-3;
    delta_y = delta_x;
    delta_z = delta_x;

    param.delta = {delta_x,delta_y,delta_z};

% Simulation length
    param.N_t_max = 500;

% Model

    N = m_to_n(M_x,M_y,M_z,unit,param.delta);
    
    param.N = N;

    grid = ones(N{1},N{2},N{3});
    
    grid = add_cuboid(grid,param.delta,unit,200,280,30,280,30,280,3);
    
    source.coord = m_to_n(150,150,150,unit,param.delta);
end

function source_signal = source_func(t)

    freq = 10e9;
    source_signal = sin(2*pi*freq*t);
end


function grid = add_cuboid(grid,delta,unit,xmin,xmax,ymin,ymax,zmin,zmax,material)
   
    min = m_to_n(xmin,ymin,zmin,unit,delta);
    max = m_to_n(xmax,ymax,zmax,unit,delta);

    grid(min{1}:max{1},min{2}:max{2},min{3}:max{3}) = material;

end

function point = m_to_n(x,y,z,unit,delta)
%     align_to_grid = false;
    ii = x*unit/delta{1}+1;
    jj = y*unit/delta{2}+1;
    kk = z*unit/delta{3}+1;
%     if (align_to_grid == true && (~is_int(ii) || ~is_int(jj) || is_int(kk)))
%         fprintf('Model could not be aligned to grid')
%         return;
%     end

    point = {floor(ii),floor(jj),floor(kk)};
end

function result = is_int(num)
    error = 1e-3;
    fl = floor(num);
    diff = fl-num;
    result = false;
    if (diff<0 && diff>error) || (diff>0 && diff<error)
        result = true;
    end
end

