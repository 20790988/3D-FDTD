function [param, grid, source, monitor] = model_waveguide()
    % Superconductor TM line

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
    sigma = [0 0 0];
    sigma_m = [0 0 0];    
    epsilon_r = [1 4.6 1];
    mu_r = [1 1 1];
    
    PEC = 0;
    FREE_SPACE = 1;
    DIELECTRIC = 2;
    SUPERCONDUCTOR = 3;

    % Material at Border
    param.border_material_index = 1;

    % Superconductivity
        %0 = no superconductivity
        %1 = two-fluid only
        param.sc_model_level = 1;

        %lambda in m
        param.lambda_L = 90e-9;
        %sigma in S/m
        param.sigma_n = 6.7e6;
        %temps in K
        param.T_op = 4.2;
        param.T_c = 9.3;


    
    % Cell size in units 
    unit = 1e-6;
    
    delta_x = 0.05;
    delta_y = delta_x;
    delta_z = delta_x;
    
    %Grid size and variables

    l_ = 10;

    M_x = l_;
    M_y = 12;
    M_z = 2;
    
    % Grid alignment behaviour
    grid_pause_on_unaligned = true;
    grid_error_tolerance = 1;

    % Simulation length in seconds
    param.M_t_max = 80e-15;
    
%============================================================%

    param.material = [sigma;sigma_m;epsilon_r;mu_r];
    delta = {delta_x,delta_y,delta_z};
    param.delta = {delta_x*unit,delta_y*unit,delta_z*unit};

    N = m_to_n(M_x,M_y,M_z,delta,{0,0,0});   
    param.N = {N{1},N{2},N{3}};

    grid = ones(N{1},N{2},N{3});
    
%====================MODEL SETUP====================%
    grid(:,:,:) = DIELECTRIC;
    
    origin = {0,1,0.5};

    %M2
     grid = add_cuboid(grid,delta,0,l_, ...
        0,10, ...
        0,0.2, ...
        SUPERCONDUCTOR, ...
        origin);

     %M3
     grid = add_cuboid(grid,delta,0,l_, ...
        2.75,7.25, ...
        0.4,0.6, ...
        SUPERCONDUCTOR, ...
        origin);

     %M4
     grid = add_cuboid(grid,delta,0,l_, ...
        0,10, ...
        0.8,1.0, ...
        SUPERCONDUCTOR, ...
        origin);

     grid = routeblock(grid,delta,SUPERCONDUCTOR,origin,{0,0,0});
%      grid = routeblock(grid,delta,SUPERCONDUCTOR,origin,{10,0,0});

%====================SOURCE PROPERTIES====================%
    %if t_max = 0 (default), source will continue as long as simulation
    a = 0.2+delta_z;
    b = 0.4-delta_z;
    source.coord = m_to_n(2.5, 2.75:delta_y:7.25, a:delta_z:b, delta, origin);

%     source.t_max = 0.6e-9;

%====================MONITOR SETUP====================%
    N_temp = 9;

    for ii = 1:N_temp
        str = sprintf('port_1_%dmm',ii);
        monitor(ii).name = str;
        monitor(ii).coords = m_to_n(ii, 5, 0.2:delta_z:0.4, delta, origin);
        monitor(ii).normal_direction = 1;
    end
   
%==================================================%
    if sum(grid==SUPERCONDUCTOR,'all') > 0 && param.sc_model_level == 0
        fprintf('Model contains superconductor but model level == 0\n');
        param = 0;
        grid = 0;
        source = 0;
        fprintf('Model setup unsucsessful\n');
    end

   if sum(grid==SUPERCONDUCTOR,'all') == 0 && param.sc_model_level > 0
        fprintf('Supeconducting model ~= 0 but no superconductor in model\n');
        fprintf('This will increase runtime\n');
    end
    

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
function [source_signal_E,source_signal_H] = source_func(t,delta_t,delta_x)


    epsilon_0 = 8.8542e-12;
    mu_0 = 1.2566e-6;
    e_eff = 4.6;
    eta = sqrt(mu_0./(epsilon_0*e_eff));
    a = 2e-3;

%     t0 = 36e-12;
%     T = 4.115e-12;
%     source_signal_E = -exp(-(t-t0).^2./(T^2))*a;
%     source_signal_H = exp(-(t-t0).^2./(T^2))*a/eta;
    
   
    source_signal_E = -gaus_derv(t,20e-15,5e-15)*a;
    
    source_signal_H = gaus_derv(t,20e-15-0.366*delta_t,5e-15)*a/eta;


%     t = t-35e-12;
%     source_signal_E = -gauspuls(t,40e9,1)*a;
%     
%     t = t-(0.366*delta_t);
%     source_signal_H = gauspuls(t,40e9,1)*a/eta;
end

function pulse = gaus_derv(t,mu,sigma)
    pulse = (t-mu).*exp(-((t-mu).^2)/(2*sigma^2));
    pulse =  pulse./max(pulse,[],'all');
end

%==================================================%

function grid = routeblock(grid,delta,material,global_origin,desired_placement)
    
    origin = sumcell(global_origin,desired_placement);
    
    lside = 10;
    lvia = 0.6;
    lcon = 1.25;
    loff = 0.35;
    
    %M3
    coords = {{0,0,0}, ...
            {0,lside-lcon,0}, ...
            {lside-lcon,lside-lcon,0},...
            {lside-lcon,0,0}};
    
    %M3
    for i = 1:4
        temp_origin = sumcell(origin,coords{i});

        grid = add_cuboid(grid,delta,0,lcon, ...
            0,lcon, ...
            0.4,0.6, ...
            material, ...
            temp_origin);

        if i==1 || i==3
            z1 = 0.2;
            z2 = 0.4;
        else
            z1 = 0.6;
            z2 = 0.8;
        end
        
        grid = add_cuboid(grid,delta,loff,loff+lvia, ...
            loff,loff+lvia, ...
            z1,z2, ...
            material, ...
            temp_origin);
    end
end

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

function result = sumcell(a,b)
    for i = 1:length(a)
        result{i} = a{i}+b{i};
    end
end

