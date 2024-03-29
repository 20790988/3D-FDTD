function [param, grid, source, monitor] = s11_Via()
    % Superconductor TM line

    param = struct('material',0);
%     source = struct('coord',0);
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

    % Superconductivity
        %0 = no superconductivity
        %1 = two-fluid
        param.sc_model_level = 1;

        %lambda in m
        param.lambda_L_0 = 0;
        param.lambda_L = 90e-9;
        %sigma in S/m
        param.sigma_n = 6.7e6;
        %temps in K
        param.T_op = 4.2;
        param.T_c = 9.3;
    
    unit = 1e-6;

    % Cell size in units 
    delta_x = 0.05;
    delta_y = delta_x;
    delta_z = delta_x;
    
    %Grid size and variables

    l_ = 30;
    w_ = 10;
    
    wSig = 4.4;
    wGP = 10;

    hSim = 2.2;

    tLine = 0.2;

    M_x = l_;
    M_y = w_;
    M_z = hSim;
    
    % Grid alignment behaviour
    grid_pause_on_unaligned = true;
    grid_error_tolerance = 1;

    % Simulation length in seconds
    param.M_t_max = 800e-15;

    % Field capture
    field_capture = false;
    field_cap_normal_direction = 1;
    field_cap_x = 12;
    field_cap_y = 0:delta_y:M_y-1*delta_y;
    field_cap_z = 0:delta_z:M_z-1*delta_z;
    field_cap_fields = [0,1,1,0,1,1];
%     Ex Ey Ez Hx Hy Hz
    
%============================================================%

    param.material = [sigma;sigma_m;epsilon_r;mu_r];
    delta = {delta_x,delta_y,delta_z};
    param.delta = {delta_x*unit,delta_y*unit,delta_z*unit};

    N = m_to_n(M_x-delta_x,M_y-delta_y,M_z-delta_z,delta,{0,0,0});   
    param.N = {N{1},N{2},N{3}};

    grid = ones(N{1},N{2},N{3});
    
%====================MODEL SETUP====================%
    grid(:,:,:) = DIELECTRIC;
    
    origin = {0,0,4*delta_x};

     grid = routeblock(grid,delta,SUPERCONDUCTOR,DIELECTRIC,origin,{0,0,0});
     grid = routeblock(grid,delta,SUPERCONDUCTOR,DIELECTRIC,origin,{10,0,0});
     grid = routeblock(grid,delta,SUPERCONDUCTOR,DIELECTRIC,origin,{20,0,0});

     grid = viablock(grid,delta,SUPERCONDUCTOR,DIELECTRIC,origin,{10,0,0});


     %M3
     grid = add_cuboid(grid,delta,0,15, ...
        wGP/2-wSig/2,wGP/2+wSig/2, ...
        6*tLine,7*tLine, ...
        SUPERCONDUCTOR, ...
        origin);

     %M3
     grid = add_cuboid(grid,delta,15,30, ...
        wGP/2-wSig/2,wGP/2+wSig/2, ...
        2*tLine,3*tLine, ...
        SUPERCONDUCTOR, ...
        origin);
    
     %PEC top & bottom
     grid = add_cuboid(grid,delta,0,30, ...
        0,wGP, ...
        0,-1*tLine, ...
        PEC, ...
        origin);

      grid = add_cuboid(grid,delta,0,30, ...
        0,wGP, ...
        9*tLine,10*tLine, ...
        PEC, ...
        origin);



%====================SOURCE PROPERTIES====================%
    %if t_max = 0 (default), source will continue as long as simulation
    source.t_max = 0;

    source_x{1} = 2;
    source_y{1} = [wGP/2-wSig/2, wGP/2+wSig/2];
    source_z{1} = [5 6].*tLine;

    source_x{2} = 2;
    source_y{2} = [wGP/2-wSig/2, wGP/2+wSig/2];
    source_z{2} = [7 8].*tLine;

    is_microstrip = false;

    d = 0.3;
    w = 4;
    param.e_eff_0 = epsilon_r(2);
    

%============================================================%
    if is_microstrip
        param.e_eff_0 = (param.e_eff_0+1)/2+(param.e_eff_0-1)/(2*sqrt(1+12*d/w));
    end
    
    for ind = 1:length(source_x)
        source.coord{ind} = m_to_n(source_x{ind}, ...
            (source_y{ind}(1)):delta_y:(source_y{ind}(2)-1*delta_y), ...
            (source_z{ind}(1)):delta_z:(source_z{ind}(2)-1*delta_z), delta, origin);
    end

%====================MONITOR SETUP====================%
    if field_capture 

        str = sprintf('field_capture');
        monitor(1).name = str;
        monitor(1).coords = m_to_n(field_cap_x, field_cap_y, field_cap_z, ...
            delta, {0,0,0});
        monitor(1).normal_direction = field_cap_normal_direction;
        monitor(1).fields_to_monitor = field_cap_fields;
    else

        N_temp = 11; %#ok<*UNRCH> 
        
        monitor_y = [wGP/2,wGP/2+delta_y];
        monitor_z = [7 8].*tLine;
    
        for ii = 1:N_temp
            str = sprintf('port_%dmm',ii);
            monitor(ii).name = str;
            monitor(ii).coords = m_to_n(ii, ...
                (monitor_y(1)):delta_y:(monitor_y(2)-1*delta_y), ...
                (monitor_z(1)):delta_z:(monitor_z(2)-1*delta_z), delta, origin);
            monitor(ii).normal_direction = 1;
            monitor(ii).fields_to_monitor = [0,0,1,0,1,0];
        end


        monitor_z = [3 4].*tLine;
        for ii = 1:N_temp
            x_ = ii+14;
            jj = ii+N_temp;
            str = sprintf('port_%dmm',x_);
            monitor(jj).name = str;
            monitor(jj).coords = m_to_n(x_, ...
                (monitor_y(1)):delta_y:(monitor_y(2)-1*delta_y), ...
                (monitor_z(1)):delta_z:(monitor_z(2)-1*delta_z), delta, origin);
            monitor(jj).normal_direction = 1;
            monitor(jj).fields_to_monitor = [0,0,1,0,1,0];
        end
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
function [source_signal_E,source_signal_H] = source_func(t,delta_t,delta_x,e_eff_0,source_index)
    
    EXP = 1;
    GAUSDERV = 2;
    GAUSPULSE = 3;

    signal_type = EXP;
    t0 = 20e-15;
    T = 5e-15;
    

    %========================================%

    epsilon_0 = 8.8542e-12;
    mu_0 = 1.2566e-6;
    eta = sqrt(mu_0./(epsilon_0*e_eff_0));
    a = 1;
    
    td = sqrt(mu_0.*(epsilon_0*e_eff_0))*0.5*delta_x-0.5*delta_t;

    if signal_type == EXP
        source_signal_E = -exp(-(t-t0).^2./(T^2))*a;
        source_signal_H = exp(-(t-t0- td).^2./(T^2))*a/eta;
    elseif signal_type == GAUSDERV
        source_signal_E = -gaus_derv(t,t0,T)*a;
        source_signal_H = gaus_derv(t-td,t0,T)*a/eta;
    elseif signal_type == GAUSPULSE
        t = t-t0;
        source_signal_E = -gauspuls(t,f0,bw)*a;
        source_signal_H = gauspuls(t- td,f0,bw)*a/eta;
    end

    if source_index == 2
        source_signal_E = -source_signal_E;
        source_signal_H = -source_signal_H;
    end
end

function pulse = gaus_derv(t,mu,sigma)
    pulse = (t-mu).*exp(-((t-mu).^2)/(2*sigma^2));
    pulse =  pulse./max(pulse,[],'all');
end

function grid = routeblock(grid,delta,material_cond,material_die,global_origin,desired_placement)
    
    origin = sumcell(global_origin,desired_placement);
    %local origin is in left bottom corner

    delta_x = delta{1};

    s_side = 10;

    s_i_via = 0.6;
    l_i_off = 0.35;

    s_m_via = 1.25;

    s_m_cutout = 2.8;
    l_m_off = 0.4;

    t_height = 0.2;
    
    %cross fill
    % fill entire block
    % cut holes
    % fill side pillars


    %M0
    % cross fill

    %I0
    % thin side pillars top left bottom right

    %M1
    % fill side pillars

    %ii
    % thin side pillars top r bot l

    %M2
    % cross fill

    %I2
    % thin side pillars top l bot r

    %M3
    %  fill side pillars

    %I3
    % thin side pillars top r bot l

    %M4
    % cross fill

    
    

    %M0 = M2 = M4
    z_vals = [0 1; 4 5; 8 9].*t_height;
    
    for ii = 1:3
        % fill entire block
        grid = add_cuboid(grid,delta,0,s_side,0,s_side,z_vals(ii,1),z_vals(ii,2),material_cond,origin);
       
        % cut holes
        edges_vals = [l_m_off+delta_x, s_m_cutout-delta_x];
        
        grid = add_cuboid(grid,delta, ...
            edges_vals(1),edges_vals(2), ...
            edges_vals(1),edges_vals(2), ...
            z_vals(ii,1),z_vals(ii,2),material_die,origin);

         grid = add_cuboid(grid,delta, ...
            edges_vals(1),edges_vals(2), ...
            s_side-edges_vals(1),s_side-edges_vals(2), ...
            z_vals(ii,1),z_vals(ii,2),material_die,origin);

         grid = add_cuboid(grid,delta, ...
            s_side-edges_vals(1),s_side-edges_vals(2), ...
            edges_vals(1),edges_vals(2), ...
            z_vals(ii,1),z_vals(ii,2),material_die,origin);

        grid = add_cuboid(grid,delta, ...
            s_side-edges_vals(1),s_side-edges_vals(2), ...
            s_side-edges_vals(1),s_side-edges_vals(2), ...
            z_vals(ii,1),z_vals(ii,2),material_die,origin);
     
        % fill side pillars
        edges_vals = [0, s_m_via];
        
        grid = add_cuboid(grid,delta, ...
            edges_vals(1),edges_vals(2), ...
            edges_vals(1),edges_vals(2), ...
            z_vals(ii,1),z_vals(ii,2),material_cond,origin);

         grid = add_cuboid(grid,delta, ...
            edges_vals(1),edges_vals(2), ...
            s_side-edges_vals(1),s_side-edges_vals(2), ...
            z_vals(ii,1),z_vals(ii,2),material_cond,origin);

         grid = add_cuboid(grid,delta, ...
            s_side-edges_vals(1),s_side-edges_vals(2), ...
            edges_vals(1),edges_vals(2), ...
            z_vals(ii,1),z_vals(ii,2),material_cond,origin);

        grid = add_cuboid(grid,delta, ...
            s_side-edges_vals(1),s_side-edges_vals(2), ...
            s_side-edges_vals(1),s_side-edges_vals(2), ...
            z_vals(ii,1),z_vals(ii,2),material_cond,origin);
    end
    
    %M1 = M3
    z_vals = [2 3; 6 7].*t_height;
    
    for ii = 1:2
        edges_vals = [0, s_m_via];
        
        grid = add_cuboid(grid,delta, ...
            edges_vals(1),edges_vals(2), ...
            edges_vals(1),edges_vals(2), ...
            z_vals(ii,1),z_vals(ii,2),material_cond,origin);

         grid = add_cuboid(grid,delta, ...
            edges_vals(1),edges_vals(2), ...
            s_side-edges_vals(1),s_side-edges_vals(2), ...
            z_vals(ii,1),z_vals(ii,2),material_cond,origin);

         grid = add_cuboid(grid,delta, ...
            s_side-edges_vals(1),s_side-edges_vals(2), ...
            edges_vals(1),edges_vals(2), ...
            z_vals(ii,1),z_vals(ii,2),material_cond,origin);

        grid = add_cuboid(grid,delta, ...
            s_side-edges_vals(1),s_side-edges_vals(2), ...
            s_side-edges_vals(1),s_side-edges_vals(2), ...
            z_vals(ii,1),z_vals(ii,2),material_cond,origin);
        
    end

    %I layers
    
    z_vals = [1 2; 3 4; 5 6; 7 8].*t_height;
    edges_vals = [l_i_off l_i_off+s_i_via, ...
                s_side-l_i_off s_side-(l_i_off+s_i_via)];
    for ii = 1:4
        if ii == 1 || ii == 3
             grid = add_cuboid(grid,delta, ...
                edges_vals(1),edges_vals(2), ...
                edges_vals(1),edges_vals(2), ...
                z_vals(ii,1),z_vals(ii,2),material_cond,origin);

             grid = add_cuboid(grid,delta, ...
                edges_vals(3),edges_vals(4), ...
                edges_vals(3),edges_vals(4), ...
                z_vals(ii,1),z_vals(ii,2),material_cond,origin);
        else
            grid = add_cuboid(grid,delta, ...
                edges_vals(1),edges_vals(2), ...
                edges_vals(3),edges_vals(4), ...
                z_vals(ii,1),z_vals(ii,2),material_cond,origin);

             grid = add_cuboid(grid,delta, ...
                edges_vals(3),edges_vals(4), ...
                edges_vals(1),edges_vals(2), ...
                z_vals(ii,1),z_vals(ii,2),material_cond,origin);
        end

       
    end
end

function grid = viablock(grid,delta,material_cond,material_die,global_origin,desired_placement)
    
    origin = sumcell(global_origin,desired_placement);
    %local origin is in left bottom corner

    delta_x = delta{1};

    s_side = 10;

    s_i_via = 0.6;
    l_i_off = 0.35;

    s_m_via = 1.25;

    s_m_cutout = 2.8;
    l_m_off = 0.4;

    t_height = 0.2;
    
    s_center = 4.4;
    s_center_via = 1.2;
    l_center_off = 0.5;

    l_g_stitch = 2.3;

    %M1
    % center block
    
    %I1
    % t r b l
    
    %M2
    % center block
    % GP cutaway
    
    %I1
    % t l b l
    
    %M3
    % center block
    
    %I3
    % t r b l
    
    %M2 GP cutaway
    z_vals = [4 5].*t_height;
    edges_vals = [l_g_stitch+delta_x, s_side-l_g_stitch-delta_x];
    

    grid = add_cuboid(grid,delta, ...
        edges_vals(1),edges_vals(2), ...
        edges_vals(1),edges_vals(2), ...
        z_vals(1,1),z_vals(1,2),material_die,origin);

    %M1 = M2 = M3
    z_vals = [2 3; 4 5; 6 7].*t_height;
    edges_vals = [s_side/2-s_center/2, s_side/2+s_center/2];
    
    for ii = 1:3
        grid = add_cuboid(grid,delta, ...
        edges_vals(1),edges_vals(2), ...
        edges_vals(1),edges_vals(2), ...
        z_vals(ii,1),z_vals(ii,2),material_cond,origin);
    end
    
    %I layers
    z_vals = [3 4; 5 6].*t_height;
    edges_vals = [s_side/2-s_center/2+l_center_off, ...
        s_side/2-s_center/2+l_center_off+s_center_via, ...
        s_side/2+s_center/2-l_center_off, ...
        s_side/2+s_center/2-(l_center_off+s_center_via)];

    for ii = 1:2
        if ii == 2
             grid = add_cuboid(grid,delta, ...
                edges_vals(1),edges_vals(2), ...
                edges_vals(1),edges_vals(2), ...
                z_vals(ii,1),z_vals(ii,2),material_cond,origin);

             grid = add_cuboid(grid,delta, ...
                edges_vals(3),edges_vals(4), ...
                edges_vals(3),edges_vals(4), ...
                z_vals(ii,1),z_vals(ii,2),material_cond,origin);
        else
            grid = add_cuboid(grid,delta, ...
                edges_vals(1),edges_vals(2), ...
                edges_vals(3),edges_vals(4), ...
                z_vals(ii,1),z_vals(ii,2),material_cond,origin);

             grid = add_cuboid(grid,delta, ...
                edges_vals(3),edges_vals(4), ...
                edges_vals(1),edges_vals(2), ...
                z_vals(ii,1),z_vals(ii,2),material_cond,origin);
        end
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

    grid(min{1}:(max{1}-1),min{2}:(max{2}-1),min{3}:(max{3}-1)) = material;

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

