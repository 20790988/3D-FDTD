%
%  3D-FDTD
%  K Meyer
%  2023-04-20
% 

clearvars

fprintf('Start\n')
main_timer = tic;


%====================SIMULATION SETTINGS=====================%        
    [param, material, source, monitor] = s09_SFQLine;

    TEMP_BOOTSTRAP_OFFSET = 300;
    window_bootstrap = false;

    gpu_accel = true;
    should_plot_output = false;

    use_bootstrapped_fields = false;
    bootstrap_field_name = 'field_cap_SC_.mat';


% Space around E-field in grid cells
    offset = 3;
    bc_offset = 2;

% corner angle
    cos_a = cos(pi/4);
    
%====================MODEL IMPORT=====================%    
    sigma = param.material(1,:);
    sigma_m = param.material(2,:);
    
    epsilon_0 = 8.8542e-12;
    epsilon_r = param.material(3,:);
    epsilon = epsilon_0*epsilon_r;
    
    mu_0 = 1.2566e-6;
    mu_r = param.material(4,:);
    mu = mu_0*mu_r;
    
    superconducting_model = param.sc_model_level;
    NONE = 0;
    TWOFLUID = 1;

    lambda_L = param.lambda_L;
    lambda_L_0 = param.lambda_L_0;
    sigma_n = param.sigma_n;
    T_op = param.T_op;
    T_c = param.T_c;

    N = param.N;

    N_x = N{1};
    N_y = N{2};
    N_z = N{3};
    
    delta = param.delta;
    delta_x = delta{1};
    delta_y = delta{2};
    delta_z = delta{3};
    
    c = 1./sqrt(epsilon.*mu);
    delta_t = delta_x/(max(c)*sqrt(3));
    
% Simulation length
    N_t_max = floor(param.M_t_max/delta_t);

% Fix sigma
%     sigma = sqrt(sigma./(2*pi*10e9*mu))./delta_x;


% Source plane
    source_coords = source.coord;
    sftf_x = source_coords{1}{1};

    for ind = 1:length(source_coords)
        source_y{ind} = source_coords{ind}{2};
        source_z{ind} = source_coords{ind}{3};
    end

source_N_t_max = floor(source.t_max/delta_t);

%====================SIMULATION SETUP AND INITIALISE=====================%   

t = (0:N_t_max)*delta_t;

for ind = 1:length(source_coords);
    [source_val_E{ind},source_val_H{ind}] = source.value{3}(t,delta_t,delta_x,param.e_eff_0,ind);
end

% load bootstrap fields
if use_bootstrapped_fields
    fprintf('Loading bootstrap fields...\n') %#ok<*UNRCH> 

    bootstrap = load(bootstrap_field_name);
%     source_field_Ex = permute(cell2mat(bootstrap.monitor_values{1}(1),[3 1 2]){};
    source_field_Ey = permute(cell2mat(bootstrap.monitor_values{1}(2)),[3 1 2]);
    source_field_Ez = permute(cell2mat(bootstrap.monitor_values{1}(3)),[3 1 2]);

%     source_field_Hx = permute(bootstrap.monitor_values{1}(4),[3 1 2]);
    source_field_Hy = permute(cell2mat(bootstrap.monitor_values{1}(5)),[3 1 2]);
    source_field_Hz = permute(cell2mat(bootstrap.monitor_values{1}(6)),[3 1 2]);

    fprintf('Loaded successfully\n')
end

% check stability
S = max(c)*delta_t/delta_x;

if S>1
    fprintf('Simulation unstable');
    return
end

%superconducting constants
if lambda_L == 0
    lambda_L = lambda_L_0/(sqrt(1-(T_op/T_c)^4));
end

alp_s = 1;
gam_s = delta_t/(2*mu_0*lambda_L^2);

sigma(3) = sigma_n*(T_op/T_c)^4;

%FDTD parameters
if superconducting_model == NONE

    C_a_single = (1-(sigma.*delta_t)./(2.*epsilon)) ...
        ./(1+(sigma.*delta_t)./(2.*epsilon));
    C_b_single = (delta_t./(epsilon.*delta_x)) ...
        ./(1+(sigma.*delta_t)./(2.*epsilon));

elseif superconducting_model == TWOFLUID
    
    half_sum_gam = 0.5*[0 0 gam_s];
    
    C_a_single = (epsilon/delta_t - half_sum_gam- sigma/2) ...
        ./(epsilon/delta_t + half_sum_gam + sigma/2);

    C_b_single = (1/delta_x)./(epsilon/delta_t + half_sum_gam + sigma/2);
    C_c_single = -1./(epsilon/delta_t + half_sum_gam + sigma/2);

    C_c_single(1:2) = 0;
    C_c = C_c_single(material);
    
end

D_a_single = (1-(sigma_m.*delta_t)./(2.*mu))./(1+(sigma_m.*delta_t)./(2.*mu));
D_b_single = (delta_t./(mu.*delta_x))./(1+(sigma_m.*delta_t)./(2.*mu));

pec_mask = (material ~= 0);
material(material==0) = 1;

sc_mask = (material == 3);

% Don't update supercurrent inside boundary - redundant
% sc_mask([1:bc_offset,(N_x-bc_offset+1):N_x],:,:) = 0;
% sc_mask(:,[1:bc_offset,(N_y-bc_offset+1):N_y],:) = 0;
% sc_mask(:,:,[1:bc_offset,(N_z-bc_offset+1):N_z]) = 0;

% expand material matrices
C_a = C_a_single(material);
C_b = C_b_single(material);

D_a = D_a_single(material);
D_b = D_b_single(material);

c_ = c(material);

% matrix declaration
Ex_old_old = zeros(N_x,N_y,N_z);
Ey_old_old = zeros(N_x,N_y,N_z);
Ez_old_old = zeros(N_x,N_y,N_z);

Ex_old = zeros(N_x,N_y,N_z);
Ex_new = zeros(N_x,N_y,N_z);
Ey_old = zeros(N_x,N_y,N_z);
Ey_new = zeros(N_x,N_y,N_z);
Ez_old = zeros(N_x,N_y,N_z);
Ez_new = zeros(N_x,N_y,N_z);
Hx_old = zeros(N_x,N_y,N_z);
Hx_new = zeros(N_x,N_y,N_z);
Hy_old = zeros(N_x,N_y,N_z);
Hy_new = zeros(N_x,N_y,N_z);
Hz_old = zeros(N_x,N_y,N_z);
Hz_new = zeros(N_x,N_y,N_z);

sftf_x = floor(sftf_x);
% source_y = floor(source_y);
% source_z = floor(source_z);

Hz_inc = zeros(1,N_y,N_z);
Hy_inc = zeros(1,N_y,N_z);
Ey_inc = zeros(1,N_y,N_z);
Ez_inc = zeros(1,N_y,N_z);
  
%superconducting
J_sx_new = zeros(N_x,N_y,N_z);
J_sy_new = zeros(N_x,N_y,N_z);
J_sz_new = zeros(N_x,N_y,N_z);

% J_nx_new = zeros(N_x,N_y,N_z);
% J_ny_new = zeros(N_x,N_y,N_z);
% J_nz_new = zeros(N_x,N_y,N_z);

J_sx_old = zeros(N_x,N_y,N_z);
J_sy_old = zeros(N_x,N_y,N_z);
J_sz_old = zeros(N_x,N_y,N_z);

% J_nx_old = zeros(N_x,N_y,N_z);
% J_ny_old = zeros(N_x,N_y,N_z);
% J_nz_old = zeros(N_x,N_y,N_z);

% setup and preallocate monitors
num_monitors = length(monitor);

for num = 1:num_monitors
    [iii,jjj,kkk] = unpack_coords(monitor(num).coords);
    if monitor(num).normal_direction == 1
        rows = length(jjj);
        cols = length(kkk);
    elseif monitor(num).normal_direction == 2
        rows = length(iii);
        cols = length(kkk);
    else
        rows = length(iii);
        cols = length(jjj);
    end
    field_selector = monitor(num).fields_to_monitor;

    for count = 1:6
        if field_selector(count)
            monitor_values{num}{count} = zeros(rows,cols,N_t_max);
        else
            monitor_values{num}{count} = 0;
        end
    end
    
    monitor_names{num} = monitor(num).name;
end

step = 0;
stop_cond = false;

print_text_at_sim_start(N_x*N_y*N_z,delta_x,delta_t);

if gpu_accel
    pec_mask = gpuArray(pec_mask);
    sc_mask = gpuArray(sc_mask);

    C_a = gpuArray(C_a);
    C_b = gpuArray(C_b);
    
    if superconducting_model ~= 0
        C_c = gpuArray(C_c);
    end

    D_a = gpuArray(D_a);
    D_b = gpuArray(D_b);
    
    c_ = gpuArray(c_);
    
    Ex_old_old = gpuArray(Ex_old_old);
    Ey_old_old = gpuArray(Ey_old_old);
    Ez_old_old = gpuArray(Ez_old_old);
    
    Ex_old = gpuArray(Ex_old);
    Ex_new = gpuArray(Ex_new);
    Ey_old = gpuArray(Ey_old);
    Ey_new = gpuArray(Ey_new);
    Ez_old = gpuArray(Ez_old);
    Ez_new = gpuArray(Ez_new);
    Hx_old = gpuArray(Hx_old);
    Hx_new = gpuArray(Hx_new);
    Hy_old = gpuArray(Hy_old);
    Hy_new = gpuArray(Hy_new);
    Hz_old = gpuArray(Hz_old);
    Hz_new = gpuArray(Hz_new);

    Hz_inc = gpuArray(Hz_inc);
    Hy_inc = gpuArray(Hy_inc);
    Ey_inc = gpuArray(Ey_inc);
    Ez_inc = gpuArray(Ez_inc);
      
    %superconducting
    J_sx_new = gpuArray(J_sx_new);
    J_sy_new = gpuArray(J_sy_new);
    J_sz_new = gpuArray(J_sz_new);
    
    % J_nx_new = gpuArray();
    % J_ny_new = gpuArray();
    % J_nz_new = gpuArray();
    
    J_sx_old = gpuArray(J_sx_old);
    J_sy_old = gpuArray(J_sy_old);
    J_sz_old = gpuArray(J_sz_old);
  
    % J_nx_old = gpuArray();
    % J_ny_old = gpuArray();
    % J_nz_old = gpuArray();
end

%====================LOOP START=====================%

while stop_cond == false
    text_update(step,N_t_max,delta_t)

    if use_bootstrapped_fields
        Ey_inc(1,:,:) = (source_field_Ey(step+TEMP_BOOTSTRAP_OFFSET,:,:));
        Ez_inc(1,:,:) = (source_field_Ez(step+TEMP_BOOTSTRAP_OFFSET,:,:));
          
        Hy_inc(1,:,:) = (source_field_Hy(step+TEMP_BOOTSTRAP_OFFSET,:,:));
        Hz_inc(1,:,:) = (source_field_Hz(step+TEMP_BOOTSTRAP_OFFSET,:,:));
    else
        for ind = 1:length(source_val_E)
            Ez_inc(1,source_y{ind},source_z{ind}) = source_val_E{ind}(step+1);
            Hy_inc(1,source_y{ind},source_z{ind}) = source_val_H{ind}(step+1);
        end
    end
  
    %====================PLOTTING START=====================%
    if  should_plot_output
        
%         H_tot = sqrt(Hx_old.^2+Hy_old.^2+Hz_old.^2);
%          plot_field(H_tot,N_x/2,N_y/2,14,step,delta,delta_t);

%         tempz = floor(N_z/2);
%         tempy = floor(N_y/2);
%         H_tot_line = (H_tot(:,tempy,tempz));
%         plot_line(H_tot_line,delta_x*(0:N_x-1),step,'|H_{tot}| (A/m)',1,1/eta);
        
        tempz = floor(17);
        tempy = floor(N_y/2);

        E_tot = sqrt(Ex_old.^2+Ey_old.^2+Ez_old.^2);
        plot_field(E_tot,N_x/2,tempy,tempz,step,delta,delta_t);
%          view([1 0 0])

        E_tot_line = (E_tot(:,tempy,tempz));
        plot_line(E_tot_line,delta_x*(0:N_x-1),step,'|E_{tot}| (V/m)',2);
        temp = 0;
    end
    %====================PLOTTING END=====================%

    % H-field calculate
    hofs = offset-1;
    hnfs = offset-2;

    ii = hofs:N_x-hnfs;
    jj = hofs:N_y-hnfs;
    kk = hofs:N_z-hnfs;
      
    Hx_new(ii,jj,kk) = D_a(ii,jj,kk).*Hx_old(ii,jj,kk) ...
        + D_b(ii,jj,kk).*(Ey_old(ii,jj,kk+1)-Ey_old(ii,jj,kk) ...
        + Ez_old(ii,jj,kk) - Ez_old(ii,jj+1,kk));

    Hy_new(ii,jj,kk) = D_a(ii,jj,kk).*Hy_old(ii,jj,kk) ...
        + D_b(ii,jj,kk).*(Ez_old(ii+1,jj,kk)-Ez_old(ii,jj,kk) ...
        + Ex_old(ii,jj,kk) - Ex_old(ii,jj,kk+1));

    Hz_new(ii,jj,kk) = D_a(ii,jj,kk).*Hz_old(ii,jj,kk) ...
        + D_b(ii,jj,kk).*(Ex_old(ii,jj+1,kk)-Ex_old(ii,jj,kk) ...
        + Ey_old(ii,jj,kk) - Ey_old(ii+1,jj,kk));
   
    Hz_new(sftf_x,jj,kk) = Hz_new(sftf_x,jj,kk) ...
        + D_b(sftf_x,jj,kk).*Ey_inc(1,jj,kk);
    Hy_new(sftf_x,jj,kk) = Hy_new(sftf_x,jj,kk) ...
        - D_b(sftf_x,jj,kk).*Ez_inc(1,jj,kk);

    % H-field increment
    Hx_old = Hx_new;
    Hy_old = Hy_new;
    Hz_old = Hz_new;

    % E-field calculate
    ofs = offset;
    nfs = offset-1;
    ii = ofs:N_x-nfs;
    jj = ofs:N_y-nfs;
    kk = ofs:N_z-nfs;

    if superconducting_model == NONE
        Ex_new(ii,jj,kk) = C_a(ii,jj,kk).*Ex_old(ii,jj,kk) ...
            + C_b(ii,jj,kk).*(Hz_old(ii,jj,kk)-Hz_old(ii,jj-1,kk) ...
            + Hy_old(ii,jj,kk-1) - Hy_old(ii,jj,kk));
    
        Ey_new(ii,jj,kk) = C_a(ii,jj,kk).*Ey_old(ii,jj,kk) ...
            + C_b(ii,jj,kk).*(Hx_old(ii,jj,kk)-Hx_old(ii,jj,kk-1) ...
            + Hz_old(ii-1,jj,kk) - Hz_old(ii,jj,kk));
    
        Ez_new(ii,jj,kk) = C_a(ii,jj,kk).*Ez_old(ii,jj,kk) ...
            + C_b(ii,jj,kk).*(Hy_old(ii,jj,kk)-Hy_old(ii-1,jj,kk) ...
            + Hx_old(ii,jj-1,kk) - Hx_old(ii,jj,kk));

    elseif superconducting_model == TWOFLUID      
        %update E n+1
        Ex_new(ii,jj,kk) = C_a(ii,jj,kk).*Ex_old(ii,jj,kk) ...
            + C_b(ii,jj,kk).*(Hz_old(ii,jj,kk)-Hz_old(ii,jj-1,kk) ...
            + Hy_old(ii,jj,kk-1) - Hy_old(ii,jj,kk)) ...
            + C_c(ii,jj,kk).*0.5.*(alp_s+1).*J_sx_old(ii,jj,kk);
    
        Ey_new(ii,jj,kk) = C_a(ii,jj,kk).*Ey_old(ii,jj,kk) ...
            + C_b(ii,jj,kk).*(Hx_old(ii,jj,kk)-Hx_old(ii,jj,kk-1) ...
            + Hz_old(ii-1,jj,kk) - Hz_old(ii,jj,kk)) ...
            + C_c(ii,jj,kk).*0.5.*(alp_s+1).*J_sy_old(ii,jj,kk);
    
        Ez_new(ii,jj,kk) = C_a(ii,jj,kk).*Ez_old(ii,jj,kk) ...
            + C_b(ii,jj,kk).*(Hy_old(ii,jj,kk)-Hy_old(ii-1,jj,kk) ...
            + Hx_old(ii,jj-1,kk) - Hx_old(ii,jj,kk)) ...
            + C_c(ii,jj,kk).*0.5.*(alp_s+1).*J_sz_old(ii,jj,kk);
           
        %update J
       J_sx_new = alp_s*J_sx_old+gam_s*(Ex_new+Ex_old);
       J_sy_new = alp_s*J_sy_old+gam_s*(Ey_new+Ey_old);
       J_sz_new = alp_s*J_sz_old+gam_s*(Ez_new+Ez_old);


        J_sx_old = J_sx_new.*sc_mask;
        J_sy_old = J_sy_new.*sc_mask;
        J_sz_old = J_sz_new.*sc_mask;
        
%         J_nx_old = J_nx_new.*sc_mask;
%         J_ny_old = J_ny_new.*sc_mask;
%         J_nz_old = J_nz_new.*sc_mask;
    end

    %SFTF source insert
    Ey_new(sftf_x,jj,kk) = Ey_new(sftf_x,jj,kk) ...
        + C_b(sftf_x,jj,kk).*Hz_inc(1,jj,kk);
    Ez_new(sftf_x,jj,kk) = Ez_new(sftf_x,jj,kk) ...
        - C_b(sftf_x,jj,kk).*Hy_inc(1,jj,kk);
    
    % apply PEC: set E-fields to zero in PEC region
    Ex_new = Ex_new.*pec_mask;
    Ey_new = Ey_new.*pec_mask;
    Ez_new = Ez_new.*pec_mask;

    % E-field Boundary Conditions
    Ex_new = assist_mur_abc_plane(c_, delta_t, delta, N, bc_offset, ...
        Ex_new, Ex_old, Ex_old_old);
    Ey_new = assist_mur_abc_plane(c_, delta_t, delta, N, bc_offset, ...
        Ey_new, Ey_old, Ey_old_old);
    Ez_new = assist_mur_abc_plane(c_, delta_t, delta, N, bc_offset, ...
        Ez_new, Ez_old, Ez_old_old);

    Ex_new = assist_mur_abc_line(c_, delta_t, delta, N, bc_offset, cos_a, ...
        Ex_new, Ex_old);
    Ey_new = assist_mur_abc_line(c_, delta_t, delta, N, bc_offset, cos_a, ...
        Ey_new, Ey_old);
    Ez_new = assist_mur_abc_line(c_, delta_t, delta, N, bc_offset, cos_a, ...
        Ez_new, Ez_old);

    % E-field increment
    Ex_old_old = Ex_old;
    Ey_old_old = Ey_old;
    Ez_old_old = Ez_old;

    Ex_old = Ex_new;
    Ey_old = Ey_new;
    Ez_old = Ez_new;

    %store values to monitors
    for num = 1:num_monitors
        [iii,jjj,kkk] = unpack_coords(monitor(num).coords);

        if monitor(num).fields_to_monitor(1)
            monitor_values{num}{1}(:,:,step+1) = Ex_old(iii,jjj,kkk);
        end
        if monitor(num).fields_to_monitor(1)
            monitor_values{num}{2}(:,:,step+1) = Ey_old(iii,jjj,kkk);
        end
        if monitor(num).fields_to_monitor(3)
            monitor_values{num}{3}(:,:,step+1) = Ez_old(iii,jjj,kkk);
        end
        if monitor(num).fields_to_monitor(4)
            monitor_values{num}{4}(:,:,step+1) = Hx_old(iii,jjj,kkk);
        end
        if monitor(num).fields_to_monitor(5)
            monitor_values{num}{5}(:,:,step+1) = Hy_old(iii,jjj,kkk);
        end
        if monitor(num).fields_to_monitor(6)
            monitor_values{num}{6}(:,:,step+1) = Hz_old(iii,jjj,kkk);
        end

    end

    step = step+1;

    if step >= N_t_max || (use_bootstrapped_fields && step >= N_t_max-TEMP_BOOTSTRAP_OFFSET)
        stop_cond = true;
    end

end
fprintf('\nSaving monitors (this might take a while)...\n');

save('monitor.mat','monitor_values','-v7.3');
save('monitor.mat','monitor_names','source_val_E','delta_t','delta_z','-append');

if use_bootstrapped_fields
    source_from_bootstrap = source_field_Ez(TEMP_BOOTSTRAP_OFFSET:end,source_y,source_z);
    save('monitor.mat','source_from_bootstrap','-append');
end

fprintf('Duration: %.0f seconds.\n',toc(main_timer));
fprintf('Simulation end.\n');

%====================Function Declarations====================%

function text_update(step,N_t_max,delta_t)
   persistent reverseStr;
   persistent time_per_iteration;
   persistent time_estimate;

   if step == 0
       tic;
       reverseStr = '';
       time_per_iteration = 1;
   end
   time_per_iteration = time_per_iteration*0.7+0.3*toc;
   
   if mod(step,10) == 0 && step ~=0
       time_estimate = ceil((N_t_max-step)*time_per_iteration/60);
   end

   msg = sprintf('Step %d/%d  t = %s s\nEst. time remaining %d mins', step, N_t_max, num2eng(delta_t*step),time_estimate);
    fprintf([reverseStr, msg]);
   reverseStr = repmat(sprintf('\b'), 1, length(msg));
   tic;
end


function [ii_,jj_,kk_] = unpack_coords(coord)
    ii_ = coord{1};
    jj_ = coord{2};
    kk_ = coord{3};
end

function W_new_ = assist_mur_abc_plane(c_, delta_t, deltas, Ns, offset, W_new, W_old, W_old_old)

    [delta_x, delta_y, delta_z] = unpack_coords(deltas);
    [N_x, N_y, N_z] = unpack_coords(Ns);

    ofs = offset;
    nfs = offset-1;

    ii = ofs+1:N_x-nfs-1;
    jj = ofs+1:N_y-nfs-1;
    kk = ofs+1:N_z-nfs-1;

    W_new_ = W_new;

    W_new_(ofs,jj,kk) = mur_abc_plane('x', c_, delta_t, delta_x, ...
        ofs, N_x, jj, kk, W_new, W_old, W_old_old);
    W_new_(N_x-nfs,jj,kk) = mur_abc_plane('x', c_, delta_t, delta_x, ...
        N_x-nfs, N_x, jj, kk, W_new, W_old, W_old_old);

    W_new_(ii,ofs,kk) = mur_abc_plane('y', c_, delta_t, delta_y, ...
        ofs, N_y, kk, ii, W_new, W_old, W_old_old);
    W_new_(ii,N_y-nfs,kk) = mur_abc_plane('y', c_, delta_t, delta_y, ...
        N_y-nfs, N_y, kk, ii, W_new, W_old, W_old_old);

    W_new_(ii,jj,ofs) = mur_abc_plane('z', c_, delta_t, delta_z, ...
        ofs, N_z, ii, jj, W_new, W_old, W_old_old);
    W_new_(ii,jj,N_z-nfs) = mur_abc_plane('z', c_, delta_t, delta_z, ...
        N_z-nfs, N_z, ii, jj, W_new, W_old, W_old_old);
end


function W_new_ = mur_abc_plane(boundary, c, delta_t, delta, ii, N, jj, kk, ...
    W_new, W_old, W_old_old)
    
    i0 = ii;

    if ii < N/2
        i1 = i0+1;
    else
        i1 = i0-1;
    end
    
    iii = {i1, ...
        i1,i0, ...
        i0,i1};
    jjj = {jj, ...
        jj,jj, ...
        jj,jj};
    kkk = {kk, ...
        kk,kk, ...
        kk,kk};

    if boundary == 'x'
        %default
        
    elseif boundary == 'y'
        tempi = iii;
        tempj = jjj;
        tempk = kkk;
        iii = tempk;
        jjj = tempi;
        kkk = tempj;
    elseif boundary == 'z'
        tempi = iii;
        tempj = jjj;
        tempk = kkk;
        iii = tempj;
        jjj = tempk;
        kkk = tempi;
    end
    
    c_ = c(iii{1},jjj{1},kkk{1});

    coeffs_1 = -1;
    coeffs_2 = (c_*delta_t-delta)./(c_*delta_t+delta);
    coeffs_3 = (2*delta)./(c_*delta_t+delta);

    W_new_ = coeffs_1.* W_old_old(iii{1},jjj{1},kkk{1}) ...
        +coeffs_2.* (W_new(iii{2},jjj{2},kkk{2})+W_old_old(iii{3},jjj{3},kkk{3})) ...
        +coeffs_3.* (W_old(iii{4},jjj{4},kkk{4})+W_old(iii{5},jjj{5},kkk{5}));
end

function W_new_ = assist_mur_abc_line(c, delta_t, deltas, Ns, offset, cos_a, ...
    W_new, W_old)

    [delta_x, ~, ~] = unpack_coords(deltas);
    [N_x, N_y, N_z] = unpack_coords(Ns);

    ofs = offset;
    nfs = offset-1;

    ii = ofs+1:N_x-nfs-1;
    jj = ofs+1:N_y-nfs-1;
    kk = ofs+1:N_z-nfs-1;
    
    W_new_ = W_new;

    i0 = {ofs, N_x-nfs,ii};
    j0 = {ofs, N_y-nfs,jj};
    k0 = {ofs, N_z-nfs,kk};

    i1 = {ofs+1, N_x-nfs-1,ii};
    j1 = {ofs+1, N_y-nfs-1,jj};
    k1 = {ofs+1, N_z-nfs-1,kk};
    
    

    for ind = 1:3
        for jnd = 1:3
            for knd = 1:3
                if sum([ind==3,jnd==3,knd==3]) == 1
                    coords_0 = {i0{ind}, j0{jnd}, k0{knd}};
                    coords_1 = {i1{ind}, j1{jnd}, k1{knd}};
                    W_new_(i0{ind}, j0{jnd}, k0{knd}) = ...
                    mur_abc_point(c, delta_t, delta_x, cos_a, coords_0, coords_1,...
                    W_new, W_old);   
                end          
            end
        end
    end
end


function W_new_ = mur_abc_point(c, delta_t, delta_x, cos_a, coords_0, coords_1, ...
    W_new, W_old)
    
    [ii_0,jj_0,kk_0] = unpack_coords(coords_0);
    [ii_1,jj_1,kk_1] = unpack_coords(coords_1);

    c_ = c(ii_0,jj_0,kk_0);
        
    W_new_ = W_old(ii_1,jj_1,kk_1) ...
        + (c_*delta_t*cos_a-delta_x)./(c_*delta_t*cos_a+delta_x) ...
        .*(W_new(ii_1,jj_1,kk_1) - W_old(ii_0,jj_0,kk_0));

end

function plot_line(field1,index,step,y_text,fig_no,y_max)
    figure(fig_no)
    field1 = reshape(field1,1,length(field1));
    plot(index,[field1]);

    name = sprintf('n = %d',step);
    title(name);
    ylabel(y_text);
    xlabel('x (m)');
    grid on   

if exist('y_max','var')
     ylim([0 y_max]);
end

end

function plot_field(field,slice_x,slice_y,slice_z,step,deltas,delta_t,c_max)

figure(3)
    colormap jet
    s = size(field);
    field(s(1)+1,s(2)+1,s(3)+1) = 0;
    
    field = permute(field,[2 1 3]);
    
    if ~exist('deltas','var')

        h = slice(field,slice_x,slice_y,slice_z);
        set(h,'edgecolor','none')
    
        name = sprintf('n = %d',step);
        title(name);
        xlabel('x (cells)');
        ylabel('y (cells)');
        zlabel('z (cells)');
        lim = max(s)+1;
        xlim([0 lim])
        ylim([0 lim])
        zlim([0 lim])

    else
        [delta_x,delta_y,delta_z] = unpack_coords(deltas);

        x_index = (0:s(1)).*delta_x;
        y_index = (0:s(2)).*delta_y;
        z_index = (0:s(3)).*delta_z;
        
        [plot_x,plot_y,plot_z] = meshgrid(x_index,y_index,z_index);
        h = slice(plot_x,plot_y,plot_z,field, ...
            slice_x*delta_x, ...
            slice_y*delta_y, ...
            slice_z*delta_z);
        set(h,'edgecolor','none')
    
        name = sprintf('n = %d, t = %s (s)',step,num2eng(step*delta_t));
        title(name);
        xlabel('x (m)');
        ylabel('y (m)');
        zlabel('z (m)');
        xlim([0 max(x_index)])
        ylim([0 max(y_index)])
        zlim([0 max(z_index)])
    end

    bar = colorbar();
    ylabel(bar,'|H_{tot}| (A/m)');
    grid on   
if exist('c_max','var')
     clim([0 c_max]);
end

end

function print_text_at_sim_start(vol,delta_x,delta_t)
    fprintf('Simulation start:\n  #Cells = %d\n  delta_x = %s m\n  delta_t = %s s\n'...
    ,vol,num2eng(delta_x),num2eng(delta_t))

    fprintf(strcat(datestr(datetime( ...
    'now','TimeZone','local','Format','d-MM-y HH:mm:ss')),'\n'));
end