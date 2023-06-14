%
%  3D-FDTD
%  K Meyer
%  2023-04-20
% 

clear

fprintf('Start\n')
%====================MODEL IMPORT=====================%

    [param, material, source, monitor] = scond_test3;
    
    sigma = param.material(1,:);
    sigma_m = param.material(2,:);
    
    epsilon_0 = 8.8542e-12;
    epsilon_r = param.material(3,:);
    epsilon = epsilon_0*epsilon_r;
    
    mu_0 = 1.2566e-6;
    mu_r = param.material(4,:);
    mu = mu_0*mu_r;
    
    superconducting_model = 1;
    NONE = 0;
    TWOFLUID = 1;

    lambda_L = param.lambda_L;
    sigma_n = param.sigma_n;
    T_op = param.T_op;
    T_c = param.T_c;

% Material at Border
    border_material_index = param.border_material_index;

% Grid and cell size
    N = param.N;

    N_x = N{1};
    N_y = N{2};
    N_z = N{3};
    
    delta = param.delta;
    delta_x = delta{1};
    delta_y = delta{2};
    delta_z = delta{3};
    
    c = 1./sqrt(epsilon.*mu);
    delta_t = delta_x/(2*max(c)*sqrt(3));
    
% Simulation length
    N_t_max = floor(param.M_t_max/delta_t);

% Fix sigma
    sigma = sqrt(sigma./(2*pi*10e9*mu))./delta_x;

% Space around E-field in grid cells
    offset = 3;

% Source plane
    source_coords = source.coord;
    sftf_x = source_coords{1};
    source_y = source_coords{2};
    source_z = source_coords{3};

%     bootstrap = load("bootstrap.mat");
%     bs_E = bootstrap.tempE;
%     bs_E = bs_E./bs_E(1,N_y/2-0.5,N_z/2-0.5);
% 
%     bs_H = bootstrap.tempH;
%     bs_H = bs_H./bs_H(1,N_y/2-0.5,N_z/2-0.5);


%     source_val_x = source.value{1}((0:N_t_max)*delta_t);
%     source_val_y = source.value{2}((0:N_t_max)*delta_t);

    t = (0:N_t_max)*delta_t;
    [source_val_E,source_val_H] = source.value{3}(t,delta_t,delta_x);

    figure(1);
    plot(t,source_val_E),
    clf;

    source_N_t_max = floor(source.t_max/delta_t);

%====================SIMULATION SETUP AND INITIALISE=====================%

S = max(c)*delta_t/delta_x;

if S>1
    fprintf('Simulation unstable');
    return
end

%superconducting constants
lambda_L = lambda_L/(sqrt(1-(T_op/T_c)^4));

w_n2 = ((T_op/T_c)^4)/(lambda_L^2*mu_0*epsilon_0);
w_s2 = (1-(T_op/T_c)^4)/(lambda_L^2*mu_0*epsilon_0);


tau_n = sigma_n*lambda_L^2*mu_0;
gamma_n = 1/tau_n;

alp_s = 2;
xii_s = -1;
gam_s = (delta_t)^2*epsilon_0*w_s2;

alp_n = 4/(delta_t/tau_n+2);
xii_n = (delta_t/tau_n-2)/(delta_t/tau_n+2);
gam_n = (epsilon_0*w_n2*2*(delta_t^2))/(delta_t/tau_n+2);


%FDTD parameters
if superconducting_model == NONE

    C_a_single = (1-(sigma.*delta_t)./(2.*epsilon)) ...
        ./(1+(sigma.*delta_t)./(2.*epsilon));
    C_b_single = (delta_t./(epsilon.*delta_x)) ...
        ./(1+(sigma.*delta_t)./(2.*epsilon));
    C_c_single = 0;

elseif superconducting_model == TWOFLUID
    gam_s_ = [0 0 gam_s];
    gam_n_ = [0 0 gam_n];

    half_sum_gam = 0.5*(gam_s_+gam_n_);

    C_a_single = half_sum_gam./(2*epsilon+half_sum_gam+sigma*delta_t);
    C_b_single = (2*epsilon-sigma*delta_t) ...
        ./(2*epsilon+half_sum_gam+sigma*delta_t);
    C_c_single = (2*delta_t)./(2*epsilon+half_sum_gam+sigma*delta_t);
   
end

D_a_single = (1-(sigma_m.*delta_t)./(2.*mu))./(1+(sigma_m.*delta_t)./(2.*mu));
D_b_single = (delta_t./(mu.*delta_x))./(1+(sigma_m.*delta_t)./(2.*mu));

pec_mask = (material ~= 0);
material(material==0) = 1;

sc_mask = (material == 3);

C_a = C_a_single(material);
C_b = C_b_single(material);
C_c = C_c_single(material);

D_a = D_a_single(material);
D_b = D_b_single(material);

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
source_y = floor(source_y);
source_z = floor(source_z);

Hz_inc = zeros(1,N_y,N_z);
Hy_inc = zeros(1,N_y,N_z);
Ey_inc = zeros(1,N_y,N_z);
Ez_inc = zeros(1,N_y,N_z);
  
%superconducting
J_sx_new = zeros(N_x,N_y,N_z);
J_sy_new = zeros(N_x,N_y,N_z);
J_sz_new = zeros(N_x,N_y,N_z);

J_nx_new = zeros(N_x,N_y,N_z);
J_ny_new = zeros(N_x,N_y,N_z);
J_nz_new = zeros(N_x,N_y,N_z);

J_sx_old = zeros(N_x,N_y,N_z);
J_sy_old = zeros(N_x,N_y,N_z);
J_sz_old = zeros(N_x,N_y,N_z);

J_nx_old = zeros(N_x,N_y,N_z);
J_ny_old = zeros(N_x,N_y,N_z);
J_nz_old = zeros(N_x,N_y,N_z);

J_sx_old_old = zeros(N_x,N_y,N_z);
J_sy_old_old = zeros(N_x,N_y,N_z);
J_sz_old_old = zeros(N_x,N_y,N_z);

J_nx_old_old = zeros(N_x,N_y,N_z);
J_ny_old_old = zeros(N_x,N_y,N_z);
J_nz_old_old = zeros(N_x,N_y,N_z);

%monitors
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
    monitor_values{num} = zeros(rows,cols,N_t_max);
    monitor_names{num} = monitor(num).name;
end

step = 0;
stop_cond = false;

fprintf('Simulation start:\n  #Cells = %d\n  delta_x = %s m\n  delta_t = %s s\n'...
    ,N_x*N_y*N_z,num2eng(delta_x),num2eng(delta_t))
fprintf(strcat(datestr(datetime( ...
    'now','TimeZone','local','Format','d-MM-y HH:mm:ss')),'\n'));

%====================LOOP START=====================%

while stop_cond == false
    text_update(step,N_t_max,delta_t)

   
    Ez_inc(1,source_y,source_z) = (source_val_E(step+1));
    Hy_inc(1,source_y,source_z) = (source_val_H(step+1));

  
    %====================PLOTTING START=====================%
    if  step>0
        
        H_tot = sqrt(Hx_old.^2+Hy_old.^2+Hz_old.^2);
%          plot_field(H_tot,N_x/2,N_y/2,14,step,delta,delta_t);

        tempz = floor(N_z/2);
        tempy = floor(N_y/2);
        H_tot_line = (H_tot(:,tempy,tempz));
%         plot_line(H_tot_line,delta_x*(0:N_x-1),step,'|H_{tot}| (A/m)',1,1/eta);
        
        E_tot = sqrt(Ex_old.^2+Ey_old.^2+Ez_old.^2);
        plot_field(E_tot,N_x/2,N_y/2,N_z/2,step,delta,delta_t,1);
%         view([0 0 1])

        tempz = floor(N_z/2);
        tempy = floor(N_y/2);
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

        Ex_new(ii,jj,kk) = C_a(ii,jj,kk).*Ex_old_old(ii,jj,kk) ...
            + C_b(ii,jj,kk).*Ex_old(ii,jj,kk) ...
            + C_c(ii,jj,kk).*((1/delta_x).*(Hz_old(ii,jj,kk)-Hz_old(ii,jj-1,kk) ...
            + Hy_old(ii,jj,kk-1) - Hy_old(ii,jj,kk)) ...
            + 0.5*((1+alp_s).*J_sx_old(ii,jj,kk)+xii_s.*J_sx_old_old(ii,jj,kk) ...
            + (1+alp_n).*J_nx_old(ii,jj,kk)+xii_n.*J_nx_old_old(ii,jj,kk)));
                
        Ey_new(ii,jj,kk) = C_a(ii,jj,kk).*Ey_old_old(ii,jj,kk) ...
            + C_b(ii,jj,kk).*Ey_old(ii,jj,kk) ...
            + C_c(ii,jj,kk).*((1/delta_x).*(Hx_old(ii,jj,kk)-Hx_old(ii,jj,kk-1) ...
            + Hz_old(ii-1,jj,kk) - Hz_old(ii,jj,kk)) ...
            + 0.5*((1+alp_s).*J_sy_old(ii,jj,kk)+xii_s.*J_sy_old_old(ii,jj,kk) ...
            + (1+alp_n).*J_ny_old(ii,jj,kk)+xii_n.*J_ny_old_old(ii,jj,kk)));
    
        Ez_new(ii,jj,kk) = C_a(ii,jj,kk).*Ez_old_old(ii,jj,kk) ...
            + C_b(ii,jj,kk).*Ez_old(ii,jj,kk) ...
            + C_c(ii,jj,kk).*((1/delta_x).*(Hy_old(ii,jj,kk)-Hy_old(ii-1,jj,kk) ...
            + Hx_old(ii,jj-1,kk) - Hx_old(ii,jj,kk)) ...
            + 0.5*((1+alp_s).*J_sz_old(ii,jj,kk)+xii_s.*J_sz_old_old(ii,jj,kk) ...
            + (1+alp_n).*J_nz_old(ii,jj,kk)+xii_n.*J_nz_old_old(ii,jj,kk)));
    
        %update J
        
        J_sx_new = update_J(J_sx_old,J_sx_old_old,Ex_new,Ex_old_old, ...
            alp_s,xii_s,gam_s,delta_t);
        J_sy_new = update_J(J_sy_old,J_sy_old_old,Ey_new,Ey_old_old, ...
            alp_s,xii_s,gam_s,delta_t);
        J_sz_new = update_J(J_sz_old,J_sz_old_old,Ez_new,Ez_old_old, ...
            alp_s,xii_s,gam_s,delta_t);
        
        J_nx_new = update_J(J_nx_old,J_nx_old_old,Ex_new,Ex_old_old, ...
            alp_n,xii_n,gam_n,delta_t);
        J_ny_new = update_J(J_ny_old,J_ny_old_old,Ey_new,Ey_old_old, ...
            alp_n,xii_n,gam_n,delta_t);
        J_nz_new = update_J(J_nz_old,J_nz_old_old,Ez_new,Ez_old_old, ...
            alp_n,xii_n,gam_n,delta_t);

        J_sx_old_old = J_sx_old;
        J_sy_old_old = J_sy_old;
        J_sz_old_old = J_sz_old;
        
        J_nx_old_old = J_nx_old;
        J_ny_old_old = J_ny_old;
        J_nz_old_old = J_nz_old;

        J_sx_old = J_sx_new.*sc_mask;
        J_sy_old = J_sy_new.*sc_mask;
        J_sz_old = J_sz_new.*sc_mask;
        
        J_nx_old = J_nx_new.*sc_mask;
        J_ny_old = J_ny_new.*sc_mask;
        J_nz_old = J_nz_new.*sc_mask;
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
    c_ = c(material);   
    
    bc_offset = 2;

    Ex_new = assist_mur_abc_plane(c_, delta_t, delta, N, bc_offset, ...
        Ex_new, Ex_old, Ex_old_old);
    Ey_new = assist_mur_abc_plane(c_, delta_t, delta, N, bc_offset, ...
        Ey_new, Ey_old, Ey_old_old);
    Ez_new = assist_mur_abc_plane(c_, delta_t, delta, N, bc_offset, ...
        Ez_new, Ez_old, Ez_old_old);

    Ex_new = assist_mur_abc_line(c_, delta_t, delta, N, bc_offset, pi/4, ...
        Ex_new, Ex_old);
    Ey_new = assist_mur_abc_line(c_, delta_t, delta, N, bc_offset, pi/4, ...
        Ey_new, Ey_old);
    Ez_new = assist_mur_abc_line(c_, delta_t, delta, N, bc_offset, pi/4, ...
        Ez_new, Ez_old);


    % E-field increment
    Ex_old_old = Ex_old;
    Ey_old_old = Ey_old;
    Ez_old_old = Ez_old;

    Ex_old = Ex_new;
    Ey_old = Ey_new;
    Ez_old = Ez_new;
    
%     E_tot = sqrt(Ex_old.^2+Ey_old.^2+Ez_old.^2);

    %store values to monitors
    for num = 1:num_monitors
%         num
        [iii,jjj,kkk] = unpack_coords(monitor(num).coords);
%         iii
%         jjj
%         kkk
%         monitor_values{num}(:,:,step+1)
        monitor_values{num}(:,:,step+1) = Ez_old(iii,jjj,kkk);
%         monitor_values{num}(:,:,step+1)
    end

    step = step+1;
    if step >= N_t_max
        stop_cond = true;
    end

end
fprintf('Saving monitors...\n');

save('monitor.mat','monitor_values','monitor_names','delta_t','delta_z');

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

function W_new_ = assist_mur_abc_line(c, delta_t, deltas, Ns, offset, a, ...
    W_new, W_old)

    [delta_x, delta_y, delta_z] = unpack_coords(deltas);
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

    for ind = 3:-1:1
        for jnd = 3:-1:1
            for knd = 3:-1:1
                if sum([ind==3,jnd==3,knd==3]) == 1
                    coords_0 = {i0{ind}, j0{jnd}, k0{knd}};
                    coords_1 = {i1{ind}, j1{jnd}, k1{knd}};
                    W_new_(i0{ind}, j0{jnd}, k0{knd}) = ...
                    mur_abc_point(c, delta_t, delta_x, a, coords_0, coords_1,...
                    W_new, W_old);   
                end          
            end
        end
    end
end


function W_new_ = mur_abc_point(c, delta_t, delta_x, a, coords_0, coords_1, ...
    W_new, W_old)
   
    [ii_0,jj_0,kk_0] = unpack_coords(coords_0);
    [ii_1,jj_1,kk_1] = unpack_coords(coords_1);

    c_ = c(ii_0,jj_0,kk_0);
        
    W_new_ = W_old(ii_1,jj_1,kk_1) ...
        + (c_*delta_t*cos(a)-delta_x)./(c_*delta_t*cos(a)+delta_x) ...
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

function J_new = update_J(J_old,J_old_old,E_new,E_old_old,alp,xii,gam,delta_t)
    J_new =alp*J_old+xii*J_old_old+(0.5*gam/delta_t)*(E_new-E_old_old);
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