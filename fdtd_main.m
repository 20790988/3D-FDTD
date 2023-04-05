%
%  3D-FDTD
%  K Meyer
%  2023-03-21
% 

clear

% Material specification
    [param, material, source] = import_model;

    sigma = param.material(1,:);
    sigma_m = param.material(2,:);
    
    epsilon_0 = 8.8542e-12;
    epsilon_r = param.material(3,:);
    epsilon = epsilon_0*epsilon_r;
    
    mu_0 = 1.2566e-6;
    mu_r = param.material(4,:);
    mu = mu_0*mu_r;

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
    delta_t = delta_x/(max(c)*sqrt(3));

% Simulation length
    
    N_t_max = floor(param.M_t_max/delta_t);

% Model
%     material = import_model(N_x,N_y,N_z,delta_x,delta_x,delta_x);

% Space around E-field
    offset = 3;
    
%====================SIMULATION SETUP END=====================%

fprintf('Setup start\n')



% Excitation (J)
source_coords = source.coord;
source_x = source_coords{1};
source_y = source_coords{2};
source_z = source_coords{3};


%simulation stability check
S = max(c)*delta_t/delta_x;

if S>1
    fprintf('Simulation unstable');
    return
end

% setup simulation space
C_a_single = (1-(sigma.*delta_t)./(2.*epsilon))./(1+(sigma.*delta_t)./(2.*epsilon));
C_b_single = (delta_t./(epsilon.*delta_x))./(1+(sigma.*delta_t)./(2.*epsilon));

D_a_single = (1-(sigma_m.*delta_t)./(2.*mu))./(1+(sigma_m.*delta_t)./(2.*mu));
D_b_single = (delta_t./(mu.*delta_x))./(1+(sigma_m.*delta_t)./(2.*mu));

C_a = C_a_single(material);
C_b = C_b_single(material);

D_a = D_a_single(material);
D_b = D_b_single(material);

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

Jsource_x = zeros(N_x,N_y,N_z);
Jsource_y = zeros(N_x,N_y,N_z);
Jsource_z = zeros(N_x,N_y,N_z);

source_x = cast(source_x,'int32');
source_y = cast(source_y,'int32');
source_z = cast(source_z,'int32');

step = 0;
stop_cond = false;


fprintf('Simulation start:\n  N_x = %d\n  delta_x = %.4e m\n  delta_t = %.4e s\n  ',N_x,delta_x,delta_t)
fprintf(strcat(datestr(datetime('now','TimeZone','local','Format','d-MM-y HH:mm:ss')),'\n'));
%====================LOOP START=====================%

while stop_cond == false

    text_update(step,N_t_max,delta_t)

%     Ex_old(source_x,source_y,source_z) = source_signal(step+1);
%     Ey_old(source_x,source_y,source_z) = source_signal(step+1);
%     Jsource_x(source_x,source_y,source_z) = source_signal(step+1);
%     Jsource_y(source_x,source_y,source_z) = source_signal(step+1);
    Jsource_z(source_x,source_y,source_z) = source.value{3}(step*delta_t);

    %====================PLOTTING START=====================%
    if step>0
        
        H_tot = sqrt(Hx_old.^2+Hy_old.^2+Hz_old.^2);
        plot_field(H_tot,N_x/2,N_y/2,N_z/2,step,delta,delta_t);
        
%         H_tot_line = (H_tot(:,N_y/2,N_z/2));
%         plot_line(H_tot_line,delta_x*(1:N_x),step);
        
%         E_tot = sqrt(Ex_old.^2+Ey_old.^2+Ez_old.^2);
%          plot_field(E_tot,N_x/2,[],N_z/2,step);
%         
%         E_tot_line = (E_tot(:,N_y/2,N_z/2));
%         plot_line(E_tot_line,delta_x*(1:N_x),step);
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


    Ex_new(ii,jj,kk) = C_a(ii,jj,kk).*Ex_old(ii,jj,kk) ...
        + C_b(ii,jj,kk).*(Hz_old(ii,jj,kk)-Hz_old(ii,jj-1,kk) ...
        + Hy_old(ii,jj,kk-1) - Hy_old(ii,jj,kk) ...
        + Jsource_x(ii,jj,kk).*delta_x);

    Ey_new(ii,jj,kk) = C_a(ii,jj,kk).*Ey_old(ii,jj,kk) ...
        + C_b(ii,jj,kk).*(Hx_old(ii,jj,kk)-Hx_old(ii,jj,kk-1) ...
        + Hz_old(ii-1,jj,kk) - Hz_old(ii,jj,kk) ...
        + Jsource_y(ii,jj,kk).*delta_x);

    Ez_new(ii,jj,kk) = C_a(ii,jj,kk).*Ez_old(ii,jj,kk) ...
        + C_b(ii,jj,kk).*(Hy_old(ii,jj,kk)-Hy_old(ii-1,jj,kk) ...
        + Hx_old(ii,jj-1,kk) - Hx_old(ii,jj,kk) ...
        + Jsource_z(ii,jj,kk).*delta_x);

    % E-field Boundary Conditions
    c_ = c(border_material_index);   
    
    temp_offset = 2;

    %planes
    

    Ex_new = assist_mur_abc_plane(c_, delta_t, delta, N, temp_offset, Ex_new, Ex_old, Ex_old_old);
    Ey_new = assist_mur_abc_plane(c_, delta_t, delta, N, temp_offset, Ey_new, Ey_old, Ey_old_old);
    Ez_new = assist_mur_abc_plane(c_, delta_t, delta, N, temp_offset, Ez_new, Ez_old, Ez_old_old);

    %lines
   
    Ex_new = assist_mur_abc_line(c_, delta_t, delta, N, temp_offset, pi/4, Ex_new, Ex_old);
    Ey_new = assist_mur_abc_line(c_, delta_t, delta, N, temp_offset, pi/4, Ey_new, Ey_old);
    Ez_new = assist_mur_abc_line(c_, delta_t, delta, N, temp_offset, pi/4, Ez_new, Ez_old);

% 
%     %points
%     Ex_new = assist_mur_abc_point(c_, delta_t, delta_x, N, offset, pi/4, Ex_new, Ex_old);
%     Ey_new = assist_mur_abc_point(c_, delta_t, delta_x, N, offset, pi/4, Ey_new, Ey_old);
%     Ez_new = assist_mur_abc_point(c_, delta_t, delta_x, N, offset, pi/4, Ez_new, Ez_old);

    % E-field increment
    Ex_old_old = Ex_old;
    Ey_old_old = Ey_old;
    Ez_old_old = Ez_old;

    Ex_old = Ex_new;
    Ey_old = Ey_new;
    Ez_old = Ez_new;
    
    step = step+1;
    if step >= N_t_max
        stop_cond = true;
    end

end

fprintf('Simulation end.\n');

%====================Function Declarations====================%

function text_update(step,N_t_max,delta_t)
   persistent reverseStr;
   if step == 0 
       reverseStr = '';
   end

   msg = sprintf('Step %d/%d  t = %.3e s\n', step, N_t_max, delta_t*step);
   fprintf([reverseStr, msg]);
   reverseStr = repmat(sprintf('\b'), 1, length(msg));
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

function W_new_ = mur_abc_plane(boundary, c, delta_t, delta, ii, N, jj, kk, W_new, W_old, W_old_old)
    
    i0 = ii;

    if ii < N/2
        i1 = i0+1;
    else
        i1 = i0-1;
    end
    
%     fprintf('%c %d %d\n',boundary,ii,i1);

    iii = {i1, ...
        i1,i0, ...
        i0,i1, ...
        i0,i0,i0,i1,i1,i1,i0,i0,i1,i1}; 
    jjj = {jj, ...
        jj,jj, ...
        jj,jj, ...
        jj+1,jj,jj-1,jj+1,jj,jj-1,jj,jj,jj,jj}; 
    kkk = {kk, ...
        kk,kk, ...
        kk,kk, ...
        kk,kk,kk,kk,kk,kk,kk+1,kk-1,kk+1,kk-1}; 

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
       
    coeffs = [-1 ...
        (c*delta_t-delta)/(c*delta_t+delta) ...
        (2*delta)/(c*delta_t+delta) ...
        (c*delta_t).^2/(2*delta*(c*delta_t+delta))];

    W_new_ = coeffs(1)* W_old_old(iii{1},jjj{1},kkk{1}) ...
        +coeffs(2)* (W_new(iii{2},jjj{2},kkk{2})+W_old_old(iii{3},jjj{3},kkk{3})) ...
        +coeffs(3)* (W_old(iii{4},jjj{4},kkk{4})+W_old(iii{5},jjj{5},kkk{5})) ...
        +coeffs(4)* (W_old(iii{6},jjj{6},kkk{6}) -4*W_old(iii{7},jjj{7},kkk{7}) +W_old(iii{8},jjj{8},kkk{8}) ...
        +W_old(iii{9},jjj{9},kkk{9}) -4*W_old(iii{10},jjj{10},kkk{10}) +W_old(iii{11},jjj{11},kkk{11})...
        +W_old(iii{12},jjj{12},kkk{12}) +W_old(iii{13},jjj{13},kkk{13}) +W_old(iii{14},jjj{14},kkk{14}) +W_old(iii{15},jjj{15},kkk{15}));
end

function W_new_ = mur_abc_plane_with_permute(boundary, c, delta_t, delta, ii, N, jj, kk, W_new, W_old, W_old_old)
    
    i0 = ii;

    if ii < N/2
        i1 = i0+1;
    else
        i1 = i0-1;
    end

    if boundary == 'x'
        %default
    elseif boundary == 'y'
        W_new = permute(W_new,[2 1 3]);
        W_old = permute(W_old,[2 1 3]);
        W_old_old = permute(W_old_old,[2 1 3]);
    elseif boundary == 'z'
        W_new = permute(W_new,[3 2 1]);
        W_old = permute(W_old,[3 2 1]);
        W_old_old = permute(W_old_old,[3 2 1]);
    end
    
    coeffs = [-1 ...
        (c*delta_t-delta)/(c*delta_t+delta) ...
        (2*delta)/(c*delta_t+delta) ...
        (c*delta_t).^2/(2*delta*(c*delta_t+delta))];

    W_new_ = coeffs(1)* W_old_old(i1,jj,kk) ...
        +coeffs(2)* (W_new(i1,jj,kk)+W_old_old(i0,jj,kk)) ...
        +coeffs(3)* (W_old(i0,jj,kk)+W_old(i1,jj,kk)) ...
        +coeffs(4)* (W_old(i0,jj+1,kk) -4*W_old(i0,jj,kk) +W_old(i0,jj-1,kk) ...
        +W_old(i1,jj+1,kk) -4*W_old(i1,jj,kk) +W_old(i1,jj-1,kk)...
        +W_old(i0,jj,kk+1) +W_old(i0,jj,kk-1) +W_old(i1,jj,kk+1) +W_old(i1,jj,kk-1));

    if boundary == 'x'
        %default
    elseif boundary == 'y'
        W_new_ = permute(W_new_,[2 1 3]);       
    elseif boundary == 'z'
        W_new_ = permute(W_new_,[3 2 1]);
    end
end

function W_new_ = assist_mur_abc_line(c, delta_t, deltas, Ns, offset, a, W_new, W_old)

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

    for ind = [3:-1:1]
        for jnd = [3:-1:1]
            for knd = [3:-1:1]
                if sum([ind==3,jnd==3,knd==3]) == 1
                    coords_0 = {i0{ind}, j0{jnd}, k0{knd}};
                    coords_1 = {i1{ind}, j1{jnd}, k1{knd}};
                    W_new_(i0{ind}, j0{jnd}, k0{knd}) = ...
                    mur_abc_point(c, delta_t, delta_x, a, coords_0, coords_1, W_new, W_old);   
                end          
            end
        end
    end
end

function W_new_ = assist_mur_abc_point(c, delta_t, delta_x, Ns, offset, a, W_new, W_old)
    
    [N_x, N_y, N_z] = unpack_coords(Ns);

    ofs = offset;
    nfs = offset-1;
    
    W_new_ = W_new;

    i0 = {ofs, N_x-nfs};
    j0 = {ofs, N_y-nfs};
    k0 = {ofs, N_z-nfs};

    i1 = {ofs+1, N_x-nfs-1};
    j1 = {ofs+1, N_y-nfs-1};
    k1 = {ofs+1, N_z-nfs-1};

    for ind = [1,2]
        for jnd = [1,2]
            for knd = [1,2]
                coords_0 = {i0{ind}, j0{jnd}, k0{knd}};
                coords_1 = {i1{ind}, j1{jnd}, k1{knd}};
                W_new_(i0{ind}, j0{jnd}, k0{knd}) = ...
                mur_abc_point(c, delta_t, delta_x, a, coords_0, coords_1, W_new, W_old);             
            end
        end
    end
end

function W_new_ = mur_abc_point(c, delta_t, delta_x, a, coords_0, coords_1, W_new, W_old)

    [ii_0,jj_0,kk_0] = unpack_coords(coords_0);
    [ii_1,jj_1,kk_1] = unpack_coords(coords_1);

    W_new_ = W_old(ii_1,jj_1,kk_1) ...
        + (c*delta_t*cos(a)-delta_x)/(c*delta_t*cos(a)+delta_x) ...
        *(W_new(ii_1,jj_1,kk_1) - W_old(ii_0,jj_0,kk_0));

    end

function plot_line(field1,index,step)
    figure(2)
    field1 = reshape(field1,1,length(field1));
    plot(index,[field1]);
%     theory_field = 1./(2*pi*(abs(index-5).^2));
%     fit_field = 0.01596./((abs(index-5)+2.5));
    name = sprintf('n = %d',step);
    title(name);
    ylabel('|H_{tot}| (A/m)');
    xlabel('x (m)');
    grid on   
%     ylim([0 7e-5])
    ylim([0 0.025]);
end

function plot_field_slice(field,step)
    
    figure(1)
    colormap jet
    s = size(field);
    field(s(1)+1,s(2)+1) = 0;
    h = surf(field);

%     view([0 0 1])
    set(h,'edgecolor','none')

    name = sprintf('n = %d',step);
    title(name);
    xlabel('x');
    ylabel('y');
   
%     set(gca,'ColorScale','log')

    bar = colorbar();
    ylabel(bar,'|H_{tot}| (A/m)');
        
    grid on   
    zlim([0 1]);
    clim([0 1]);
end

function plot_field(field,slice_x,slice_y,slice_z,step,deltas,delta_t)

figure(1)
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
        h = slice(plot_x,plot_y,plot_z,field,slice_x*delta_x,slice_y*delta_y,slice_z*delta_z);
        set(h,'edgecolor','none')
    
        name = sprintf('n = %d, t = %s (s)',step,num2eng(step*delta_t));
        title(name);
        xlabel('x (m)');
        ylabel('y (y)');
        zlabel('z (z)');
        xlim([0 max(x_index)])
        ylim([0 max(y_index)])
        zlim([0 max(z_index)])
    end

    bar = colorbar();
    ylabel(bar,'|H_{tot}| (A/m)');
    grid on   
    clim([0 2e-5]);

end