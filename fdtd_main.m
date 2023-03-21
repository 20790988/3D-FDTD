%
%  3D-FDTD
%  K Meyer
%  2023-03-21
% 

clear

%====================SIMULATION SETUP START====================%

% Material specification
    sigma = [0 60e6 0];
    sigma_m = [0 0 0];
    
    epsilon_0 = 8.8542e-12;
    epsilon_r = [1 60e6 3];
    epsilon = epsilon_0*epsilon_r;
    
    mu_0 = 1.2566e-6;
    mu_r = [1 1 1];
    mu = mu_0*mu_r;

% Grid and cell size
    N_x = 50;
    N_y = N_x;
    N_z = N_x;
    
    delta_x = 3e-3;
    %delta_y = delta_x;
    %delta_z = delta_x;

% Simulation length
    N_t_max = 200;

% Model
    material = import_model(N_x,N_y,N_z,delta_x,delta_x,delta_x);

% Simulation instability max value
    e_field_max = 1e6;

% Space around E-field
    offset = 3;
    
%====================SIMULATION SETUP END=====================%

c = 1./sqrt(epsilon.*mu);
delta_t = delta_x/(max(c)*sqrt(3));

% Excitation (J)
[source_x, source_y, source_z, source_signal] = import_source(N_x,N_y,N_z,N_t_max,delta_t);

%simulation stability check
S = max(c)*delta_t/delta_x;
fprintf('S=%.4f\n',S);

if S>1
    fprintf('Simulation unstable');
    return
end

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

fprintf('simulation start\n')

while stop_cond == false

    if mod(step,50) == 0 || step<=5
        fprintf('step %d\n',step)
    end

    if max(Ez_old,[],'all')>e_field_max
        fprintf('Simulation unstable')
        return
    end
%     Jsource_x(source_x,source_y,source_z) = source_signal(step+1);
%     Jsource_y(source_x,source_y,source_z) = source_signal(step+1);
    Ex_old(source_x,source_y,source_z) = source_signal(step+1);
    Ey_old(source_x,source_y,source_z) = source_signal(step+1);
    Ez_old(source_x,source_y,source_z) = source_signal(step+1);

    if step>0
        E_tot = sqrt(Ex_old.^2+Ey_old.^2+Ez_old.^2);
        plot_field(E_tot,N_x,N_y,N_z,step);
        
        E_tot_line = diag(E_tot(:,:,N_z/2));
        plot_line(E_tot_line,delta_x*(1:N_x),step);
        tempvar = 0;
    end
    
    
    % H-field calculate
    ii = 2:N_x-1;
    jj = 2:N_y-1;
    kk = 2:N_z-1;
    Hx_new(ii,jj,kk) = D_a(ii,jj,kk).*Hx_old(ii,jj,kk) ...
        + D_b(ii,jj,kk).*(Ey_old(ii,jj,kk+1)-Ey_old(ii,jj,kk) ...
        + Ez_old(ii,jj,kk) - Ez_old(ii,jj+1,kk));

    ii = 2:N_x-1;
    jj = 2:N_y-1;
    kk = 2:N_z-1;
    Hy_new(ii,jj,kk) = D_a(ii,jj,kk).*Hy_old(ii,jj,kk) ...
        + D_b(ii,jj,kk).*(Ez_old(ii+1,jj,kk)-Ez_old(ii,jj,kk) ...
        + Ex_old(ii,jj,kk) - Ex_old(ii,jj,kk+1));

    ii = 2:N_x-1;
    jj = 2:N_y-1;
    kk = 2:N_z-1;
    Hz_new(ii,jj,kk) = D_a(ii,jj,kk).*Hz_old(ii,jj,kk) ...
        + D_b(ii,jj,kk).*(Ex_old(ii,jj+1,kk)-Ex_old(ii,jj,kk) ...
        + Ey_old(ii,jj,kk) - Ey_old(ii+1,jj,kk));

   
    % H-field increment
    Hx_old = Hx_new;
    Hy_old = Hy_new;
    Hz_old = Hz_new;

    % E-field calculate
    ii = 3:N_x-2;
    jj = 3:N_y-2;
    kk = 3:N_z-2;
    Ex_new(ii,jj,kk) = C_a(ii,jj,kk).*Ex_old(ii,jj,kk) ...
        + C_b(ii,jj,kk).*(Hz_old(ii,jj,kk)-Hz_old(ii,jj-1,kk) ...
        + Hy_old(ii,jj,kk-1) - Hy_old(ii,jj,kk) ...
        + Jsource_x(ii,jj,kk).*delta_x);

    ii = 3:N_x-2;
    jj = 3:N_y-2;
    kk = 3:N_z-2;
    Ey_new(ii,jj,kk) = C_a(ii,jj,kk).*Ey_old(ii,jj,kk) ...
        + C_b(ii,jj,kk).*(Hx_old(ii,jj,kk)-Hx_old(ii,jj,kk-1) ...
        + Hz_old(ii-1,jj,kk) - Hz_old(ii,jj,kk) ...
        + Jsource_y(ii,jj,kk).*delta_x);

    ii = 3:N_x-2;
    jj = 3:N_y-2;
    kk = 3:N_z-2;
    Ez_new(ii,jj,kk) = C_a(ii,jj,kk).*Ez_old(ii,jj,kk) ...
        + C_b(ii,jj,kk).*(Hy_old(ii,jj,kk)-Hy_old(ii-1,jj,kk) ...
        + Hx_old(ii,jj-1,kk) - Hx_old(ii,jj,kk) ...
        + Jsource_z(ii,jj,kk).*delta_x);

    % E-field Boundary ConditionsW_old
    c_ = max(c);   
    
    %planes
    Ex_new = assist_mur_abc_plane(c_, delta_t, delta_x, N_x, N_y, N_z, offset, Ex_new, Ex_old, Ex_old_old);
    Ey_new = assist_mur_abc_plane(c_, delta_t, delta_x, N_x, N_y, N_z, offset, Ey_new, Ey_old, Ey_old_old);
    Ez_new = assist_mur_abc_plane(c_, delta_t, delta_x, N_x, N_y, N_z, offset, Ez_new, Ez_old, Ez_old_old);
    
    %lines
    Ex_new = assist_mur_abc_line(c_, delta_t, delta_x, N_x, N_y, N_z, offset, pi/4, Ex_new, Ex_old);
    Ey_new = assist_mur_abc_line(c_, delta_t, delta_x, N_x, N_y, N_z, offset, pi/4, Ey_new, Ey_old);
    Ez_new = assist_mur_abc_line(c_, delta_t, delta_x, N_x, N_y, N_z, offset, pi/4, Ez_new, Ez_old);

    %points
%     Ex_new = assist_mur_abc_point(c_, delta_t, delta_x, N_x, N_y, N_z, offset, pi/4, Ex_new, Ex_old);
%     Ey_new = assist_mur_abc_point(c_, delta_t, delta_x, N_x, N_y, N_z, offset, pi/4, Ey_new, Ey_old);
%     Ez_new = assist_mur_abc_point(c_, delta_t, delta_x, N_x, N_y, N_z, offset, pi/4, Ez_new, Ez_old);

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
fprintf('simulation end\n');

function W_new_ = assist_mur_abc_plane(c_, delta_t, delta_x, N_x, N_y, N_z, offset, W_new, W_old, W_old_old)
    
    n_offset = offset-1;

    ii = offset+1:N_x-n_offset-1;
    jj = offset+1:N_y-n_offset-1;
    kk = offset+1:N_z-n_offset-1;

    W_new_ = W_new;

    W_new_(offset,jj,kk) = mur_abc_plane('x', c_, delta_t, delta_x, ...
        offset, offset+1, jj, kk, W_new, W_old, W_old_old);
    W_new_(N_x-n_offset,jj,kk) = mur_abc_plane('x', c_, delta_t, delta_x, ...
        N_x-n_offset, N_x-n_offset-1, jj, kk, W_new, W_old, W_old_old);
    W_new_(ii,offset,kk) = mur_abc_plane('y', c_, delta_t, delta_x, ...
        offset, offset+1, jj, kk, W_new, W_old, W_old_old);
    W_new_(ii,N_y-n_offset,kk) = mur_abc_plane('y', c_, delta_t, delta_x, ...
        N_x-n_offset, N_x-n_offset-1, jj, kk, W_new, W_old, W_old_old);
    W_new_(ii,jj,offset) = mur_abc_plane('z', c_, delta_t, delta_x, ...
        offset, offset+1, jj, kk, W_new, W_old, W_old_old);
    W_new_(ii,jj,N_z-n_offset) = mur_abc_plane('z', c_, delta_t, delta_x, ...
        N_x-n_offset, N_x-n_offset-1, jj, kk, W_new, W_old, W_old_old);

end

function Wz_new_ = mur_abc_plane(boundary, c, delta_t, delta_x, i0, i1, jj, kk, W_new, W_old, W_old_old)
    
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
        (c*delta_t-delta_x)/(c*delta_t+delta_x) ...
        (2*delta_x)/(c*delta_t+delta_x) ...
        (c*delta_t).^2/(2*delta_x*(c*delta_t+delta_x))];

    Wz_new_ = coeffs(1)* W_old_old(i1,jj,kk) ...
        +coeffs(2)* (W_new(i1,jj,kk)+W_old_old(i0,jj,kk)) ...
        +coeffs(3)* (W_old(i0,jj,kk)+W_old(i1,jj,kk)) ...
        +coeffs(4)* (W_old(i0,jj+1,kk) -4*W_old(i0,jj,kk) +W_old(i0,jj-1,kk) ...
        +W_old(i1,jj+1,kk) -4*W_old(i1,jj,kk) +W_old(i1,jj-1,kk)...
        +W_old(i0,jj,kk+1) +W_old(i0,jj,kk-1) +W_old(i1,jj,kk+1) +W_old(i1,jj,kk-1));

    if boundary == 'x'
        %default
    elseif boundary == 'y'
        Wz_new_ = permute(Wz_new_,[2 1 3]);       
    elseif boundary == 'z'
        Wz_new_ = permute(Wz_new_,[3 2 1]);
    end
end

function W_new_ = assist_mur_abc_line(c, delta_t, delta_x, N_x, N_y, N_z, offset, a, W_new, W_old)
    ofs = offset;
    nfs = offset-1;
    iii = ofs:N_x-nfs;
    jjj = ofs:N_y-nfs;
    kkk = ofs:N_z-nfs;
    
    W_new_ = W_new;
    coords_0 = ...
        {{iii, ofs, ofs}, ...
        {iii, N_y-nfs, ofs}, ...
        {ofs, jjj, ofs}, ...
        {N_x-nfs, jjj, ofs}, ...
        ...
        {ofs, ofs, kkk}, ...
        {ofs, N_y-nfs, kkk}, ...
        {N_x-nfs, N_y-nfs, kkk}, ...
        {N_x-nfs, ofs, kkk}, ...
        ...
        {iii, ofs, N_z-nfs}, ...
        {iii, N_y-nfs, N_z-nfs}, ...
        {ofs, jjj, N_z-nfs}, ...
        {N_x-nfs, jjj, N_z-nfs}};
    
    coords_1 = coords_0;
    for ii = 1:length(coords_0)
        for jj = [1:3]
            if length(coords_0{ii}) ~= 1
                coords_1{ii}{jj} = coords_0{ii}{jj};
            elseif coords_0{ii} == ofs
                coords_1{ii}{jj} = coords_0{ii}{jj}+1;
            else
                coords_1{ii}{jj} = coords_0{ii}{jj}-1;    
            end  
        end
    end
    
    for ii = 1:length(coords_0)
        W_new_(coords_0{ii}{1}, coords_0{ii}{2}, coords_0{ii}{3}) = ...
            mur_abc_point(c, delta_t, delta_x, a, ...
            coords_0{ii}{1}, coords_0{ii}{2}, coords_0{ii}{3}, coords_1{ii}{1}, coords_1{ii}{2}, coords_1{ii}{3}, W_new, W_old);
    end
    
end

function W_new_ = assist_mur_abc_point(c, delta_t, delta_x, N_x, N_y, N_z, offset, a, W_new, W_old)
    n_offset = offset-1;

    W_new_ = W_new;

    for z_count = [1,2]
        if (z_count==1)
            z_off_0 = offset;
            z_off_1 = offset+1;
        else
            z_off_0 = N_z-n_offset;
            z_off_1 = N_z-n_offset-1;
        end

        W_new_(offset, offset, z_off_0) = mur_abc_point(c, delta_t, delta_x, a, ...
            offset, offset, z_off_0, offset+1, offset+1, z_off_1, W_new, W_old);
        W_new_(N_x-n_offset, offset, z_off_0) = mur_abc_point(c, delta_t, delta_x, a, ...
            N_x-n_offset, offset, z_off_0, N_x-n_offset-1, offset+1, z_off_1, W_new, W_old);
        W_new_(N_x-n_offset, N_y-n_offset, z_off_0) = mur_abc_point(c, delta_t, delta_x, a, ...
            N_x-n_offset, N_y-n_offset, z_off_0, N_x-n_offset-1, N_y-n_offset-1, z_off_1, W_new, W_old);
        W_new_(offset, N_y-n_offset, z_off_0) = mur_abc_point(c, delta_t, delta_x, a, ...
            offset, N_y-n_offset, z_off_0, offset+1, N_y-n_offset, z_off_1, W_new, W_old);
    end
end

function W_new_ = mur_abc_point(c, delta_t, delta_x, a, ii_0, jj_0, kk_0, ii_1, jj_1, kk_1, W_new, W_old)
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
     ylim([0 1]);
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

function plot_field(field,N_x,N_y,N_z,step)
    
    figure(1)
    colormap jet
    xslice = N_x/2;   
    yslice = N_y/2;
    zslice = N_z/2;
    s = size(field);
    field(s(1)+1,s(2)+1,s(3)+1) = 0;

    h = slice(field,xslice,yslice,zslice);
    set(h,'edgecolor','none')

    name = sprintf('n = %d',step);
    title(name);
    xlabel('x');
    ylabel('y');
    zlabel('z');
    xlim([0 N_x+1])
    ylim([0 N_x+1])
    zlim([0 N_x+1])

%     set(gca,'ColorScale','log')

    bar = colorbar();
    ylabel(bar,'|H_{tot}| (A/m)');
        
    grid on   

    clim([0 0.05]);
%     view([1 0 0])
end