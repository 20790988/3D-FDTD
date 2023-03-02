%
%  3D-FDTD
%  K Meyer
%  2023-02-15
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
    N_x = 100;
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
    
    Jsource_z(source_x,source_y,source_z) = source_signal(step+1);

    if (step == 50) || (step == 100) || (step == 150)
        H_tot = sqrt(Hx_old.^2+Hy_old.^2+Hz_old.^2);
        plot_field(H_tot,N_x,N_y,N_z,step);
        
        H_tot_line = H_tot(1:N_x,N_y/2,N_z/2);
        plot_line(H_tot_line,delta_x*(1:N_x),step);
        tempvar = 0;
    end
    
    
    % H-field calculate
    ii = 1:N_x;
    jj = 1:N_y-1;
    kk = 1:N_z-1;
    Hx_new(ii,jj,kk) = D_a(ii,jj,kk).*Hx_old(ii,jj,kk) ...
        + D_b(ii,jj,kk).*(Ey_old(ii,jj,kk+1)-Ey_old(ii,jj,kk) ...
        + Ez_old(ii,jj,kk) - Ez_old(ii,jj+1,kk));

    ii = 1:N_x-1;
    jj = 1:N_y;
    kk = 1:N_z-1;
    Hy_new(ii,jj,kk) = D_a(ii,jj,kk).*Hy_old(ii,jj,kk) ...
        + D_b(ii,jj,kk).*(Ez_old(ii+1,jj,kk)-Ez_old(ii,jj,kk) ...
        + Ex_old(ii,jj,kk) - Ex_old(ii,jj,kk+1));

    ii = 1:N_x-1;
    jj = 1:N_y-1;
    kk = 1:N_z;
    Hz_new(ii,jj,kk) = D_a(ii,jj,kk).*Hz_old(ii,jj,kk) ...
        + D_b(ii,jj,kk).*(Ex_old(ii,jj+1,kk)-Ex_old(ii,jj,kk) ...
        + Ey_old(ii,jj,kk) - Ey_old(ii+1,jj,kk));

   
    % H-field increment
    Hx_old = Hx_new;
    Hy_old = Hy_new;
    Hz_old = Hz_new;

    % E-field calculate
    ii = 1:N_x;
    jj = 2:N_y;
    kk = 2:N_z;
    Ex_new(ii,jj,kk) = C_a(ii,jj,kk).*Ex_old(ii,jj,kk) ...
        + C_b(ii,jj,kk).*(Hz_old(ii,jj,kk)-Hz_old(ii,jj-1,kk) ...
        + Hy_old(ii,jj,kk-1) - Hy_old(ii,jj,kk) ...
        + Jsource_x(ii,jj,kk).*delta_x);

    ii = 2:N_x;
    jj = 1:N_y;
    kk = 2:N_z;
    Ey_new(ii,jj,kk) = C_a(ii,jj,kk).*Ey_old(ii,jj,kk) ...
        + C_b(ii,jj,kk).*(Hx_old(ii,jj,kk)-Hx_old(ii,jj,kk-1) ...
        + Hz_old(ii-1,jj,kk) - Hz_old(ii,jj,kk) ...
        + Jsource_y(ii,jj,kk).*delta_x);

    ii = 2:N_x;
    jj = 2:N_y;
    kk = 1:N_z;
    Ez_new(ii,jj,kk) = C_a(ii,jj,kk).*Ez_old(ii,jj,kk) ...
        + C_b(ii,jj,kk).*(Hy_old(ii,jj,kk)-Hy_old(ii-1,jj,kk) ...
        + Hx_old(ii,jj-1,kk) - Hx_old(ii,jj,kk) ...
        + Jsource_z(ii,jj,kk).*delta_x);

     % E-field Boundary Condition
     ii = 2:N_x-1;
     jj = 2:N_y-1;
     kk = 2:N_z-1;
     c_ = max(c);   

     Ex_new(ii,1,kk) = mur_abc('y', c_, delta_t, delta_x, 1, 2, jj, kk, Ex_new, Ex_old, Ex_old_old);
     Ex_new(ii,N_y,kk) = mur_abc('y', c_, delta_t, delta_x, N_x, N_x-1, jj, kk, Ex_new, Ex_old, Ex_old_old);
     Ex_new(ii,jj,1) = mur_abc('z', c_, delta_t, delta_x, 1, 2, jj, kk, Ex_new, Ex_old, Ex_old_old);
     Ex_new(ii,jj,N_z) = mur_abc('z', c_, delta_t, delta_x, N_x, N_x-1, jj, kk, Ex_new, Ex_old, Ex_old_old);
         
     Ey_new(1,ii,jj) = mur_abc('x', c_, delta_t, delta_x, 1, 2, jj, kk, Ey_new, Ey_old, Ey_old_old);
     Ey_new(N_x,ii,jj) = mur_abc('x', c_, delta_t, delta_x, N_x, N_x-1, jj, kk, Ey_new, Ey_old, Ey_old_old);
     Ey_new(ii,jj,1) = mur_abc('z', c_, delta_t, delta_x, 1, 2, jj, kk, Ey_new, Ey_old, Ey_old_old);
     Ey_new(ii,jj,N_z) = mur_abc('z', c_, delta_t, delta_x, N_x, N_x-1, jj, kk, Ey_new, Ey_old, Ey_old_old);
    
     
     Ez_new(1,ii,jj) = mur_abc('x', c_, delta_t, delta_x, 1, 2, jj, kk, Ez_new, Ez_old, Ez_old_old);
     Ez_new(N_x,ii,jj) = mur_abc('x', c_, delta_t, delta_x, N_x, N_x-1, jj, kk, Ez_new, Ez_old, Ez_old_old);
     Ez_new(ii,1,kk) = mur_abc('y', c_, delta_t, delta_x, 1, 2, jj, kk, Ez_new, Ez_old, Ez_old_old);
     Ez_new(ii,N_y,kk) = mur_abc('y', c_, delta_t, delta_x, N_x, N_x-1, jj, kk, Ez_new, Ez_old, Ez_old_old);
     
%       mur_abc('x', c_, delta_t, delta_x, 1, 2, jj, kk, W_new, W_old, W_old_old);
%       mur_abc('x', c_, delta_t, delta_x, N_x, N_x-1, jj, kk, W_new, W_old, W_old_old);


    Ex_new([1,N_x],[1,N_y],[1,N_z]) = 0;
    Ey_new([1,N_x],[1,N_y],[1,N_z]) = 0;
    Ez_new([1,N_x],[1,N_y],[1,N_z]) = 0;

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

function Wz_new_ = mur_abc(boundary,c, delta_t, delta_x, i0, i1, jj, kk, W_new, W_old, W_old_old)
    
    if boundary == 'x'
        %default
    elseif boundary == 'y'
        permute(W_new,[2 1 3]);
        permute(W_old,[2 1 3]);
        permute(W_old_old,[2 1 3]);
    elseif boundary == 'z'
        permute(W_new,[3 2 1]);
        permute(W_old,[3 2 1]);
        permute(W_old_old,[3 2 1]);
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
        permute(Wz_new_,[2 1 3]);       
    elseif boundary == 'z'
        permute(Wz_new_,[3 2 1]);
    end
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
    ylim([0 1e-3]);
end

function plot_field_slice(field)
    
    figure(1)
    colormap jet
    h = surf(field);
    set(h,'edgecolor','none')
%     set(gca,'ColorScale','log')
    view([0 0 1])
    colorbar();
    clim([-0.01 0.01]);
end

function plot_field(field,N_x,N_y,N_z,step)
    
    figure(1)
    colormap gray
    xslice = N_x/2;   
    yslice = N_y/2;
    zslice = N_z/2;
    h = slice(field,xslice,yslice,zslice);
    set(h,'edgecolor','none')

    name = sprintf('n = %d',step);
    title(name);
    xlabel('x');
    ylabel('y');
    zlabel('z');
%     set(gca,'ColorScale','log')

    bar = colorbar();
    ylabel(bar,'|H_{tot}| (A/m)');
        
    grid on   

    clim([0 0.05e-3]);
%     view([1 0 0])
end

function plot_fields(field1,field2,field3,N_x,N_y,N_z,fig_no)
    
    figure(fig_no)
    subplot(3,1,1);
    colormap jet
    xslice = N_x/2;   
    yslice = N_y/2;
    zslice = [];
    h = slice(field1,xslice,yslice,zslice);
    set(h,'edgecolor','none')
    colorbar();
    clim([-0.01 0.01]);

    subplot(3,1,2);
    colormap jet
    xslice = N_x/2;   
    yslice = N_y/2;
    zslice = [];
    h = slice(field2,xslice,yslice,zslice);
    set(h,'edgecolor','none')
    colorbar();
    clim([-0.01 0.01]);

    subplot(3,1,3);
    colormap jet
    xslice = N_x/2;   
    yslice = N_y/2;
    zslice = [];
    h = slice(field3,xslice,yslice,zslice);
     set(h,'edgecolor','none')
    colorbar();
    clim([-0.01 0.01]);
end
