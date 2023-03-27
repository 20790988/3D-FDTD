function model = import_model(N_x,N_y,N_z,delta_x,delta_y,delta_z)
    
    size = [N_x,N_y,N_z];
    model = ones(size);

    % mid = N_x/2;
    for ii = 60:80
        for jj = 25:75
            for kk = 40:60
               model(ii,jj,kk) = 3;
            end
        end
    end
end