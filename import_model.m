function model = import_model(N_x,N_y,N_z,delta_x,delta_y,delta_z)
    
    size = [N_x,N_y,N_z];
    model = ones(size);

    % mid = N_x/2;
% for ii = 1:N_x
%     for jj = 1:N_y
%         
%         if sqrt((ii-mid)^2+(jj-mid)^2) <= 30
%             material(ii,jj) = 2;
%         end
%     end
% end
end