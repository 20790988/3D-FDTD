% %Cohn 1954
% 
% %strip width
% w = 4.5e-6;
% % distance between top & bottom plates
% b = 0.6e-6;
% % strip thickness
% t = 0.2e-6;
% 
% e_r = 4.6;
% epsilon_0 = 8.8542e-12;
% 
%     
% mu_0 = 1.2566e-6;
% 
% if ~(w/(b-t) >= 0.35)
%     sprintf('Calulation not valid!')
%     return;
% end
% 
% c_f = (0.0885*e_r/pi) * (2/(1-t/b)* log(1/(1-t/b) + 1) -...
%     (1/(1-t/b)* log(1/(1-t/b) - 1))); 
% 
% Z_0 = 94.15/(sqrt(e_r)*((w/b)/(1-t/b) +c_f/(0.0885*e_r)));
% 
% v = 1/sqrt(mu_0*epsilon_0*e_r);
% 
% b_cm = b*1e-2/1e-6;
% fc = 15/(b_cm*sqrt(e_r)*(w/b +pi/4));

source{1} = @one;
source{2} = @two;

fprintf('%d\n',source{1}());

function answ = one(~)
    answ = 1;
end

function answ = two(~)
    answ = 2;
end