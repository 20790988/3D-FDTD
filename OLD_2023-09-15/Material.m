classdef Material
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here

    properties
        sigma_e = 0;
        sigma_m = 0;
        epsilon_r = 1;
        mu_r = 0;
    end

    methods
        function obj = Material(n_x,n_y)
            %UNTITLED Construct an instance of this class
            %   Detailed explanation goes here
            obj.N_x = n_x;
            obj.N_y = n_y;
        end

        function outputArg = method1(obj,inputArg)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            outputArg = obj.Property1 + inputArg;
        end
    end
end