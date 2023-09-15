classdef Simulation < handle
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here

    properties
        materials = [0; 0; 1; 1]
        N
        M = {100,100,100}
        delta = {1e-3,1e-3,1e-3}
        border_material = 1;
        align_with_grid = true;
    end

    methods
        function obj = Simulation(M,delta)
            obj.M = M;
            obj.delta = delta;          
        end
        
        function obj = setMaterials(obj,sigma,sigma_m,epsilon_r,mu_r)
            obj.materials = [sigma;sigma_m;epsilon_r;mu_r];
        end
    end
end