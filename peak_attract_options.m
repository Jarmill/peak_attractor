classdef peak_attract_options < handle
    %attributes of peak estimation routine for attractor
    %a default set of options
    properties
        
        %% properties of run
        %discount factor > 0 
        discount(1,1) double{mustBePositive}  = 1;           
        
        %discrete (true) or continuous (false) system
        discrete(1,1) logical = false;
        
        %function dynamics
        %f: dynamics
        %X: space on which dynamics are valid (arbitrary switching)
        dynamics = struct('f', [], 'X', {})
        
        
        %objective to minimize
        %could be a function (single objective)
        %or a cell of functions (minimum of entries)        
        obj = [];
        
        %% Variables and descriptors
        %variables (array of symbolic variables)
%         x mpol = [];     %states
%         w mpol= [];     %uncertain parameters
        var = struct('x', [], 'w', []);
        
        %% support sets
        %type @mpol/supcon       
        state_supp = []; %states to consider, set X
         
        %Coordinate ranges for variables for scaling
        %state variables lie in a box (utils/box_process)
        
        %box is scalar B:           -B      <= x_i <= B
        %box is [Bmin, Bmax]:       -Bmin   <= x_i <= Bmax
        %box is [B_i]:              -Bi     <= x_i <= B_i
        %box is [Bmin_i, Bmax_i]    -Bmin_i <= x_i <= Bmax_i
        box = 1; 
        
        scale = 0; %should variables be scaled to [-1,1] (state)
                        
        param = [];  %range of parameters w  
        
        %% additional options
        solver = 'mosek';
        
        %what is the tolerance for identifying a matrix as rank-1?
        rank_tol(1,1) double{mustBePositive} = 1e-4; 
        
        
    end
end