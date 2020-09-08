function [out_sim] = attractor_sim(dynamics, x0, Tmax, Tstart, nw, odefcn)
%ATTRACTOR_SIM Simulate the path of a trajectory in an attractor

if nargin < 4
    Tstart = 0.2*Tmax;
end

if nargin < 5
    nw = 0;
end

if nargin < 6
    odefcn = @ode45;
end

%deal with uncertain parameters later, assume
if dynamics.discrete
    X = zeros(Tmax+1, length(x0));
    X(1, :) = x0;
    if nw > 0
        %includes uncertain fixed parameters
        curr_f = @(t, x) dynamics.f{1}(t, x, w0);        
    else
        %no uncertain fixed parameters
        curr_f = dynamics.f{1};        
    end
    for i = 1:Tmax
        xnext = curr_f(i, X(i, :)');
        X(i+1, :) = xnext;
    end
    T = 0:Tmax;
    Tstart = ceil(Tstart);
else
    if length(Tmax)==2
        [T,X] = odefcn(dynamics.f{1},Tmax,x0);
    else
        [T,X] = odefcn(dynamics.f{1},[0,Tmax],x0);
    end
end

%time shift to get on the attractor
X = X(T >= Tstart, :);
T = T(T >= Tstart);
T = T - Tstart;

out_sim = struct;
out_sim.t = T;
out_sim.x = X;
if nw > 0
    out_sim.w = w0;
else
    out_sim.w = [];
end
out_sim.cost = dynamics.cost(X');
%out_sim.nonneg = dynamics.nonneg(X');

if isfield(dynamics, 'nonneg')    
    if nw > 0        
        out_sim.nonneg = dynamics.nonneg(X', w0);        
    else
        out_sim.nonneg = dynamics.nonneg(X');        
    end
end

end

