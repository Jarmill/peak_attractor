function [out_sim] = attractor_sim(dynamics, x0, Tmax, nw, odefcn)
%ATTRACTOR_SIM Simulate the path of a trajectory in an attractor

if nargin < 4
    nw = 0;
end

if nargin < 5
    odefcn = @ode45;
end

%deal with uncertain parameters later, assume
if dynamics.discrete
    X = zeros(length(x0), Tmax+1);
    X(1, :) = x0;
    for i = 1:Tmax
        xnext = out.dynamics(X(i, :));
        X(i+1, :) = xnext;
    end
else
    if length(Tmax)==2
        [T,X] = odefcn(dynamics.f{1},Tmax,x0);
    else
        [T,X] = odefcn(dynamics.f{1},[0,Tmax],x0);
    end
end



out_sim = struct;
out_sim.t = T;
out_sim.x = X;
out_sim.cost = dynamics.cost(X');
out_sim.nonneg = dynamics.nonneg(X');
end

