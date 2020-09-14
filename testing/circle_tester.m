SOLVE = 1;
PLOT = 1;

if SOLVE
%van der pol attractor
mpol('x', 2, 1);

%support space
Xsupp = [];

%dynamics
f0 = [-0.5*x(2); 0.5*x(1)];

%f0 = [ 2*x(2,:) ; -0.8*x(1,:) - 10*(x(1,:).^2-0.21).*x(2,:) ];
X = [];

% x0 = [2; 0.5];
% x0 = [-1.3; 1];
x0 = [-0.5; -0.5];

%objective

f = f0;
% objective = x'*x;
% objective = x(1)^2;
objective = x(1);
% objective = x(1)^2 + 2*x(2)^2 + 0.1*x(1)*x(2);

% objective = [x(1)^2; x(2)^2];

p_opt = peak_attract_options;
p_opt.var.x = x;
p_opt.dynamics = struct;
p_opt.dynamics.f = f;
p_opt.dynamics.X = X;

p_opt.state_supp = Xsupp;
p_opt.state_init = (x == x0);
p_opt.box = 2;

p_opt.discount = 1;
p_opt.obj = objective;


order = 8;
out = peak_attract(p_opt, order);
peak_val = out.peak_val;
xp = out.xp;
end

if PLOT
    Tmax_sim = 100;
    %x0 = [0.486;-0.2608];
%     x0 = [2, 0.5];
    out_sim = {attractor_sim(out.dynamics, x0, Tmax_sim)};
    if out.recover == 1
        %approximate initial condition recovered
        out_sim_peak = {attractor_sim(out.dynamics, out.xp, Tmax_sim)};
        out_sim_peak{1}.tp = 0;

        nplot = nonneg_plot(out, out_sim, out_sim_peak);
        cplot = cost_plot(out, out_sim, out_sim_peak);
        splot = state_plot_2(out, out_sim, out_sim_peak);
    else
        %just plot sample trajectories
        nplot = nonneg_plot(out, out_sim);
        cplot = cost_plot(out, out_sim);
        splot = state_plot_2(out, out_sim);
    end
end