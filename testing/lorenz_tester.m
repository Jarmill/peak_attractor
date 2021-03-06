SOLVE = 0;
PLOT = 1;


if SOLVE
%lorenz attractor
mpol('x', 3, 1);

%parameters of Lorenz system
rho = 28;
sigma = 10;
beta = 8/3;

%support space
Xsupp = [];

%dynamics
f0 = [10*(x(2)-x(1)) ; x(1).*(28-x(3)) - x(2) ; x(1).*x(2) - (8/3)*x(3) ];
X = [];

%objective
p0 = x(1);
% p0 = x'*x;
% p0 = x(2)^2;
% p0 = [x(1)^2; x(3)^2];
% p0 = [x(1:2)'*x(1:2)];

q = 1./[25; 30; 50];
f = q.*subs(f0, x, x.*(1./q));
objective = subs(p0, x, x.*(1./q));
% objective = x(1:2)'*x(1:2);


p_opt = peak_attract_options;
p_opt.var.x = x;
p_opt.dynamics = struct;
p_opt.dynamics.f = f;
p_opt.dynamics.X = X;

p_opt.box = [-0.85, 0.85; -0.9, 0.9; 0.01, 0.91];

p_opt.state_supp = Xsupp;

p_opt.discount = 1;
p_opt.obj = objective;


% order = 4;
order = 6;
% order = 8;
out = peak_attract(p_opt, order);
peak_val = out.peak_val;
xp = out.xp ./ q;
end

if PLOT
    Tmax_sim = 100;
    x0 = [0.5;0.5;0.6];

    out_sim = {attractor_sim(out.dynamics, x0, Tmax_sim)};
    if out.recover == 1
        %approximate initial condition recovered
        out_sim_peak = {attractor_sim(out.dynamics, out.xp, Tmax_sim)};
        out_sim_peak{1}.tp = 0;
        
%         dyn_rev= out.dynamics;
%         dyn_rev.f = {@(t, x) -out.dynamics.f{1}(t, x)};
%         out_sim_peak_rev = {attractor_sim(dyn_rev, out.xp, Tmax_sim, 0)};
%         out_sim_peak_rev{1}.tp = 0;

        nplot = nonneg_plot(out, out_sim, out_sim_peak);
        cplot = cost_plot(out, out_sim, out_sim_peak);
%         splot = state_plot_3(out, out_sim, out_sim_peak, out_sim_peak_rev);
        splot = state_plot_3(out, out_sim, out_sim_peak);
        
        xlim([-1, 1])
        ylim([-1, 1])
        zlim([0, 1])
        view(-152, 7)
    else
        %just plot sample trajectories
        nplot = nonneg_plot(out, out_sim);
        cplot = cost_plot(out, out_sim);
        splot = state_plot_3(out, out_sim);
        view(-168, 20);
    end
end
% syms y [3, 1];