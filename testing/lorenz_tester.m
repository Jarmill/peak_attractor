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

q = 1./[25; 30; 50];
f = q.*subs(f0, x, x.*(1./q));
objective = subs(p0, x, x.*(1./q));


p_opt = peak_attract_options;
p_opt.var.x = x;
p_opt.dynamics = struct;
p_opt.dynamics.f = f;
p_opt.dynamics.X = X;

p_opt.state_supp = Xsupp;

p_opt.obj = objective;


order = 4;
out = peak_attract(p_opt, order);
peak_val = out.peak_val;
xp = out.xp ./ q;
end

if PLOT
    Tmax_sim = 100;
    x0 = [0.5;0.5;0.6];
%     [T,X] = ode45(f_ode,[0,10000],x0); X = X(round(size(X,1) / 5):end,:);
    out_sim = {attractor_sim(out.dynamics, x0, Tmax_sim)};
    if out.recover == 1
        %out_sim_peak = attractor_sim(out.dynamics, out.xp, Tmax_sim/2 * [-1, 1]);
        out_sim_peak = {attractor_sim(out.dynamics, out.xp, Tmax_sim)};
%         out_sim_peak.t = out_sim_peak.t + Tmax_sim/2;
        %out_sim_peak.tp = Tmax_sim/2;
        out_sim_peak{1}.tp = 0;
%         out_sim_peak = switch_sampler(out.dynamics, out.xp, 1, Tmax_sim);


        nplot = nonneg_plot(out, out_sim, out_sim_peak);
        cplot = cost_plot(out, out_sim, out_sim_peak);
        splot = state_plot_3_cont(out, out_sim, out_sim_peak);
    else
        %nothing 
        nplot = nonneg_plot(out, out_sim);
        cplot = cost_plot(out, out_sim);
        splot = state_plot_3_cont(out, out_sim);
        view(-168, 20);
    end
end
% syms y [3, 1];