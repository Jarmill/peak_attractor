SOLVE = 1;
PLOT = 1;
% rng(40, 'twister')

if SOLVE
%lorenz attractor
mpol('x', 2, 1);

%parameters of Lorenz system
% rho = 28;
% sigma = 10;
% beta = 8/3;

%support space
Xsupp = [];

%dynamics
f0 = [ 2*x(2,:) ; -0.8*x(1,:) - 10*(x(1,:).^2-0.21).*x(2,:) ];
X = [];

%objective

f = f0;
% objective = x'*x;
% objective = x(1)^2;
% objective = x(1)^2 + 2*x(2)^2 + 0.1*x(1)*x(2);

objective = [x(1)^2; x(2)^2];

p_opt = peak_attract_options;
p_opt.var.x = x;
p_opt.dynamics = struct;
p_opt.dynamics.f = f;
p_opt.dynamics.X = X;

p_opt.state_supp = Xsupp;

p_opt.discount = 1;
p_opt.obj = objective;


order = 6;
out = peak_attract(p_opt, order);
peak_val = out.peak_val;
xp = out.xp;
end

if PLOT
    Tmax_sim = 150;
    x0 = [0.486;-0.2608];
%     [T,X] = ode45(f_ode,[0,10000],x0); X = X(round(size(X,1) / 5):end,:);
    out_sim = {attractor_sim(out.dynamics, x0, Tmax_sim)};
    if out.recover == 1
        %out_sim_peak = attractor_sim(out.dynamics, out.xp, Tmax_sim/2 * [-1, 1]);
        out_sim_peak = {attractor_sim(out.dynamics, out.xp, Tmax_sim)};
        out_sim_peak{1}.tp = 0;

        nplot = nonneg_plot(out, out_sim, out_sim_peak);
        cplot = cost_plot(out, out_sim, out_sim_peak);
        splot = state_plot_2_cont(out, out_sim, out_sim_peak);
    else
        %nothing 
        nplot = nonneg_plot(out, out_sim);
        cplot = cost_plot(out, out_sim);
        splot = state_plot_2_cont(out, out_sim);
    end
end