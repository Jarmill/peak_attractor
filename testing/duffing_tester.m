SOLVE = 1;
PLOT = 1;

if SOLVE
    %henon attractor
    mpol('x', 2, 1);

    %parameters of duffing system
    a = 2.75; b = 0.2;

    %support space
    Xsupp = [];

    %dynamics
    f0 = [x(2); -b*x(1) + a*x(2) - x(2)^3];
    X = [];

    %objective
%     p0 = x(1);
%     p0 = -x'*x;
p0 = -x(2);
%     p0 = [x(1); -x(2)];
    % p0 = x'*x;
    % p0 = [x(1)^2; x(3)^2];
    % p0 = [x(1:2)'*x(1:2)];

    q = 1./[1.75;1.75];
    f = q.*subs(f0, x, x.*(1./q));
    objective = subs(p0, x, x.*(1./q));
    % objective = x(1:2)'*x(1:2);


    p_opt = peak_attract_options;
    p_opt.var.x = x;
    p_opt.dynamics = struct;
    p_opt.dynamics.f = f;
    p_opt.dynamics.X = X;

    p_opt.state_supp = Xsupp;

    p_opt.obj = objective;
    p_opt.discount = 0.5;

    p_opt.discrete = 1;

    order = 4;
    out = peak_attract(p_opt, order);
    peak_val = out.peak_val;
    xp = out.xp ./ q;
end

if PLOT
    Tmax_sim = 40000;
    x0 = q.*[1; 0];

    out_sim = {attractor_sim(out.dynamics, x0, Tmax_sim)};
    if out.recover == 1
        %approximate initial condition recovered
        out_sim_peak = {attractor_sim(out.dynamics, out.xp, 0.8*Tmax_sim, 0)};
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
% syms y [3, 1];