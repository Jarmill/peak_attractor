SOLVE = 1;
PLOT = 1;


if SOLVE
    %lorenz attractor
    mpol('x', 2, 1);

    %parameters of Lorenz system
    a = 1.4; b = 0.3;

    %support space
    Xsupp = [];

    %dynamics
    f0 = [  1 - a*x(1)^2 + x(2) ; b*x(1)  ];
    X = [];

    %objective
%     p0 = x(1);
    p0 = -x'*x;
    % p0 = x'*x;
    % p0 = [x(1)^2; x(3)^2];
    % p0 = [x(1:2)'*x(1:2)];

    q = [1/1.5;1];
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
    p_opt.discount = 0.2;

    p_opt.discrete = 1;

    order = 4;
    out = peak_attract(p_opt, order);
    peak_val = out.peak_val;
    xp = out.xp ./ q;
end

if PLOT
    Tmax_sim = 4000;
    x0 = q.*[1; 0];
%     [T,X] = ode45(f_ode,[0,10000],x0); X = X(round(size(X,1) / 5):end,:);
    out_sim = {attractor_sim(out.dynamics, x0, Tmax_sim)};
    if out.recover == 1
        %out_sim_peak = attractor_sim(out.dynamics, out.xp, Tmax_sim/2 * [-1, 1]);
        out_sim_peak = {attractor_sim(out.dynamics, out.xp, 0.8*Tmax_sim, 0)};
%         out_sim_peak.t = out_sim_peak.t + Tmax_sim/2;
        %out_sim_peak.tp = Tmax_sim/2;
        out_sim_peak{1}.tp = 0;
%         out_sim_peak = switch_sampler(out.dynamics, out.xp, 1, Tmax_sim);


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
% syms y [3, 1];