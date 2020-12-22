SOLVE = 1;
PLOT = 1;

if SOLVE
%van der pol attractor
mpol('x', 2, 1);

%support space
Xsupp = [];

%dynamics
f0 = [ 2*x(2,:) ; -0.8*x(1,:) - 10*(x(1,:).^2-0.21).*x(2,:) ];
X = [];

%objective

f = f0;
% objective = x'*x;
objective = x(1)^2;
% objective = x(1)^2 + 2*x(2)^2 + 0.1*x(1)*x(2);

% objective = [x(1)^2; x(2)^2];

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

Mp_mix = out.Mp(2:3, 2:3);
rank_sym_p = rank(Mp_mix, p_opt.rank_tol);

if rank_sym_p==1
    sym_signp = sign(Mp_mix(2,1));    
    xp_abs = sqrt([Mp_mix(1,1); Mp_mix(2,2)]);
    xp = xp_abs * [1 sym_signp];
    
    out.recover = 1;
    out.xp = xp;
end



end

if PLOT

    
    Tmax_sim = 150;
    %x0 = [0.486;-0.2608];
    x0 = [0.3; 0.0];
    out_sim = {attractor_sim(out.dynamics, x0, Tmax_sim)};
    if out.recover == 1
        %approximate initial condition recovered        
        out_sim_peak = {attractor_sim(out.dynamics, out.xp(:, 1), Tmax_sim), ...
                        attractor_sim(out.dynamics, out.xp(:, 2), Tmax_sim)};
        out_sim_peak{1}.tp = 0;
        out_sim_peak{2}.tp = 0;
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