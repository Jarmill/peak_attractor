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


order = 3;
out = peak_attract(p_opt, order);
peak_val = out.peak_val;
xp = out.xp ./ q;