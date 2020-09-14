function [out] = peak_attract(options, order)
%PEAK_ESTIMATE  Find the maximum value of a function p(x) along all
%trajectories inside an attractor with  polynomial dynamics x' = f(x) 
%(continuous) or x+ = f(x) (discrete)
%Input: options structure (peak_options) with fields:
%   var:        Structure of symbolic variables (mpol)
%       x:      state
%       w:      parametric uncertainty (default empty)
%
%   dynamics:   Structure with fields (f, X)
%       f:      dynamics x' = f(t,x) over the space X.
%                   Each entry is a polynomial                
%       X:      over what space do the dynamics evolve
%
%   obj:        Functions to maximize along trajectories
%       Scalar:     Single function
%       List:       Maximize the minimum of all objectives
%                      Each entry is a polynomial       
%   
%   state_supp: Support set of X  (by default the box in `box')
%   state_init: Support set of X0 (by default is X)
%   param:      Support of parametric uncertainty w (if var.w nonempty)
%       
%   rank_tol:   Rank tolerance for moment matrix to be rank-1
%
%Input (separate argument for convenience)
%   order:      Order of relaxation
%
%Output: out structure with fields
%   peak_val:   Value of estimated peak optimum
%   optimal:    Is peak_val an approximate optimum with rank-tolerance rank_tol
%
%   xp:         Optimal point
%   Mp:         Moment Matrix of Peak Measure
%   Mfw:        Moment Matrices of Forward Occupation Measure
%   Mbk:        Moment Matrices of Backward Occupation Measure
%
%   v:          Polynomial from dual program, v(t, x) - gamma >=0 invariant
%   Lv:         Lie derivatives along trajectories of system for v
%   theta:       multiple cost tradeoff at optimal point

if nargin < 2
    order = 3;
end
d = 2*order;

%% Process option data structure
mset clearmeas

%state variable
nx = length(options.var.x);
nvar = nx;

%parameter variables
nw = length(options.var.w);
nvar = nvar + nw;

%compact support (artificial)
%XR = (sum(options.var.x.^2) <= options.R^2);
%scaling of the box
xp = options.var.x;
Xsupp = options.state_supp;
%TODO: deal with parameters nx -> nx + nw;
if isempty(options.box) || (length(options.box)==1 && options.box == 0)
%     XR_scale = [];
%     XR_unscale = [];
    
%     xp_scale = xp;
%     xp_inv_scale = xp;
    
    XR = [];
else    
    [box, box_center, box_half] = box_process(nx, options.box);

%     XR_scale = (xp.^2) <= 1;
%     XR_unscale = (xp - box_center).^2 <= box_half.^2;
    XR = (xp - box_center).^2 <= box_half.^2;
%     Xsupp = [Xsupp; XR];
%     xp_scale = box_half.*xp + box_center;
%     xp_inv_scale = (xp - box_center) .* (1./box_half);
end
% Xsupp = options.state_supp;
% 
% 
% if options.scale
%     XR = XR_scale;
% else
%     XR = XR_unscale;
% end



%number of switching subsystems
f = options.dynamics.f;
X = options.dynamics.X;

if iscell(options.dynamics.f)
    nsys = length(options.dynamics.f);
else
    nsys = 1;
    options.dynamics.f = {options.dynamics.f};
    options.dynamics.X = {options.dynamics.X};    
    f = {f};
    X = {X};
end

%now scale everything
% if options.scale    
%     for i = 1:nsys
%         f{i} = subs((1./box_half) .*f {i}, xp, xp_scale);
%         
%         X{i} = [subs([X{i}; Xsupp], xp, xp_scale); XR_scale];
%     end
% else
    for i = 1:nsys
       X{i} = [X{i}; Xsupp; XR];
    end
% end

%deal with initial measure
INITIAL = ~isempty(options.state_init);
if INITIAL
    mpol('x0', nx, 1);

    X0 = subs([options.state_init; Xsupp; XR], options.var.x, x0);
    %TODO: deal with W param scaling
else
    x0 = [];
    X0 = [];
end


%number of objectives (1 standard, 2+ minimum)
nobj = length(options.obj);

%% Measures and variables
%measures in the problem:
%start with no time breaks, one initial measure
%one occupation measure per switched systems

%one peak measure

% if options.scale
% %     Xp = [subs(options.state_supp, xp, xp_scale); Xsupp; XR];
%     Xp = [subs(options.state_supp, xp, xp_scale)]; %that might be the bug (hopefully)
% else
    
    %this is where dual_rec breaks
%     Xp = [Xsupp; XR];
    %adding equality constraints here breaks the invariant function
    %dual_rec'monomials. Since X0 and X_occ are constrained within Xsupp,
    %this should not be a problem. If dynamics stay within Xsupp, then Xp
    %should be supported in Xsupp at optimality.
    Xp = [Xsupp; XR];
% end
%deal with hanging variables and measures by letting the original x be the
%peak measure

% %deal with parameters
if nw > 0
    %parameters
    wp = options.var.w;    
    
    W= options.param;
    
    if INITIAL
        mpol('w0', nw, 1);
        W0 = subs(W, wp, w0);
    else
        w0 = [];
        W0 = [];
    end
    
    %W = [Wp; W0; wp == w0] 
    
else
    W = [];    
    wp = [];
    w0 = [];
    W0 = [];
end


%occupation measures
mu_fw = cell(nsys, 1);  %forward
mu_bk = cell(nsys, 1);  %backward

% mon_occ = cell(nsys, 1);
Ay_fw = 0;
Ay_bk = 0;
%X = cell(nsys, 1);
% X_fw = []; %support
% X_bk = [];
X_occ = [];



%occupation measure
mpol('x_fw', nx, nsys);
mpol('x_bk', nx, nsys);
if nw > 0
    mpol('w_fw', nw, nsys);
    mpol('w_bk', nw, nsys);
else
    w_bk = zeros(0, nsys);
    w_fw = zeros(0, nsys);
end

%measure information         
mup = meas([xp; wp]);
monp = mmon([xp; wp], d);
yp = mom(monp);


for i = 1:nsys
    
    %xcurr = x_occ(:, i);
    
    %new variables
    var_fw = struct('x', x_fw(:, i), 'w', w_fw(:, i));
    var_bk = struct('x', x_bk(:, i), 'w', w_bk(:, i));

    %occupation_measure(f, X, var, var_new, d, forward, discrete, discount)
    [mu_fw{i}, mon_fw{i}, X_fw_curr, Ay_fw_curr] = occupation_measure(f{i},...
        X{i}, options.var, var_fw, d, 1, options.discrete, options.discount);
    
    [mu_bk{i}, mon_bk{i}, X_bk_curr, Ay_bk_curr] = occupation_measure(f{i},...
        X{i}, options.var, var_bk, d, 0, options.discrete, options.discount);
    
    %[mu_fw{i}, mon_occ{i}, X_occ_curr, Ay_curr] = occupation_measure(f{i}, ...
    %    X{i}, options.var, var_new, d);

    %mu_occ_sum = mu_occ_sum + mu_occ{i};
    X_occ = [X_occ; X_fw_curr; X_bk_curr];
    Ay_fw = Ay_fw + Ay_fw_curr;
    Ay_bk = Ay_bk + Ay_bk_curr;

    if nw > 0
        W_fw = subs(Wp, wp, w_fw(:, i));
        W_bk = subs(Wp, wp, w_bk(:, i));
        W = [W; W_fw; W_bk];
    end
end


%% Form Constraints and Solve Problem

supp_con = [Xp; X_occ; W];
%careful of monic substitutions ruining dual variables
Liou_fw = - yp + Ay_fw;
Liou_bk = - yp + Ay_bk;

if INITIAL
    mu0 = meas([x0; w0]);
    mon0 = mmon([x0; w0], d);
    y0 = mom(mon0);
    Liou_bk = Liou_bk + y0;
    supp_con = [supp_con; X0; W0];
end

mom_con = [mass(mup) == 1; Liou_fw == 0; Liou_bk == 0];

obj = options.obj;
if nobj == 1
    cost = obj;

%     if options.scale
%         cost = subs(obj, xp, xp_scale);
%     end
%     %cost = subs(options.obj, [options.var.t; options.var.x], [tp; xp]);
else
    %add a new measure for the cost
    %honestly I only want a free variable, so bounds on support are not
    %necessary. The free variable is the off-diagonal entry of a 2x2 psd
    %matrix with top corner 1. Inefficient, but that's how gloptipoly is
    %interfaced.
    
    mpol('c');
    muc = meas(c);
    momc = mom(c);
    cost = momc;
    
%     cost_mom_con = (mass(muc) == 1);
    cost_mom_con = (mass(muc) == 1);
%     cost_mom_con = [];
    for i = 1:nobj
        %curr_obj = subs(options.obj(i), var.x, xp);
        curr_obj = options.obj(i);
%         if options.scale
%             curr_obj = subs(options.obj(i), xp, xp_scale);
%         end
%         
        cost_mom_con = [cost_mom_con; (momc <= mom(curr_obj))];
        
    end 
    
    mom_con = [cost_mom_con; mom_con];
end

objective = max(cost);

%% Solve program and extract results
%set up problem
mset('yalmip',true);
mset(sdpsettings('solver', options.solver));

P = msdp(objective, ...
    mom_con, supp_con);

%solve LMI moment problem    

%for the flow problem, peak_test_alt has 1 substitution
%this script has 2 moment substitutions, so dual_rec does not have the same
%size as the vernoese map vv. I can't see a difference. What is going on?
[status,obj_rec, m,dual_rec]= msol(P);
% 
% M0 = double(mmat(mu0));
Mp = double(mmat(mup));
if nobj > 1
    Mc = double(mmat(muc));
end
%occupation measure
M_fw = cellfun(@(m) double(mmat(m)), mu_fw, 'UniformOutput', false);
M_bk = cellfun(@(m) double(mmat(m)), mu_bk, 'UniformOutput', false);

%can an approximately-optimal initial condition be recovered?
Mp_1 = Mp(1:(nvar+1), 1:(nvar+1));
rankp = rank(Mp_1, options.rank_tol);

xp_rec = double(mom(xp));

if INITIAL
    M0 = double(mmat(mu0));
    M0_1 = M0(1:(nvar+1), 1:(nvar+1));
    rank0 = rank(M0_1, options.rank_tol);
    x0_rec = double(mom(x0))/M0(1,1);
else
    M0 = [];
    x0_rec = [];
end


%how does this work with moment substitutions?
% %necessary with c^2 + s^2 = 1
% if options.scale
%     x0_rec = eval(xp_inv_scale, xp, x0_rec);
%     xp_rec = eval(xp_inv_scale, xp, xp_rec);
%     monp_unscale = subs(monp, xp, xp_inv_scale);
% else
%     monp_unscale = monp;
% end
% 

% %with equality constraints, dual_rec is bigger than just v
% %find the ordering of dual_rec, and which entries correspond to v
% %v = dual_rec'*monp_unscale;
dual_rec_v = dual_rec{1};

len_liou = length(Liou_fw);
dual_rec_v_fw = dual_rec_v(1:len_liou);
dual_rec_v_bk = dual_rec_v((len_liou+1):(2*len_liou));
 
if nobj == 1
    theta = 1;
else
    %this is an artifact of gloptipoly not having standard scalar
    %variables. A degree-1 measure is required.
    theta = dual_rec{2};
end

v_fw = dual_rec_v_fw'*monp;
v_bk = dual_rec_v_bk'*monp;


Lv_fw = [];
Lv_bk = [];
for i = 1:nsys
    if options.discrete
        Lv_fw_curr = subs(v_fw, [xp; wp], [options.dynamics.f{i}; wp]);
        Lv_bk_curr = subs(v_bk, [xp; wp], [options.dynamics.f{i}; wp]);
    else
        Lv_fw_curr = diff(v_fw, xp)*options.dynamics.f{i};
        Lv_bk_curr = diff(v_bk, xp)*options.dynamics.f{i};
    end
    Lv_fw = [Lv_fw; Lv_fw_curr];
    Lv_bk = [Lv_bk; Lv_bk_curr];
end

%% Output results to data structure

out = struct;

%recover optima
out.order = order;
out.peak_val = obj_rec;
out.x0 = x0_rec;
out.xp = xp_rec;
out.recover = (rankp == 1);
%moment matrices
out.M0 = M0;
out.Mp = Mp;
out.M_fw = M_fw;
out.M_bk = M_bk;

%objective
if nobj > 1
    out.Mc = double(mmat(muc));
end

if nw > 0
    out.w = double(mom(wp));
else
    out.w = [];
end

out.var = struct('x', xp, 'w', wp);

%Functions for evaluating system dynamics
out.func.fval = cell(nsys, 1);  %dynamics
out.func.Xval = cell(nsys, 1);  %support set
out.func.event = cell(nsys, 1); %Modification for ode45 event

for i = 1:nsys
    Xval_curr = @(x) all(eval([options.dynamics.X{i}; XR], xp, x));
    if nw > 0
        fval_curr = @(t, x, w) eval(options.dynamics.f{i}, [xp; wp], [x; w]);
        
        %space is inside support, time between Tmin and Tmax
        %if event_curr=0 stop integration, and approach from either
        %direction (only going to be in negative direction, given that
        %event is a 0/1 indicator function
        event_curr = @(t,x,w) support_event(t, x, Xval_curr, 0, Inf); 
        %should w be present here?
    else       
        fval_curr = @(t,x) eval(options.dynamics.f{i}, xp, x);
        
        event_curr = @(t, x) support_event(t, x, Xval_curr, 0, Inf);
    end
    
    out.func.fval{i} = fval_curr;
    out.func.Xval{i} = Xval_curr;        
    out.func.event{i} = event_curr;
end

out.dynamics = struct;
out.dynamics.f = out.func.fval;
out.dynamics.event = out.func.event;
out.dynamics.discrete = options.discrete;
out.dynamics.time_indep = 1;
event_all = @(tt, xt) cell2mat(cellfun(@(ex) ex(tt, xt), out.func.event, 'UniformOutput', false));
% %% functions and dual variables
out.func = struct;
out.func.dual_rec = dual_rec;
out.func.v_fw  = v_fw;
out.func.v_bk  = v_bk;
out.func.theta = theta;
out.func.Lv_fw = Lv_fw;
out.func.Lv_bk = Lv_bk;
 
% %functions that should be nonnegative along valid trajectories
out.func.cost_cell = cell(nobj, 1);
for i = 1:nobj
    out.func.cost_cell{i} = @(x) (eval(options.obj(i), xp, x));
end


out.func.cost_all = @(x) eval(options.obj, xp, x);

if nobj > 1
    out.func.cost = @(x) min(eval(options.obj, xp, x));
else
    out.func.cost = @(x) (eval(options.obj, xp, x));
end
% 
% %TODO: missing the 'theta' weights for cost function in this representation
% %under multiple cost functions. Check that out later, proper dual
% %representation and verification of nonnegativity
% 
if nw > 0    
    %some uncertain parameters
    out.func.vval_fw = @(x, w) eval(v_fw, [xp; wp], [x; w*ones(1, size(x, 2))]);     %forward value function
    out.func.vval_bk = @(x, w) eval(v_bk, [xp; wp], [x; w*ones(1, size(x, 2))]);     %backward value function
    
    out.func.vval = @(x, w) eval(v_fw + v_bk, [xp; wp], [x; w*ones(1, size(x, 2))]); %v_fw + v_bk (total value function)
    
    out.func.Lvval_fw = @(x) eval(Lv_fw, [xp; wp], [x; w*ones(1, size(x, 2))]);   %Forward Lie derivative or pushforward
    out.func.Lvval_bk = @(x) eval(Lv_bk, [xp; wp], [x; w*ones(1, size(x, 2))]);   %Backward Lie derivative or pushforward
    
    if options.discrete
         out.func.nonneg = @(x, w) [obj_rec - out.func.vval(x, w) - out.func.theta'*out.func.cost_all(x, w); ...
            options.discount*out.func.vval_fw(x, w) - out.func.Lvval_fw(x, w);...
            options.discount*out.func.vval_bk(x, w) + out.func.Lvval_bk(x, w)];
    else
        out.func.nonneg = @(x, w) [obj_rec - out.func.vval(x, w) - out.func.theta'*out.func.cost_all(x, w); ...
            out.func.vval_fw(x, w) - options.discount*out.func.Lvval_fw(x, w);...
            out.func.Lvval_bk(x, w) + options.discount*out.func.vval_bk(x, w)];
    end
else
    %no uncertain parameters
    out.func.vval_fw = @(x) eval(v_fw, xp, x);     %forward value function
    out.func.vval_bk = @(x) eval(v_bk, xp, x);     %backward value function
    
    out.func.vval = @(x) eval(v_fw + v_bk, xp, x); %v_fw + v_bk (total value function)
    
    out.func.Lvval_fw = @(x) eval(Lv_fw, xp, x);   %Forward Lie derivative or pushforward
    out.func.Lvval_bk = @(x) eval(Lv_bk, xp, x);   %Backward Lie derivative or pushforward
    
    if options.discrete
        %discrete
        out.func.nonneg = @(x) [obj_rec - out.func.vval(x) - out.func.theta'*out.func.cost_all(x); ...
            out.func.vval_fw(x) - options.discount*out.func.Lvval_fw(x);...
            out.func.Lvval_bk(x) + options.discount*out.func.vval_bk(x)];
    else
        %continuous
         out.func.nonneg = @(x) [obj_rec - out.func.vval(x) - out.func.theta'*out.func.cost_all(x); ...
            options.discount*out.func.vval_fw(x) - out.func.Lvval_fw(x);...
            options.discount*out.func.vval_bk(x) + out.func.Lvval_bk(x)];
    
    end
            
end




%
out.dynamics.cost = out.func.cost;
out.dynamics.nonneg = out.func.nonneg;
% %% Done!

end

function [mu, mon, X_occ, Ay] = occupation_measure(f, X, var, var_new, d, forward, discrete, discount)
%form the occupation measure
%Input:
%   f:          Dynamics (function)
%   X:          Support set upon which dynamics take place
%   var:        variables of problem
%   var_new:    New variables for forward occupation measure
%   forward:    forward occupation measure (1) or backwards (0)
%   d:          Degree of measure
%   discrete:   discrete system (true) or continuous system (false)
%
%Output:
%   mu:     Measure
%   mon:    monomoials
%   X:      Support of measure
%   Ay:     Equivalent of Liouville Equation

    if nargin < 6
        forward = 1;
    end

    if nargin < 7
        discrete = 0;
    end
    
    if nargin < 8
        discount = 1;
    end


    x_occ = var_new.x;
    
    if ~isempty(var_new.w)
        w_occ = var_new.w;
    else
        w_occ = [];
    end
    
    vars_prev = [var.x; var.w];
    vars_new = [x_occ; w_occ];
    
    f_occ = subs(f, vars_prev, vars_new);
   
    mu = meas([x_occ; w_occ]);
    mon = mmon([x_occ; w_occ], d);

    if discrete
        %discrete system
        pushforward = subs(mon, vars_new, [f_occ; w_occ]);
        if forward
            Ay = mom(mon) - discount * mom(pushforward);             
        else
            Ay = mom(pushforward) - discount * mom(mon);
        end
    else
        %continuous system
        lie_derivative = mom(diff(mon, x_occ)*f_occ);
        if forward
            Ay = discount * mom(mon) - lie_derivative;
        else
            Ay = discount * mom(mon) + lie_derivative;            
        end
    end
    %Ay = mom(diff(mon, x_occ)*f_occ);
    
    
    X_occ = subs(X, var.x, x_occ);

end

