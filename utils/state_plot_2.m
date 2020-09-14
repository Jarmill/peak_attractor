function [fig] = state_plot_2_cont(out, out_sim, out_sim_peak)
%STATE_PLOT_2 Plot states of system trajectories if x has 2 states
%out: information about the recovered solution
%out_sim: function evaluation on random sample trajectories
%out_sim_peak{1}: optional argument, function evaluation if peak is found at
%x0


Tmax = out_sim{1}.t(end);
nx = size(out_sim{1}.x, 2);

box_margin = 1.2;


assert(nx==2);

fig = figure;
clf
for i = 1:length(out_sim)
    t_curr = out_sim{i}.t;
    x_curr = out_sim{i}.x;
    if i == 1  
        axis square
        hold on
        if out.dynamics.discrete
            scatter(out_sim{i}.x(:, 1), out_sim{i}.x(:, 2), '.c', 'DisplayName', 'Trajectories');
        else
            plot(out_sim{i}.x(:, 1), out_sim{i}.x(:, 2), 'c', 'DisplayName', 'Trajectories');
        end
%         title(['Phase Plane of Trajectories, order = ', num2str(out.order)])
        
        xlabel('x_1')
        ylabel('x_2')
        
        legend('location', 'northwest') 

        %title formatting
        peak_str = ['Peak Value for Trajectories = ', num2str(out.peak_val, 4), ' order = ', num2str(out.order)];       
        title(peak_str)

        
        
        legend('location', 'northwest') 
        
    %end
    else 
        if out.dynamics.discrete
            scatter(out_sim{i}.x(:, 1), out_sim{i}.x(:, 2), '.c', 'HandleVisibility', 'off');
        else
            plot(out_sim{i}.x(:, 1), out_sim{i}.x(:, 2), 'c', 'HandleVisibility', 'off');
        end
    end
end


for i = 1:length(out_sim)
    if i == 1    
        if ~isempty(out.x0)
            scatter(out.x0(1), out.x0(2), 100, 'k', 'DisplayName', 'Initial Points');        
        else
            scatter(out_sim{i}.x(1, 1), out_sim{i}.x(1, 2), 100, 'k', 'DisplayName', 'Initial Points');                    
        end

    else
       
       scatter(out_sim{i}.x(1, 1), out_sim{i}.x(1, 2), 100, 'k', 'HandleVisibility', 'off');        
    
    end
end

%implicit curves
syms y [2 1]

xlim(stretch(xlim, box_margin))
ylim(stretch(ylim, box_margin))

if isempty(out.var.w)
    vy = out.func.vval(y);

    fimplicit(vy, [xlim, ylim], ':k', 'DisplayName', 'Invariant Set', 'LineWidth', 3);
end

nobj = length(out.func.theta);
for i = 1:nobj % multiple objectives
    curr_cost = @(y) out.func.cost_cell{i}(y);
    cy = curr_cost(y) - out.peak_val;
    
    if i == 1
        fimplicit(cy + 1e-8*sum(y), [xlim, ylim],  '--r', 'DisplayName', 'Cost Bound', 'LineWidth', 2)
    else
        fimplicit(cy + 1e-8*sum(y), [xlim, ylim],  '--r', 'HandleVisibility', 'Off', 'LineWidth', 2)
    end
end
% cy = out.func.cost(y) - out.peak_val;
% fimplicit(cy + 1e-8*sum(y), [xlim, ylim],  '--r', 'DisplayName', 'Cost Bound', 'LineWidth', 2);


if out.recover && (nargin == 3)
    %plot the peak functions too
    npeak_traj = length(out_sim_peak);
    for k = 1:npeak_traj
        if k == 1
            if out.dynamics.discrete
                scatter(out_sim_peak{k}.x(:, 1), out_sim_peak{k}.x(:, 2), '.b', 'HandleVisibility', 'off');
            else
                plot(out_sim_peak{k}.x(:, 1), out_sim_peak{k}.x(:, 2), 'b', 'HandleVisibility', 'off');
            end
            scatter(out.xp(1, k), out.xp(2, k), 200, '*b', 'DisplayName', 'Peak Achieved', 'LineWidth', 2);        

        else
            if out.dynamics.discrete
                scatter(out_sim_peak{k}.x(:, 1), out_sim_peak{k}.x(:, 2), '.b', 'HandleVisibility');
            else
                plot(out_sim_peak{k}.x(:, 1), out_sim_peak{k}.x(:, 2), 'b', 'HandleVisibility', 'off');
            end
            scatter(out.xp(1, k), out.xp(2, k), 200, '*b', 'HandleVisibility', 'off','LineWidth', 2);        
       end
    end
             
end

end

