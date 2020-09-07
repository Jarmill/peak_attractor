function [fig] = state_plot_3_cont(out, out_sim, out_sim_peak)
%STATE_PLOT_2 Plot states of system trajectories if x has 2 states
%out: information about the recovered solution
%out_sim: function evaluation on random sample trajectories
%out_sim_peak{1}: optional argument, function evaluation if peak is found at
%x0


Tmax = out_sim{1}.t(end);
nx = size(out_sim{1}.x, 2);

box_margin = 1.5;

assert(nx==3);

fig = figure;
clf
for i = 1:length(out_sim)
    t_curr = out_sim{i}.t;
    x_curr = out_sim{i}.x;
    if i == 1
    
        
        axis square
        hold on
        plot3(out_sim{i}.x(:, 1), out_sim{i}.x(:, 2), out_sim{i}.x(:, 3),  'c', 'DisplayName', 'Trajectories');
        title('Phase Plane of Trajectories')
        
        xlabel('x_1')
        ylabel('x_2')
        zlabel('x_3')
        
        legend('location', 'northwest') 
        
        if out.recover
            title(['Peak Value for Trajectories = ', num2str(out.peak_val, 5), ' (approx),  order = ', num2str(out.order)])
        else
            title(['Peak Value for Trajectories = ', num2str(out.peak_val, 5), ', order = ', num2str(out.order)])
        end
        
        
        
        
    %end
    else
        
        
        plot3(out_sim{i}.x(:, 1), out_sim{i}.x(:, 2), out_sim{i}.x(:, 3),  'c', 'HandleVisibility', 'Off');        
    end
end

%initial conditions
for i = 1:length(out_sim)
    if i == 1
       scatter3(out_sim{i}.x(1, 1), out_sim{i}.x(1, 2), out_sim{i}.x(1, 3), 100, 'k', 'DisplayName', 'Initial Points');       
    else        
       scatter3(out_sim{i}.x(1, 1), out_sim{i}.x(1, 2), out_sim{i}.x(1, 3), 100, 'k', 'HandleVisibility', 'Off');           
    end
end

%implicit curves
syms y [3 1]

% MD = 120;
MD = 80;

vy = out.func.vval(y);    
fimplicit3(vy, [stretch(xlim, box_margin), stretch(ylim, box_margin), stretch(zlim, box_margin)], 'EdgeColor', 'None','FaceColor', 'k', 'FaceAlpha', 0.3, ...
        'DisplayName', 'Invariant Set', 'MeshDensity', MD)
camlight; lighting phong;

%deal with multiple costs later

nobj = length(out.func.theta);
for i = 1:nobj % multiple objectives
    curr_cost = @(y) out.func.cost_cell{i}(y);
    cy = curr_cost(y) - out.peak_val;
    
    if i == 1
        fimplicit3(cy + 1e-8*sum(y), [xlim, ylim, zlim], 'EdgeColor', 'None','FaceColor', 'r', 'FaceAlpha', 0.3, ...
            'DisplayName', 'Cost Bound', 'MeshDensity', MD)

    else
        fimplicit3(cy + 1e-8*sum(y), [xlim, ylim, zlim], 'EdgeColor', 'None','FaceColor', 'r', 'FaceAlpha', 0.3, ...
            'HandleVisibility', 'Off', 'MeshDensity', MD)

    end
end

% cy = out.func.cost(y) - out.peak_val;
% 
% fimplicit3(cy + 1e-8*sum(y), [xlim, ylim, zlim], 'EdgeColor', 'None','FaceColor', 'r', 'FaceAlpha', 0.3, ...
%             'DisplayName', 'Cost Bound', 'MeshDensity', MD)


view(62, 17)

if out.recover && nargin == 3
    %plot the peak functions too
    npeak_traj = length(out_sim_peak);
    for k = 1:npeak_traj
        plot3(out_sim_peak{k}.x(:, 1), out_sim_peak{k}.x(:, 2), out_sim_peak{k}.x(:, 3), 'b',  'HandleVisibility', 'off', 'LineWidth', 2);
            
        if k == 1
%             plot3(out_sim_peak{1}.x(:, 1), out_sim_peak{1}.x(:, 2), out_sim_peak{1}.x(:, 3), 'b', 'DisplayName', 'Peak Traj.', 'LineWidth', 2);
            
            %initial condition
%             scatter3(out.x0(1, k), out.x0(2, k), out.x0(3, k), 200, 'ob', 'DisplayName', 'Peak Initial', 'LineWidth', 2);        
            scatter3(out.xp(1, k), out.xp(2, k), out.xp(3, k), 200, '*b', 'DisplayName', 'Peak Achieved', 'LineWidth', 2);        
        else
            
            
            %initial condition
%             scatter3(out.x0(1, k), out.x0(2, k), out.x0(3, k), 200, 'ob',  'HandleVisibility', 'off', 'LineWidth', 2);        
            scatter3(out.xp(1, k), out.xp(2, k), out.xp(3, k), 200, '*b',  'HandleVisibility', 'off', 'LineWidth', 2);        
        end
    end
end

end

