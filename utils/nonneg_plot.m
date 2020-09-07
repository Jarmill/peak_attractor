function [fig] = nonneg_plot(out, out_sim, out_sim_peak)
%NONNEG_PLOT Plot nonnegative functions (value, Lie derivatives, cost
%comparision)
%out:   output from peak estimation routine
%out_sim: function evaluation on random sample trajectories
%out_sim_peak: optional argument, function evaluation if peak is found at
%x0


%if iscell(f)
%    nsys = length(f);
%else
%    nsys = 1;
%end

nplt = size(out_sim{1}.nonneg, 1);

fig = figure;
clf
for i = 1:length(out_sim)
    t_curr = out_sim{i}.t;
    nonneg_curr = out_sim{i}.nonneg;
    for k = 1:(nplt)
        subplot(nplt, 1, k)
        hold on
        
        
        if i == 1
            %change this to allow for switching
            if k == 1
                title('Cost Comparision')
                xlabel('time')
                ylabel('$p^* -\theta^T cost(x) - v_+(x) - v_-(x)$', 'Interpreter', 'Latex')
                legend('location', 'east')
                plot(t_curr, nonneg_curr(k, :), 'c', 'DisplayName', 'Trajectories')
            elseif k == 2    
                plot(t_curr, nonneg_curr(k, :), 'c','DisplayName', 'Trajectories')
                title('Forward Discounted Change in Value Function')
                xlabel('time')
                if out.dynamics.discrete
                    ylabel('$v_+ (x) - \alpha v_+ \circ f(x)$', 'Interpreter', 'Latex')                
                else
                    ylabel('$\beta v_+ (x) - L_f v_+(x)$', 'Interpreter', 'Latex')                
                end
                legend('location', 'east')
            else
%                 title(['Forward Change in Value Function', num2str(k-1)])
                title('Backward Discounted Change in Value Function')
                xlabel('time')
%                 ylabel(['L_{f', num2str(k-1), '} v(t,x)'])
                if out.dynamics.discrete                
                    ylabel('$v_- \circ f(x) + \alpha v_-(x)$', 'Interpreter', 'Latex')
                else
                    ylabel('$\beta v_- (x) + L_f v_-(x)$', 'Interpreter', 'Latex')
                end
                legend('location', 'east')
                plot(t_curr, nonneg_curr(k, :), 'c', 'DisplayName', 'Trajectories')
            end
            
        else
            plot(t_curr, nonneg_curr(k, :), 'c','HandleVisibility', 'off')                
        end
        if i == length(out_sim)
            plot(xlim, [0, 0], ':k', 'HandleVisibility', 'off')
        end
    end       
end

if nargin == 3  
    for k = 1:(nplt)
        subplot(nplt, 1, k)
        hold on
        plot(out_sim_peak{1}.t, out_sim_peak{1}.nonneg(k, :), 'b', 'DisplayName', 'Peak Traj.')           
    end       
    
%     if ~out.dynamics.time_indep
        if ~isempty(out.w)
            peak_nonneg = out.func.nonneg(out.xp, out.w);
        else
            peak_nonneg = out.func.nonneg(out.xp);
        end
        for k = 1:nplt
            subplot(nplt, 1, k)
            scatter(out_sim_peak{1}.tp, peak_nonneg(k), 300, '*b', 'Linewidth', 2, 'HandleVisibility', 'Off')
        end
%     end 
end

end

