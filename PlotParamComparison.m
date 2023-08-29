function new_params = PlotParamComparison(OLS_params,RLS_params, param_names,deleteOutliers,showPlot)
    
    if all([(size(OLS_params)~=size(RLS_params)) (showPlot)])
        error('Size of the Parameter Arrays do not match.')
    end


    N = size(OLS_params,1);
    new_params = zeros(size(OLS_params));
    
    figure_handle = figure();
    grid on;
    xlim([0 size(OLS_params,2)+1])
    sgtitle('Parameter Comparison of OLS and RLS Estimates')
    for i = 1:N
        OLS = OLS_params(i,:);
        RLS = RLS_params(i,:);
        new_params(i,:) = OLS;
        if (deleteOutliers)
            idx = unique([find(abs(OLS)>1e3) find(abs(OLS)<1e-9)]);
            OLS(idx) = [];
            RLS(idx) = [];
            for k = 1:length(idx)
                new_params(i,idx(k)) = 0;
            end
        end
        if(showPlot)
            subplot(2,3,i)
            boxplot([OLS' RLS'],'Labels',{'OLS','RLS'}); 
            hold on; yline(mean(OLS),'--magenta');
            % x_vals = ones(length(OLS),1);
            % plot(x_vals,OLS ,'rd',2*x_vals, RLS,"b+");
            % set(gca,'xtick', 1:N ,'xticklabel', {'OLS', 'RLS',''});
            title(param_names(i))
            ylabel('Value'); xlim([0 3])
            grid on
            if(i==N); legend('Mean','Location','best'); end
        end
    end
    if(~showPlot); close(figure_handle); end
end