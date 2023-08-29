function PlotVariances(var_CX, CX_param_names, var_CY, CY_param_names,var_CZ, CZ_param_names,...
    var_Cl, Cl_param_names, var_Cm, Cm_param_names, var_Cn, Cn_param_names, datasetCompare)
    
    if nargin<13
        datasetCompare = true;
    end

    figure(); 
    sgtitle('Parameter Variances');
    datasets = {'da3211_2.mat','dr3211_1.mat','de3211_1.mat'};
    for i = 1:size(var_CX,1)
        subplot(2,17,i);
        bar(var_CX(i,:)','r')
        grid on; hold on;
        title(CX_param_names(i),' ')
        hold off
        if datasetCompare
            set(gca,'xtick', 1:3,'xticklabel', datasets);
        end
    end

    for j = i+1:i+size(var_CY,1)
        subplot(2,17,j);
        bar(var_CY(j-i,:)','g')
        grid on; hold on;
        title(CY_param_names(j-i),' ')
        hold off
        if datasetCompare
            set(gca,'xtick', 1:3,'xticklabel', datasets);
        end
    end

    for i = j+1:j+size(var_CZ,1)
        subplot(2,17,i);
        bar(var_CZ(i-j,:)','b')
        grid on; hold on;
        title(CZ_param_names(i-j),' ')
        hold off
        if datasetCompare
            set(gca,'xtick', 1:3,'xticklabel', datasets);
        end
    end

    
    for j = i+1:i+size(var_Cl,1)
        subplot(2,17,j);
        bar(var_Cl(j-i,:)','c')
        grid on; hold on;
        title(Cl_param_names(j-i),' ')
        hold off
        if datasetCompare
            set(gca,'xtick', 1:3,'xticklabel', datasets);
        end
    end

    for i = j+1:j+size(var_Cm,1)
        subplot(2,17,i);
        bar(var_Cm(i-j,:)','m')
        grid on; hold on;
        title(Cm_param_names(i-j),' ')
        hold off
        if datasetCompare
            set(gca,'xtick', 1:3,'xticklabel', datasets);
        end
    end

    for j = i+1:i+size(var_Cn,1)
        subplot(2,17,j);
        bar(var_Cn(j-i,:)','k')
        grid on; hold on;
        title(Cn_param_names(j-i),' ')
        hold off
        if datasetCompare
            set(gca,'xtick', 1:3,'xticklabel', datasets);
        end
    end

end

