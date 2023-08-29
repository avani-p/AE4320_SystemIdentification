function PlotAerodynamicModel(t, data_C_true, data_C_KF, data_C_OLS, coeff_name, val_dataset,secondStruct)
    
    N = size(data_C_true,2);
    figure();
    for i = 1:N
        subplot(2,3,i)
        plot(t, data_C_KF(:,i), '-s',  Color = '#0000FF', MarkerIndices=100:300:size(data_C_true,1),LineWidth=1.5)
        hold on;
        plot(t, data_C_OLS(:,i),'-o', Color = '#FF0000', MarkerIndices=200:300:size(data_C_true,1),LineWidth=1.5)
        plot(t, data_C_true(:,i), Color = '#77AC30', LineWidth=1.2, Marker='*',MarkerIndices=1:300:size(data_C_true,1))
        hold off; xlim([1 t(end)]);
        title(coeff_name(i))
    end
    if (~secondStruct)
        sgtitle(append('Aerodynamic Model Parameters: ',val_dataset));
        legend('Dimensionless','OLS Est.','True')
    else
        sgtitle(append('Aerodynamic Model Parameters (New Model Validation) : ',val_dataset));
        legend('OLS Est.','OLS Est. 2nd Model', 'True')
    end
end