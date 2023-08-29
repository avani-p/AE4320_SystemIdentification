function PlotResidual(t,C1_res,label1, C2_res, label2, param_name,dataset_name)
    N=length(t);
    figure()
    sgtitle(append('Residual Analysis: Residual ',dataset_name))
    for i = 1:size(C1_res,2)
        subplot(2,3,i)
        grid on;
        plot(t,C1_res(:,i),'-s', Color = '#027148',MarkerIndices=100:300:N,LineWidth=1.5); hold on;
        plot(t,C2_res(:,i),'-^', Color = '#904180',MarkerIndices=200:300:N,LineWidth=1.5);
        yline(mean(C1_res(101:end,i)),'--', Color='#0000FF',LineWidth=1.2);
        title(append('Residual: ', param_name(i)))
        ylabel('Value');
        xlim([1 t(end)]);
        grid on
    end
    legend(label1,label2,'Mean Residual','Location','best')
    % 
    % figure()
    % sgtitle(append('Residual Analysis: Correlation ',dataset_name))
    % for i = 1:size(C1_res,2)
    %     subplot(2,3,i)
    %     grid on;
    %     plot(xcorr(C1_res(:,i)),'-', Color = '#027148',LineWidth=1.2); hold on;
    %     plot(xcorr(C2_res(:,i)),'-o', Color = '#904180',LineWidth=1.2);
    %     title(append('Auto-Correlation: ', param_name(i)))
    %     ylabel('Value');
    %     grid on
    % end
    % legend(label1,label2,'Location','best')
end