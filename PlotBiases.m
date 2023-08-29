function PlotBiases(t, X1, X1_label, X2, X2_label)
    N = length(t);
    state_order = {'\lambda_{x_r}', '\lambda_{y_r}', '\lambda_{z_r}', ...
        '\lambda_{p_r}', '\lambda_{q_r}', '\lambda_{r_r}'}';
    
    state_units = {'m/s^2', 'm/s^2', 'm/s^2', 'rad', 'rad', 'rad'}';
    
    state_title = {'Accelerometer bias \lambda_{x_r}',...
        'Accelerometer bias \lambda_{y_r}', 'Accelerometer bias \lambda_{z_r}',...
        'Rate gyro bias \lambda_{p_r}', 'Rate gyro bias \lambda_{q_r}', 'Rate gyro bias \lambda_{r_r}'}';
    
    ylabel_vector = append(state_order, {' ['},state_units, {']'}); 
    legend_vector_1 = append(X1_label, state_order);
    legend_vector_2 = append(X2_label, state_order);
    range = [0.025 0.025 0.025 1e-4 1e-4 1e-4];
    true_val  = [0.02  0.02  0.02  deg2rad(0.003) deg2rad(0.003) deg2rad(0.003)];
    
    for i = 1:size(X1, 1)
        if(rem(i,6)==1); figure(); sgtitle('Biases Comparison'); end
        if(rem(i,6)==0); subplot(2,3,6); else; subplot(2,3,rem(i,6)); end
        plot(t, X2(i, :),LineWidth=1.50,Marker='s',MarkerIndices=1:500:N)
        hold on ; grid on
        plot(t, X1(i, :),LineWidth=1.50,Marker='o',MarkerIndices=200:500:N)
        yline(true_val(i),"b--",LineWidth=1.2)
        xlabel('Time [s]')
        ylabel(ylabel_vector(i))
        ylim([-(1e-4) range(i)])
        title( state_title(i))
        legend(string(legend_vector_2(i)), string(legend_vector_1(i)),Location='best')
        hold off
    end
end

