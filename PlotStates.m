function PlotStates(type, t, X1, X1_label, X2, X2_label, X3, X3_label)
    if nargin < 7
        dont_skip_3 = 0;
        X3 = 1;
        X3_label = 1;
    else 
        dont_skip_3 = 1;
    end
    N = length(t);
    states = 0;
    if (type(1:6) == 'States')
        states = 1;
    end

    if states ==1
        state_order = {'x', 'y', 'z', 'u_{air}^B', 'v_{air}^B', 'w_{air}^B',... 
            '\phi', '\theta', '\psi', 'W_{x_E}', 'W_{y_E}', 'W_{z_E}', '\lambda_{x_r}',...
            '\lambda_{y_r}', '\lambda_{z_r}', '\lambda_{p_r}', '\lambda_{q_r}', '\lambda_{r_r}'}';
        
        state_units = {'m', 'm', 'm', 'm/s', 'm/s', 'm/s', 'rad', 'rad', 'rad',...
            'm/s', 'm/s', 'm/s', 'm/s^2', 'm/s^2', 'm/s^2', 'rad', 'rad', 'rad'}';
        
        state_title = {'Position x', 'Position y', 'Position z', 'Velocity u',...
            'Velocity v', 'Velocity w', 'Euler angle \phi', 'Euler angle \theta',...
            'Euler angle \psi', 'Wind velocity W_{x_E}', 'Wind velocity W_{y_E}',...
            'Wind velocity W_{z_E}', 'Accelerometer bias \lambda_{x_r}',...
            'Accelerometer bias \lambda_{y_r}', 'Accelerometer bias \lambda_{z_r}',...
            'Rate gyro bias \lambda_{p_r}', 'Rate gyro bias \lambda_{q_r}', 'Rate gyro bias \lambda_{r_r}'}';
    else
        state_order = {'x', 'y', 'z', 'u_{air}', 'v_{air}', 'w_{air}',... 
            '\phi', '\theta', '\psi', 'V', '\alpha', '\beta'}';
        state_units = {'m', 'm', 'm', 'm/s', 'm/s', 'm/s', 'rad', 'rad', 'rad',...
            'm/s', 'rad', 'rad'}';
        state_title = {'Position x', 'Position y', 'Position z', 'Velocity u',...
            'Velocity v', 'Velocity w', 'Euler angle \phi', 'Euler angle \theta',...
            'Euler angle \psi', 'Airspeed V', 'Angle of Attack \alpha', 'Sideslip Angle \beta'}';
    end
    
    ylabel_vector = append(state_order, {' ['},state_units, {']'}); 
    if(dont_skip_3)
        legend_vector_1 = append(state_order, X3_label);
    end
    legend_vector_2 = append(state_order, X2_label);
    legend_vector_3 = append(state_order, X1_label);
    
    for i = 1:size(X1, 1)
        if(rem(i,6)==1); figure(); sgtitle(append(type,' Value Comparison')); end
        if(rem(i,6)==0); subplot(2,3,6); else; subplot(2,3,rem(i,6)); end
        if (dont_skip_3)
            plot(t, X3(i, :),Color = '#FF0000', LineWidth=1.50,Marker='s',MarkerIndices=1:500:N)
        end
        hold on ; grid on
        plot(t, X2(i, :), Color = '#0000FF', LineWidth=1.50,Marker='*',MarkerIndices=175:500:N)
        plot(t, X1(i, :), Color = '#77AC30', LineWidth=1.50,Marker='o',MarkerIndices=325:500:N)
        xlabel('Time [s]')
        ylabel(ylabel_vector(i))
        title(state_title(i))
        if(dont_skip_3)
            legend(string(legend_vector_1(i)), string(legend_vector_2(i)), ...
            string(legend_vector_3(i)),Location='best')
        else
            legend(string(legend_vector_2(i)), ...
            string(legend_vector_3(i)),Location='best')
        end
        hold off
    end
end

