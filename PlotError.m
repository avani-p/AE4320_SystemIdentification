function PlotError(type, t, EstErr_x)
    N = length(t);
    states = 0;
    if (type(1:6) == 'Conver')
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
   for i = 1:size(EstErr_x,1)
        if(rem(i,6)==1); figure(); sgtitle(append('Estimation Error: ',type)); end
        if(rem(i,6)==0); subplot(2,3,6); else; subplot(2,3,rem(i,6)); end
        plot(t, EstErr_x(i,:),LineWidth=1.50,Marker='s',MarkerIndices=1:500:N)
        hold on; grid on
        yline(mean(EstErr_x(i,:)),'--r',LineWidth=1.5)
        xlabel('Time [s]')
        ylabel(ylabel_vector(i))
        title(state_title(i))
        hold off
        if(rem(i,6)==0); legend('Error','Mean Error'); end
   end
end

