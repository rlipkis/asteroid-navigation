function display_sol(s, u, dt)
    n = 6;
    m = 3;
    n_steps = length(s)/6;
    ts = (1:n_steps)*dt;
    
    r = zeros(n_steps, 3);
    v = zeros(n_steps, 3);
    for i = 1:n_steps
         r(i,:) = s((i-1)*n+1:(i-1)*n+3);
         v(i,:) = s((i-1)*n+4:(i-1)*n+6);
    end
    figure
    subplot(3,1,1)
    hold on
    plot(ts, r(:,1))
    plot(ts, r(:,2))
    plot(ts, r(:,3))
    hold off
    grid on
    ylabel('Position [m]')
    
    subplot(3,1,2)
    hold on
    plot(ts, v(:,1))
    plot(ts, v(:,2))
    plot(ts, v(:,3))
    hold off
    grid on
    ylabel('Velocity [m/s]')
    
    subplot(3,1,3)
    hold on
    plot(ts, u(1:m:end))
    plot(ts, u(2:m:end))
    plot(ts, u(3:m:end))
    hold off
    grid on
    xlabel('Time [s]')
    ylabel('Control input [N]')
    
    figure
    hold on
    plot_asteroid()
    plot3(r(:,1), r(:,2), r(:,3), 'r', 'linewidth', 2)
    axis equal
    xlabel('x [m]')
    ylabel('y [m]')
    zlabel('z [m]')
    view([16 11])
end