function display_sol(s, u)
    n = 6;
    m = 3;
    n_steps = length(s)/6;
    r = zeros(n_steps, 3);
    v = zeros(n_steps, 3);
    for i = 1:n_steps
         r(i,:) = s((i-1)*n+1:(i-1)*n+3);
         v(i,:) = s((i-1)*n+4:(i-1)*n+6);
    end
    figure
    subplot(3,1,1)
    hold on
    plot(r(:,1))
    plot(r(:,2))
    plot(r(:,3))
    hold off
    grid on
    subplot(3,1,2)
    hold on
    plot(v(:,1))
    plot(v(:,2))
    plot(v(:,3))
    hold off
    grid on
    subplot(3,1,3)
    hold on
    plot(u(1:m:end))
    plot(u(2:m:end))
    plot(u(3:m:end))
    hold off
    grid on
    
    figure
    plot3(r(:,1), r(:,2), r(:,3))
    hold on
    plot3(0,0,0,'*')
    axis equal
end