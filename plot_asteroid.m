function plot_asteroid()
    p = asteroid_params();
    
    [X,Y,Z] = sphere(100);
    C1 = rand(size(X));
    C2 = rand(size(X));

    surf(p.r1(1) + p.R1*X, p.r1(2) + p.R1*Y, p.r1(3) + p.R1*Z, 0.5 + 0.1*C1)
    surf(p.r2(1) + p.R2*X, p.r2(2) + p.R2*Y, p.r2(3) + p.R2*Z, 0.5 + 0.1*C2)
%     alpha 0.75
    colormap gray
    caxis([0 1])
    shading interp
    axis equal
    set(gca,'Color','k')
end