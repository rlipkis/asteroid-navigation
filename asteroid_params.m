function p = asteroid_params()
    % asteroid parameters
    rho = 2000; % [kg/m^3]
    R1 = 2000; % [m]
    R2 = 1000; % [m]
    axis = [1; 0; 0];
    
    % asteroid masses, positions (COM at origin)
    G = 6.674e-11; 
    M1 = (4*pi/3)*R1^3*rho;
    M2 = (4*pi/3)*R2^3*rho;
    r1 = R1*axis;
    r2 = -R2*axis;
    r_cm = (M1*r1 + M2*r2)/(M1 + M2);
    r1 = r1 - r_cm;
    r2 = r2 - r_cm;
    
    % build struct
    p = struct;
    p.mu1 = G*M1;
    p.mu2 = G*M2;
    p.r1 = r1;
    p.r2 = r2;
    p.R1 = R1;
    p.R2 = R2;
end