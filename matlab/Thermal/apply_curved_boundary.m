function g_new_val = apply_curved_boundary(i, j, ksi, ...
    Ksi, w, g, g_eq, T, Rho, Zone_ID, ...
    R, Tau, T_w, radius, x, y, x_circ, y_circ)
    
    % Get neighbor indices
    i_n = i + Ksi(1, ksi);
    j_n = j + Ksi(2, ksi);
    
    % Get opposite direction
    opp = [1, 4, 5, 2, 3, 8, 9, 6, 7];
    opp_k = opp(ksi);
    
    % Check if neighbor is a solid boundary
    if Zone_ID(j_n, i_n) == 1
        % Fluid node
        x2 = x(i);
        y2 = y(j);
        
        % Neighbor node (approach point)
        x1 = x(i_n);
        y1 = y(j_n);
        
        % Find wall point on curved surface
        C_w = find_the_wall_point(x1, y1, x2, y2, radius, x_circ, y_circ);
        
        % Compute delta (fractional distance to wall)
        delta = sqrt((C_w(1)-x2)^2 + (C_w(2)-y2)^2) / sqrt((x1-x2)^2 + (y1-y2)^2);
        
        % Macroscopic temperatures
        T_f  = T(j,i);              % Fluid node
        T_ff = T(j - Ksi(2,ksi), i - Ksi(1,ksi)); % Opposite direction (2nd neighbor)
        
        % Interpolate wall temperature
        Tb1 = ((delta - 1)*T_f + T_w) / delta;
        Tb2 = ((delta - 1)*T_ff + 2*T_w) / delta;
        
        if delta >= 0.75
            Tb = Tb1;
        else
            Tb = delta*Tb1 + (1-delta)*Tb2;
        end
        
        % Equilibrium and non-equilibrium components
        g_eq_alpha = w(ksi) * Rho(j,i) * R * Tb;
        g_non_eq = g(j,i,ksi) - g_eq(j,i,ksi);
        if delta < 0.75
            g_non_eq = delta * g_non_eq + (1 - delta) * g_non_eq;
        end
        
        % Final reflected PDF
        g_new_val = g_eq_alpha + (1 - 1/Tau) * g_non_eq;
    else
        % Simple streaming if no boundary
        g_new_val = g(j_n, i_n, opp_k);
    end
end