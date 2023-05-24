function [Y] = ObjectiveFcn(X)

    pho_x = X(1);
    pho_gr = X(2);
    g_swath = X(3);
    H = X(4);
    theta0 = X(5);
    
    [Lx, Lz, pho_h, G_ant_db, Ptx, N_coverage] = Ant_design2(pho_x, pho_gr, g_swath, H, theta0);
    
    Aa = Lx*Lz;
    p_threshold = 500;
    penalty = 1e10;
    
    if Ptx>p_threshold
        Y = penalty; 
    else
        Y = Aa;
    end
    
    
end
