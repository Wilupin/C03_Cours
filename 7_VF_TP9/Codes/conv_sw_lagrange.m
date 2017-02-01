function Vnew_t = conv_sw_lagrange(mesh, Vold_t, dt)

% Procede a une iteration par la methode des volumes finis

g = 9.81;

% ---- Normales aux arretes
N1_KL = mesh.fac_nor(:,1);
N2_KL = mesh.fac_nor(:,2);

dsol_t = zeros(3, mesh.nbt);



% ---- Assemblage de la matrice A
for ia=1:mesh.nba

    ie = mesh.fac_elm(ia, :); 
    nu = [N1_KL(ia); N2_KL(ia)];
    
    
    % ---- Flux pour les arretes interieures
    if(mesh.fac_zon(ia) == 0) 
              
        [ pstar, unustar, UKL] = varstar( nu, [Vold_t(:,ie(1)),  Vold_t(:, ie(2))] );

        phi = UKL + [ 0; pstar*nu(1); pstar*nu(2)];  
        
        phi = phi*mesh.fac_mes(ia);
        
        dsol_t(:,ie(1)) = dsol_t(:,ie(1)) + phi;
        dsol_t(:,ie(2)) = dsol_t(:,ie(2)) - phi;
        
    end 
    
    
    % ----Flux pour les arretes du mur
    if(mesh.fac_zon(ia) == 1)
        hp = Vold_t(1,ie(1)); 
    
        phi = [ 0 ; ...
            0.5*g*(hp^2)*nu(1); ...
            0.5*g*(hp^2)*nu(2)];
        
        phi = phi*mesh.fac_mes(ia);
        
        dsol_t(:,ie(1)) = dsol_t(:,ie(1)) + phi;
        
    end
    
    % ---- Flux pour les arretes du bord artificiel
    if(mesh.fac_zon(ia) == 2)
        
        ps_vnu =  Vold_t(2:3,ie(1))'*nu;
        hp = Vold_t(1,ie(1)); 
        
        phi = [ ps_vnu; ...
            (Vold_t(2, ie(1))/hp)*ps_vnu + 0.5*g*(hp^2)*nu(1); ...
            (Vold_t(3, ie(1))/hp)*ps_vnu + 0.5*g*(hp^2)*nu(2)];   
        
        phi = phi*mesh.fac_mes(ia);
        
        dsol_t(:,ie(1)) = dsol_t(:,ie(1)) + phi;
        
    end
    
    
end

mes = [mesh.elm_mes'; mesh.elm_mes'; mesh.elm_mes'];

Vnew_t = Vold_t - dt*dsol_t./mes;



end
