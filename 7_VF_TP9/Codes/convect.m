function ynew = convect( mesh, u, yold , Gammam, g, dt )


deja_traite = zeros(mesh.nbs);

% Initialisation de ynew
for i = Gammam 
    xi = mesh.som_coo(i,1); 
    yi = mesh.som_coo(i,2); 
    ynew(i) = g(xi, yi);
    deja_traite(i) = 1; 
end 


% Boucle sur les triangles K
for ie = 1:mesh.nbt
    
    is = mesh.elm_som(ie, :);
    
    xi = mesh.som_coo(is(:), 1);
    yi = mesh.som_coo(is(:), 2);
    
    
    A = [ xi(1) xi(2) xi(3);  yi(1) yi(2) yi(3); 1 1 1]; 

    for i = 1:3
        if (deja_traite(is(i)) < 0.5)
            
            % Second membre pour les fonctions barycentriques
            B = [ xi(i) - dt*u(is(i), 1) ;  yi(i) - dt*u(is(i), 2) ; 1 ];
            lambda = A\B;
            
            % Si on est dans le triangle alors on met a jour
            if (min(lambda) >= 0) 
                ynew(is(i)) = lambda(1)*yold(is(1)) + ...
                              lambda(2)*yold(is(2)) + ...
                              lambda(3)*yold(is(3));
                deja_traite(is(i)) = 1;
            end 
        end
    end
end




end

