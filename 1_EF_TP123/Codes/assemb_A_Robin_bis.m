function A = assemb_A_Robin_bis(K, alpha, mesh)

A = sparse(mesh.nbs,mesh.nbs); 


% Boucle sur les triangles
for ie = 1:mesh.nbt
    
    is = mesh.elm_som(ie,:);
   
    x = mesh.som_coo(is,1);
    y = mesh.som_coo(is,2);
    
    a = circshift(mesh.som_coo(is,2),2)-circshift(mesh.som_coo(is,2),1);
    b = circshift(mesh.som_coo(is,1),1)-circshift(mesh.som_coo(is,1),2);
    
    aire = 0.5*abs(det([mesh.som_coo(is,:), ones(3,1)]));
    
    Alm = (K(ie)/(4*aire))*(a*a' + b*b');
    
    A(is,is) = A(is,is) + Alm;
    
   
end

% Boucle sur les segements du bord
for ie = 1:mesh.nbab
   
    is = mesh.abd_som(ie, :);
    
    x = mesh.som_coo(is,1);
    y = mesh.som_coo(is,2);
    
    d = sqrt((x(1)-x(2)).^2 + (y(1)-y(2)).^2);
    A(is,is) = A(is,is) + (alpha*d/6)*(ones(2,2)+eye(2,2));
   
end


