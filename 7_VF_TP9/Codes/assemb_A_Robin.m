function A = assemb_A_Robin(K, alpha, mesh)


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
    
    id = find(mesh.som_zon(is) ~= 0);
    
    if (max(size(id)) == 2)
        d = sqrt((x(id(1))-x(id(2))).^2 + (y(id(1))-y(id(2))).^2);
        Alm_bord = ((alpha*d)/3)*(0.5*ones(2,2) + 0.5*eye(2,2));
        A(is(id),is(id)) = A(is(id), is(id)) + Alm_bord;
    end
    
    
end


