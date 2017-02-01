function sol_s = tri_to_som( mesh, sol_t )

% Entree : sol_t (solution aux triangles) size : mesh.nbt x 1
% Sortie : sol_s (solution aux sommets)   size : mesh.nbs x 1

% Initialisation de tableaux
sol_s = zeros(mesh.nbs,1); 
num   = zeros(mesh.nbs,1); 
den   = zeros(mesh.nbs,1);

% Boucle sur les triangles
for ie=1:mesh.nbt
   
    som_K = mesh.elm_som(ie,:); 
    
    num(som_K) = num(som_K) + mesh.elm_mes(ie).*sol_t(ie).*ones(3,1); 
    den(som_K) = den(som_K) + mesh.elm_mes(ie).*ones(3,1);
end

sol_s =  num./den; 

end

