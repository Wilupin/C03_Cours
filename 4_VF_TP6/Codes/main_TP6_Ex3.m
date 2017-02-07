% Script de test de la fonction tri_to_sum
%
% Gaetan Facchinetti
%
% Cours C03

clear all; 


% ---- Lecture du maillage et enrichissement volumes finis
mesh = lect_mesh('../Meshs/disq0');
mesh = face_number(mesh);


% ---- Defintion de la fonction de test
f = @(x,y) (x.*y);


% ---- Calcul de f aux centres des triangles du maillage
sol_t = f(mesh.elm_gra(:,1), mesh.elm_gra (:,2));


% Utilisation de tri_to_sum
sol_s = tri_to_sum(mesh, sol_t);

% ---- Coordonnees des triangles
tri = mesh.elm_som;
x = mesh.som_coo(:,1); 
y = mesh.som_coo(:,2); 


% ---- Affichage de la solution
figure(1)
trisurf(tri, x, y, f(x,y));
figure(2)
trisurf(tri, x, y, sol_s); 
figure(3)
trisurf(tri, x, y, f(x,y) - sol_s);


% NB : On trouve une erreur plus importante aux bords puisque l'on a une
% formule tri_to_sum qui est decentree et donc gere plus mal les sommets
% aux bords. 

% ---- Erreur en norme 2 et en norme infinie
norm_2   = sqrt((f(x,y)-sol_s)'*(f(x,y)-sol_s));
norm_inf = max(f(x,y)-sol_s);
