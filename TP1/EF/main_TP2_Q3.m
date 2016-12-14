clear all;

% Lecture du mesh
mesh = lect_mesh('car0');

% choix de kappa
kappa = ones(mesh.nbt,1);

% Assemblage de la matrice de rigidité
A = assemb_A(kappa, mesh);
spy(A);

% Assemblage du second membre
F = assemb_F(@(x,y) fun(x,y), mesh); 
M = assemb_M(mesh);

% Initialisation de l'inconnue
u = zeros(mesh.nbs,1);

% Recuperation des donnees au bord
%dir = find(mesh.som_zon ~= 0);
dir = find(mesh.som_zon == 1 | mesh.som_zon == 2);
inconnues = setdiff(1:mesh.nbs, dir);
u(dir) = function_g(mesh.som_zon(dir),mesh.som_coo(dir,1), mesh.som_coo(dir,2));

% Autre façon de dire 
% g = @(z,x,y) ((z==1)*(0.0) +(z==2)*(2.0))
% u(dir) = g(mesh.som_zon(dir),mesh.som_coo(dir,1), mesh.som_coo(dir,2))
% inconnues = setdiff(1:mesh.nbs, dir);

% Pseudo elimination
F = F-A*u;

% Resolution du système lineaire
u(inconnues) = A(inconnues, inconnues)\F(inconnues); 


tri = mesh.elm_som;
x = mesh.som_coo(:,1);
y = mesh.som_coo(:,2);

trimesh(tri, x, y, u);


% Verification de l'assemblage de M 
% Res doit valoir l'aire du domaine
v = ones(mesh.nbs,1);
Res = v'*M*v;
