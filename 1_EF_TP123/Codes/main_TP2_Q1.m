% ------------------------------
%
% Ga?tan Facchinetti
%
% Cours C03 : Methodes num?riques
% Resolution de l'equation de la chaleur par EF. 
% Fiche TP1 et TP2 
%
% ------------------------------

clear all;


% ---- Lecture du maillage
mesh = lect_mesh('../Meshs/DOM1');


% ---- Choix de kappa dans l'?quation
kappa = ones(mesh.nbt,1);


% ---- Assemblage des matrices et vecteurs
A = assemb_A(kappa, mesh);                % Matrice de rigidite
F = assemb_F(@(x,y) fun(x,y), mesh);      % Vecteur second membre


% ---- Initialisation du vecteur inconnu
u = zeros(mesh.nbs,1);


% ---- Recuperation des donnees au bord
dir = find(mesh.som_zon == 1 | mesh.som_zon == 2);
inconnues = setdiff(1:mesh.nbs, dir);
g = @(z,x,y) ((z==1)*(0.0) +(z==2)*(2.0));
u(dir) = g(mesh.som_zon(dir),mesh.som_coo(dir,1), mesh.som_coo(dir,2));


% ---- Pseudo elimination
F = F-A*u;


% ---- Resolution du systeme lineaire Au=F
u(inconnues) = A(inconnues, inconnues)\F(inconnues); 


% ---- Affichage de la solution
tri = mesh.elm_som;
x = mesh.som_coo(:,1);
y = mesh.som_coo(:,2);

trimesh(tri, x, y, u);
