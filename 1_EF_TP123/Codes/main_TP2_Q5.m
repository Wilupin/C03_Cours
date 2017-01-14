clear all;

% Lecture du mesh
mesh = lect_mesh('CDR');

mesh = raf_mesh(mesh);
mesh = raf_mesh(mesh);
%mesh = raf_mesh(mesh);

nu = 0.02; 

% choix de kappa
kappa = ones(mesh.nbt,1)*nu;

% Assemblage de la matrice de rigidite
A = assemb_A(kappa, mesh);
spy(A);

% Assemblage du second membre
F = assemb_F(@(x,y) 1, mesh); 

% Assemblage de la matrice de masse 
M = assemb_M(mesh);

% Assemblage de la matrce de convection
C = assemb_C([1,0],mesh);


% Initialisation de l'inconnue
u = zeros(mesh.nbs,1);


% Recuperation des donnees au bord
g = @(z,x,y) ((z==3)*(0.0) +(z==4)*(1.0));
dir = find(mesh.som_zon == 3 | mesh.som_zon == 4);
u(dir) = g(mesh.som_zon(dir),mesh.som_coo(dir,1), mesh.som_coo(dir,2));
inconnues = setdiff(1:mesh.nbs, dir);

% Pseudo elimination
F = F-(A+C+M)*u;

% Resolution du systeme lineaire
u(inconnues) = (A(inconnues, inconnues) + C(inconnues, inconnues) + ...
    M(inconnues, inconnues))\F(inconnues); 


tri = mesh.elm_som;
x = mesh.som_coo(:,1);
y = mesh.som_coo(:,2);

trimesh(tri, x, y, u);


% Calcul de la solution excate 

figure;


l1 = (+1+ sqrt(4*nu+1))/2*nu;
l2 = (-1+ sqrt(4*nu+1))/2*nu;

u = 1 + (exp(l1*x)/(exp(l1+l2)-1)) - (exp(-l2*x)/(1-exp(-l1-l2)));

trimesh(tri, x, y, u); 


