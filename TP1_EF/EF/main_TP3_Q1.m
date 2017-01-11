clear all;

% Lecture du mesh
mesh = lect_mesh('DOM2');

% choix de kappa
kappa = ones(mesh.nbt,1);
alpha = 10^8;

% Assemblage de la matrice de rigidité
A = assemb_A_Robin(kappa,alpha,mesh);
spy(A);

ua = @(z,x,y) ((z==1)*(0.0) + (z==2)*(2.0) + (z==3)*(-1.0));
% u(dir) = g(mesh.som_zon(dir),mesh.som_coo(dir,1), mesh.som_coo(dir,2))
F = assemb_F_Robin(@(x,y) -1, alpha, ua ,@(x,y) 0,mesh) 

u = zeros(mesh.nbs,1);

u = A\F;

tri = mesh.elm_som;
x = mesh.som_coo(:,1);
y = mesh.som_coo(:,2);

trimesh(tri, x, y,u);


v = ones(mesh.nbs,1);
Res = v'*A*v;
