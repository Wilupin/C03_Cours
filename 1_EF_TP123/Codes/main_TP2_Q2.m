% Script_Ex1.m
%
% Lecture du maillage 'DOM1.amdba'
mesh = lect_mesh('../Meshs/DOM1');
% 1 raffinement ...
mesh = raf_mesh(mesh);
tri = mesh.elm_som;
x   = mesh.som_coo(:,1);
y   = mesh.som_coo(:,2);
% Trace du maillage
triplot(tri, x, y);
%
kappa = ones(mesh.nbt, 1); % kappa(ie)=1 
A = assemb_A_bis(kappa, mesh);
spy(A);
% Declaration de f
f = @(x,y) (0*x+0*y); % fonction anonyme
F = assemb_F(f, mesh);
% Declaration de g
% 
g = @(z,x,y) ( (z==1)*(0.0)+(z==2)*(2.0) );
% Initialisation de u
u = zeros(mesh.nbs, 1); % de taill nbs;
dir = find(mesh.som_zon==1 | mesh.som_zon==2);
% Chargement de u sur Gamma_Dirichlet
u(dir) = g(mesh.som_zon(dir), x(dir), y(dir));
inconnues = setdiff((1:mesh.nbs)', dir);
% Prise en compte des CL de Dirichlet dans 2nd membre
F = F - A * u;
% Resolution du systeme lineaire !
u(inconnues) = A(inconnues, inconnues) \ F(inconnues);
% Visu
figure(1); clf; colormap(cool);
trimesh(tri, x, y, u);