clear all;

% Note : dans la pratique il est plus souvnt utilise assemb_A_Robin_bis,
% qui est un codage plus intuitif de la solution. De m?me, on preferera
% utilise assemb_F_Robine_bis. 


% ---- Lecture du maillage
mesh = lect_mesh('../Meshs/DOM2');


% ---- Choix de kappa et de alpha
kappa = ones(mesh.nbt,1);
alpha = 10^8;


% ---- Assemblage de la matrice de rigidite
%A = assemb_A_Robin(kappa,alpha,mesh);
A = assemb_A_Robin_bis(kappa,alpha,mesh);


% ---- Definition de ua 
ua = @(z,x,y) ((z==1)*(0.0) + (z==2)*(2.0) + (z==3)*(-1.0));


% ---- Assemblage du second membre
%F = assemb_F_Robin(@(x,y) -1, alpha, ua ,@(x,y) 0,mesh);
F = assemb_F_Robin_bis(@(x,y) -1, alpha, ua ,@(x,y) 0,mesh);


% ---- Initialisation de l'inconnue
u = zeros(mesh.nbs,1);


% ---- Resolution du systeme lineaire
% Ici on a besoin de resoudre sur l'ensemble des noeuds donc pas de
% pseudo-elimination avant
u = A\F;


% ---- Affichage de la solution
tri = mesh.elm_som;
x = mesh.som_coo(:,1);
y = mesh.som_coo(:,2);

trimesh(tri, x, y,u);
