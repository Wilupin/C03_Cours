clear all;

% Lecture du mesh
mesh = lect_mesh('L2');
mesh = raf_mesh(mesh);

% Variables de temps
Delta_t = 0.1;
T = 5; 
tau = 1;

maxTemps = floor(T/Delta_t);


% choix de kappa
kappa = ones(mesh.nbt,1);

% Assemblage de la matrice de rigidite
A = assemb_A(kappa, mesh);

% Assemblage de la matrice de masse
M = assemb_M(mesh);

% Recuperation des inconnues
dir = find(mesh.som_zon ~= 0);
inconnues = setdiff(1:mesh.nbs, dir);

% Nouvelle matrice A_bis
A_bis = A+(1/Delta_t)*M;


% Creation de la figure
h_fig = figure(1);

% Fonction f
fun = @(x,y,t) sin(2*pi*(x+y-t/tau)); 

% Initialisation du temps et de u 
t = 0; 

u = zeros(mesh.nbs,1); 
 

for n = 1:maxTemps
    
 F = assemb_F_bis(fun,mesh,t);
 F = F + (1/Delta_t)*M*u;
 
 u(inconnues) = A_bis(inconnues, inconnues)\F(inconnues);

 tri = mesh.elm_som;
 x = mesh.som_coo(:,1);
 y = mesh.som_coo(:,2);

 trimesh(tri, x, y, u);
 view(4.5,70);
 
 pause(0.01);
 
 t = t + Delta_t; 
 
end