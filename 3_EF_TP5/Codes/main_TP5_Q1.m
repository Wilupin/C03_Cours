clear all;

% ---- Lecture du mesh
mesh = lect_mesh('../Meshs/L2');
mesh = raf_mesh(mesh);


% ---- Variables de temps
Delta_t = 0.1;
T = 5; 
tau = 1;
maxTemps = floor(T/Delta_t);


% ---- Choix de kappa
kappa = ones(mesh.nbt,1);


% ---- Assemblage des matrices
A = assemb_A(kappa, mesh);
M = assemb_M(mesh);


% ---- Recuperation des inconnues
dir = find(mesh.som_zon ~= 0);
inconnues = setdiff(1:mesh.nbs, dir);


% ---- Nouvelle matrice A_bis
A_bis = A+(1/Delta_t)*M;


% ---- Creation de la figure
h_fig = figure(1);


% ---- Fonction f
fun = @(x,y,t) sin(2*pi*(x+y-t/tau)); 


% ---- Initialisation du temps et de u 
t = 0; 
u = zeros(mesh.nbs,1); 
 
% ---- Boucle sur le temps
for n = 1:maxTemps
    
 F = assemb_F_bis(fun,mesh,t);
 F = F + (1/Delta_t)*M*u;
 
 u(inconnues) = A_bis(inconnues, inconnues)\F(inconnues);

 
 % ---- Affichage de la solution
 tri = mesh.elm_som;
 x = mesh.som_coo(:,1);
 y = mesh.som_coo(:,2);

 clf(); 
 h = trisurf(tri, x, y, u);
 
 % Option de visualisation de de la solution
 shading interp
 light;
 lighting gouraud
 material dull
 set(h, 'EdgeColor', 'none');
 view(4.5,70);

 drawnow();
 
 
 % Incrementation du temps
 t = t + Delta_t; 
 
end