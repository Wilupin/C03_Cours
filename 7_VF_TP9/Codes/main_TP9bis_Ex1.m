clear all; 


% INITTIALISATION DES PARAMETRES ET DU PROBLEME

% ---- lecture du maillage
mesh = lect_mesh('../Meshs/carre11'); 

% ---- Matrice de masse
M = assemb_M(mesh); 

% ---- Condition de bords
Gammam = zeros(0,0); 
g = @(x,y) 0; 

% ---- Coordonnees des sommets
x1 = mesh.som_coo(:,1); 
x2 = mesh.som_coo(:,2); 

% ---- Condition initiale et vitesse de propagation
u = [-sin(pi.*x2).*cos(0.5.*pi.*x1),  +sin(pi.*x1).*cos(0.5.*pi.*x2)];
y = max(0, 4*(0.25 - sqrt((x1-0.5).^2 + x2.^2)));

% ---- Variables pour l'affichage
tri = mesh.elm_som;
trisurf(tri, x1, x2, y); 

% ---- Parametres de temps
Tmax = pi; 
dt = 0.01; 
N_pt = floor(Tmax/dt);
t = 0;



% BOUCLE SUR LE TEMPS

for k = 1:N_pt    
    
    % Mise a jour de la solution
    y = convect(mesh, u, y, Gammam, g, dt); 
    
    % ---- Varification de l'evolution de la masse
    Masse = y*M*ones(mesh.nbs,1); 
    
    % ---- Affichage de la solution
    clf(); 
    h_tri = trisurf(tri, x1, x2, y); 
    light
    lighting gouraud
    material dull
    shading interp
    set(h_tri, 'EdgeColor', 'none');
    axis([-1 1 -1 1 0 1]);
    drawnow();
    
    % ---- Affichage de l'evolution de la masse
    fprintf('Iteration     : %i       | Temps        :  %f  | Masse      : %f \n', k, t, Masse);
    
    t = t + dt;
    
end 