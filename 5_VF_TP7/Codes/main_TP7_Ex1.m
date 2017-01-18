clear all;


% ---- Initialisation des constantes
delta_t = 0.015; 
T = 2*pi; 
N_T = floor(T/delta_t);

% ---- Lecture du maillage et rajout volumes finis
mesh = lect_mesh('../Meshs/disq0');
mesh = raf_mesh(mesh); 
mesh = face_number(mesh); 

% Attention il y a une condition de stabilite a respecter. Si on divise le
% diametre du maillage par 2 avec un raf_mesh il faut aussi diviser le pas
% de temps par 2 pour ne pas sortir de la condition de stabilite.


% ---- Parametres de visualisation
tri = mesh.elm_som;
x = mesh.som_coo(:,1);
y = mesh.som_coo(:,2);


% ---- Definition de la condition initiale
%u0 = @(x,y) exp(-50*(((x-0.4).^2)+((y-0.0).^2)));
u0 = @(x,y) 1.*((x-0.4).^2+ (y-0.0).^2 <= 0.04);


% ---- Initalisation de la donnee initiale
sold_t = init(u0, mesh);


figure(1)

% ---- Boucle en temps
for it=1:N_T
    
   
    snew_t = conv_sca(mesh, sold_t, delta_t); 
    sold_t = snew_t; 

    
    % ---- Visualisation de la solution
    sol_s = tri_to_som(mesh, sold_t);
  
    clf();
    
    subplot(1,2,1),
    h_tri = trisurf(tri, x, y, sol_s);
    
    % Parametres d'affichage
    light
    lighting gouraud
    material dull
    shading interp
    set(h_tri, 'EdgeColor', 'none');
    axis([-1, 1, -1, 1, 0, 1]);
    
    subplot(1,2,2),
    tri_contour(tri, x, y, sol_s, 0:0.05:1);
    axis([-1, 1, -1, 1]);
    
    drawnow;
end


