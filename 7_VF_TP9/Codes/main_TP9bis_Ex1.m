clear all; 

% ---- lecture du maillage
mesh = lect_mesh('../Meshs/carre11'); 


% ---- Matrice de masse
M = assemb_M(mesh); 


Gammam = zeros(0,0); 
g = @(x,y) 0; 

x1 = mesh.som_coo(:,1); 
x2 = mesh.som_coo(:,2); 

u = [-sin(pi.*x2).*cos(0.5.*pi.*x1),  +sin(pi.*x1).*cos(0.5.*pi.*x2)];
y = max(0, 4*(0.25 - sqrt((x1-0.5).^2 + x2.^2)));

size(y)


tri = mesh.elm_som;
trisurf(tri, x1, x2, y); 


Tmax = pi; 
dt = 0.01; 
N_pt = floor(Tmax/dt);

t = 0;


for k = 1:N_pt    
    
    y = convect(mesh, u, y, Gammam, g, dt); 
    
    Masse = y*M*ones(mesh.nbs,1); 
    
    clf(); 
    h_tri = trisurf(tri, x1, x2, y); 
    light
    lighting gouraud
    material dull
    shading interp
    set(h_tri, 'EdgeColor', 'none');
    axis([-1 1 -1 1 0 1]);
    drawnow();
    
    fprintf('Iteration     : %i       | Temps        :  %f  | Masse      : %f \n', k, t, Masse);
    
    t = t + dt;
    
end 