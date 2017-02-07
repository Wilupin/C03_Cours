clear all; 

mesh = lect_mesh('../Meshs/dam0');
mesh = face_number(mesh);


% Initialisation des variables de temps
T  = 200;

dt = 0.1; % Pour dam0
%dt = 0.1; % Pour dam1
%dt = 0.05; % Pour dam2 

N_pt = floor(T/dt); 

V_t = init_sw(mesh); 


% ---- Variables pour la visualisation
x = mesh.som_coo(:,1); 
y = mesh.som_coo(:,2); 

tri = mesh.elm_som;

t = 0;

for k=1:N_pt
    
    fprintf('Iteration     : %i       | Temps        :  %f | Hauteur max : %f \n', k, t, max(V_t(1,:))); 
    fprintf('Vitesse x max : %f  | Vitesse y max : %f | Hauteur min : %f \n', ...
        max(V_t(2,:)./V_t(1,:)), max(V_t(3,:)./V_t(1,:)), min(V_t(1,:)));  
    fprintf('Vitesse x min : %f | Vitesse y min : %f \n\n', ...
        min(V_t(2,:)./V_t(1,:)), min(V_t(3,:)./V_t(1,:)));  
    
    V_t = conv_sw_lagrange(mesh, V_t, dt);
    V_s = tri_to_som(mesh, V_t(1,:));
    
    clf();
    h_tri = trisurf(tri, x, y, V_s);
    light
    lighting gouraud
    material dull
    %shading interp
    set(h_tri, 'EdgeColor', 'none');
    view(52,26);
    axis([0, 200, 0, 200, 0, 10]);
    drawnow;
    
    t = t + dt; 
end