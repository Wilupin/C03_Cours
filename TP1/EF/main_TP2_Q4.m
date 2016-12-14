clear all;

% Lecture du mesh
mesh = lect_mesh('car0');


for k = 1:5

    id = find(mesh.som_zon ~=0 & mesh.som_coo(:,1) ~=0); 
    h(k) = min(mesh.som_coo(id,1));   
    
    % choix de kappa
    kappa = ones(mesh.nbt,1);
    
    % Assemblage de la matrice de rigidité
    A = assemb_A(kappa, mesh);
    spy(A);
    
    % Assemblage du second membre
    F = assemb_F(@(x,y) -4, mesh);
    M = assemb_M(mesh);
    
    % Initialisation de l'inconnue
    u = zeros(mesh.nbs,1);
    
    % Recuperation des donnees au bord
    dir = find(mesh.som_zon == 2);
    inconnues = setdiff(1:mesh.nbs, dir);
    u(dir) = function_g(mesh.som_zon(dir),mesh.som_coo(dir,1), mesh.som_coo(dir,2));
    
    % Autre façon de dire
    % g = @(z,x,y) ((z==2)*1.0*(x.^2 + y.^2));
    % u(dir) = g(mesh.som_zon(dir),mesh.som_coo(dir,1), mesh.som_coo(dir,2));
    % inconnues = setdiff(1:mesh.nbs, dir);
    
    % Pseudo elimination
    F = F-A*u;
    
    % Resolution du système lineaire
    u(inconnues) = A(inconnues, inconnues)\F(inconnues);
    
    
    tri = mesh.elm_som;
    x = mesh.som_coo(:,1);
    y = mesh.som_coo(:,2);
    
    %trimesh(tri, x, y, u);
    
    
    u_vrai = x.^2+y.^2;
    diff = u-u_vrai;
    
    err(k) = sqrt(diff'*M*diff);
    
    if (k~=4)
        mesh = raf_mesh(mesh);
    end

end

figure;

loglog(h, err, 'x');
hold on 
t=min(h):((max(h)-min(h))/10):max(h);
loglog(t,t*(err(1)/h(1)), '-', 'Color', 'red')
hold on 
loglog(t, t.^2*(err(1)/h(1)^2), '-', 'Color', 'blue');
hold on 
loglog(t, t.^3*(err(1)/h(1)^3), '-', 'Color', 'green');
hold on 
loglog(t, t.^4*(err(1)/h(1)^4), '-', 'Color', 'cyan');



% Verification de l'assemblage de M 
% Res doit valoir l'aire du domaine
% v = ones(mesh.nbs,1);
% Res = v'*M*v;

