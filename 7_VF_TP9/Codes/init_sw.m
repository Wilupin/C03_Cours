function Vnew_t = init_sw(mesh)

% Inialise l'inconnue V regroupant vitesse et hauteur
% Vnew_t(1,:) : hauteur du fluide
% Vnew_t(2,:) : hauteur x vitesse en x
% Vnew_t(3,:) : hauteur x vitesse en y

% Coordonnees des sommets
x = mesh.som_coo(:,1);
y = mesh.som_coo(:,1); 


% Initialisation de V 
Vnew_t = zeros(3,mesh.nbt); 

x1 = x(mesh.elm_som(:,1));
x2 = x(mesh.elm_som(:,2));
x3 = x(mesh.elm_som(:,3));

% Recherche de l'amont et de l'aval du barrage
up_bar     = find(max(max(x1,x2),x3)<=100);
down_bar   = setdiff(1:mesh.nbt, up_bar);

% Remplissage du tableau
Vnew_t(1, up_bar)    = 10.0;
Vnew_t(1, down_bar)  = 5.0;

end

