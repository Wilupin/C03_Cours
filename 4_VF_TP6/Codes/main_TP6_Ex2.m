% -------------------------------------
% Gaetan Facchinetti
% Cours C03 : methodes numeriques
%
% TP6 : Volumes finis en 1D
% -> Avec limiteur de pente
% -------------------------------------


% ---- INITIALISATION DES PARAMETRES INITIAUX

% Nombre de points, temps max
N = 1600; 
T = 4.0;

% Pas d'espace : longeur d'un intervalle
h = 1.0/N;

% Positions des centres des intervalles
xj = h/2 : h : 1-h/2;

% Parametre de flux
nu = 0.45;

% Vitesse de l'equation des ondes
a = 1; 

% Pas de temps
dt = nu*h/a; 

% Tableaux de translation des points
LEFT  = [N, 1:N];
RIGHT = [1:N, 1];

% Donnee initiale de la solution
u = max(0,sin(6*pi*xj)).*(xj<=1/3) + ...
    (3*xj-1).*(xj>1/3).*(xj<2/3) + ...
    1*(xj<=1).*(xj>2/3);


% Temps initial et pas de temps
t = 0; 
N_temps = floor(T/dt);




% ---- DEFINITION DU FLUX NUMERIQUE

varphi = @(nu) (1);      % Upwind
% varphi = @(nu) (nu);     % LW
% varphi = @(nu) (1/nu);   % LFM
% varphi = @(nu) sqrt(nu); 



% ---- BOUCLE SUR LE TEMPS

for k = 1:N_temps
    
sj = minmod(u([2:N,1])-u, u-u([N,1:N-1]));  
% sj = sweby(u([2:N,1])-u, u-u([N,1:N-1]),2); 

uplus  = u + 0.5*sj; 
umoins = u - 0.5*sj;

% ----- Reconstitution du flux numerique
phi = a*(uplus(LEFT)+umoins(RIGHT))/2 -...
    0.5*abs(a)*varphi(nu)*(umoins(RIGHT)-uplus(LEFT));

% ----- Iteration de la solution
u = u - (dt/h)*(phi(2:N+1)-phi(1:N));


% ----- Affichage de la solution
clf(); plot(xj, u, '-'); grid; 
drawnow; 

t = t+dt;
    
end
