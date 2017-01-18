clear all; 

N = 1600; 
T = 4.0;

h = 1.0/N;

% Quand on fait une reconstitution de pente l'analyse nous montre : 
% nu < 0.5. 'On travaille sur des demi cellules'
nu = 0.5;
a = 1; 
dt = nu*h/a; 

LEFT  = [N, 1:N];
RIGHT = [1:N, 1];


varphi = @(nu) (1);      % Upwind
% varphi = @(nu) (nu);     % LW
% varphi = @(nu) (1/nu);   % LFM
% varphi = @(nu) sqrt(nu); 


% Position des centres
xj = h/2 : h : 1-h/2;

% Donnée initiale
u0 = max(0,sin(6*pi*xj)).*(xj<=1/3) + ...
    (3*xj-1).*(xj>1/3).*(xj<2/3) + ...
    1*(xj<=1).*(xj>2/3);


% Initialisation de u
u = zeros(N,1); 
u = u0; 

t = 0; 

N_temps = floor(T/dt);


% ---- Boucle sur le temps 
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


% ----- Plot de la solution
clf(); plot(xj, u, '-'); grid; 
drawnow; 

t = t+dt;
    
end
