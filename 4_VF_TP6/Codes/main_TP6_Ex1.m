clear all; 

N = 1600; 
T = 4.0;

h = 1.0/N;
nu = 0.5;
a = 1; 
dt = nu*h/a; 

LEFT  = [N, 1:N];
RIGHT = [1:N, 1];


% varphi = @(nu) (1);      % Upwind
%varphi = @(nu) (nu);     % LW
varphi = @(nu) (1/nu);   % LFM
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

for k = 1:N_temps
    
    
phi = a*(u(LEFT)+u(RIGHT))/2 - 0.5*abs(a)*varphi(nu)*(u(RIGHT)-u(LEFT));
u = u - (dt/h)*(phi(2:N+1)-phi(1:N));


clf(); plot(xj, u, 'o-'); grid; 

% dim = [.2 .5 .3 .3];
% str = ['t = ' , num2str(t)];
% h_ann = annotation('textbox',dim,'String',str);
% set(h_ann, 'EdgeColor', 'none');
drawnow; 

t = t+dt;
    
end

