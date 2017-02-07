function [F] = assemb_F_Robin_bis(fun,alpha,ua,fung,mesh)
% F=ASSEMB_F(fun,mesh) assemble le second membre du pb du Laplacien
% avec conditions aux limites de Neumann Homogenes
% c'est-a-dire F_i = \int_{\Omega} fun(x,y)*\phi_i dxdy 
% pour la methode des EF P1 de Lagrange
%
% La fonction fun(x,y) est une fonction de 2 arguments 
%

% Allocation
F = zeros(mesh.nbs,1);

% Assemblage du second membre avec une formule a 1 point
for ie = 1:mesh.nbt
  is = mesh.elm_som(ie,:);            % vecteur 1x3
  x  = mesh.som_coo(is,:);            % vecteur 3x2
  a  = -x([ 3 1 2],2)+x([2 3 1],2);   % vecteur 3x1
  b  =  x([ 3 1 2],1)-x([2 3 1],1);
  mesK = 0.5*(b(2)*a(1)-a(2)*b(1));
  
  g = sum(x,1)/3;     % Centre de gravite
  
  F_elm = mesK * fun(g(1),g(2)) / 3.0 * ones(3,1);
  
  F(is) = F(is) + F_elm;
  
end  


for ie = 1:mesh.nbab
    
    is = mesh.abd_som(ie, :);
    
    % Extremites du segment formant l'arrete
    x  = mesh.som_coo(is,1);
    y  = mesh.som_coo(is,2);
    
    % Longueur du segment
    d = sqrt((x(1)-x(2)).^2 + (y(1)-y(2)).^2);
    
    % Coordonnes du centre de gravite du segment
    g2_x = sum(x(:))/2;
    g2_y = sum(y(:))/2;
    
    % Zone a laquelle appartient le segment
    z = mesh.som_zon(is(1)); 
    
    F(is) = F(is) + d*(fung(g2_x,g2_y)+alpha*ua(z,g2_x,g2_y))*0.5*ones(2,1);
    % On multiplie par 0.5*ones(2,1) car les fonctions de base valent 1/2
    % au centre de gravite de l'arr?te
    
end
  