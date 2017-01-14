function [F] = assemb_F_Robin(fun,alpha,ua,fung,mesh)
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
  
  g = sum(x,1)/3; % Centre de gravite
  
  F_elm = mesK * fun(g(1),g(2)) / 3.0 * ones(3,1);
  
  F(is) = F(is) + F_elm;
  
  
  id = find(mesh.som_zon(is) ~= 0);
  
  if (max(size(id)) == 2)
      z = mesh.som_zon(is(id(1)));
      d = sqrt((x(id(1),1)-x(id(2),1)).^2 + (x(id(1),2)-x(id(2),2)).^2);
      g2 = sum(x(id,:))/2;
      Flm_bord = d*(fung(g2(1),g2(2))+alpha*ua(z,g2(1),g2(1)))*0.5*ones(2,1);
      F(is(id)) = F(is(id)) + Flm_bord;
  end
  
  
end  
