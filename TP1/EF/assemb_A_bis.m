function [A] = assemb_A(K,mesh)
% A=ASSEMB_A(K,mesh) assemble et retourne la matrice de rigidite EF P1 de Lagrange 
% c'est a dire la matrice EF associe a l'operateur -grad (K grad U)
% sur le maillage mesh ou mesh est une structure contenant les champs 
% nbs,nbt,elm_som,som_coo,som_zon
%
% Ne tient pas compte des Conditions aux Limites
% K est suppose constant par element et transmis sous forme de tableau
% colonne à mesh.nbt lignes

% Copyright (c) 2005 Frederic Pascal, 2015 Florian De Vuyst, ENS Cachan

% Allocation memoire d'une matrice creuse
A  = sparse(mesh.nbs,mesh.nbs);

% Assemblage de la matrice
for ie = 1:mesh.nbt
  is = mesh.elm_som(ie,:);
  x  = mesh.som_coo(is,:);
  
  a  = - x([ 3 1 2],2) + x([2 3 1],2);
  b  =   x([ 3 1 2],1) - x([2 3 1],1);
  %
  mes = 0.5 * (b(2)*a(1) - a(2)*b(1));
  a = a / (2*mes);
  b = b / (2*mes);
  % Matrice elementaire 3x3
  %
  A_elmK = mes * K(ie) * (a*a' + b*b');
  
  % Contribution matrice globale
  A(is,is) = A(is,is) + A_elmK;
end

