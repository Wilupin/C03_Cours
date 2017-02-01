function [ pstar, unustar, UKL] = varstar( nu, V )

g = 9.81; 

UK = V(:,1); 
UL = V(:,2); 

% ---- Hauteur d'eau
h = zeros(1,2);
h = V(1, :); 

% ---- Celerite
cK = sqrt(g*h(1)); 
cL = sqrt(g*h(2));


% ---- Vitesse du fluide
u = zeros(2,2); 
u(1,:) = V(2, :)./h;
u(2,:) = V(3, :)./h;


% ---- Pression 
pres = zeros(1,2); 
pres = 0.5*g*(h.^2); 

pK = pres(1); 
pL = pres(2); 
hK = h(1); 
hL = h(2);
uK = u(:,1); 
uL = u(:,2);

pstar = (pK*hL + pL* hK)/(hK+hL) - (hK*hL)/(hK+hL).*max(cK, cL).*(uK'*nu - uL'*nu);

unustar = 0.5*(uK'*nu + uL'*nu) - (1/max(cK, cL)).*((pres(2)-pres(1))/(hK+hL));

UKL = UK*max(0, unustar) + UL*min(0, unustar);

end

