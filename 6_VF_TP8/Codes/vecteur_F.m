function F = vecteur_F(nu, V)

g = 9.81;

F = zeros(3,1);

u = zeros(2,1); 
u = V(2:3)./V(1); 
h = V(1);

F = [ h*(u'*nu); ...
    h*u(1)*(u'*nu) + 0.5*g*(h^2)*nu(1); ...
    h*u(2)*(u'*nu) + 0.5*g*(h^2)*nu(2)];

end

