function g = function_g(z,x,y)

g = zeros(max(size(z)),1);

ind_1 = find(z==1);
ind_2 = find(z==2); 
ind_3 = find(z==3);

% g(ind_1) = x(ind_1);
% g(ind_2) = y(ind_2);
% g(ind_3) = 1-y(ind_3);

g(ind_1) = 0;
g(ind_2) = 2;



end

