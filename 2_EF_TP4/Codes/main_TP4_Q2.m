clear all;

% ---- Lecture du maillage
mesh = lect_mesh('../Meshs/L0');
mesh = raf_mesh(mesh);


% ---- Choix de kappa
kappa = ones(mesh.nbt,1);


% ---- Assemblage des matrices et vecteurs
A = assemb_A(kappa, mesh);
M = assemb_M(mesh);
F = assemb_F(@(x,y) 1, mesh); 


% ---- Identification des bords
dir = find(mesh.som_zon == 1);
inconnues = setdiff(1:mesh.nbs, dir);


% ---- Diagonalisation du probleme
A_tronc = A(inconnues, inconnues); 
M_tronc = M(inconnues, inconnues);

[EigVec, EigVal] = eigs(A_tronc,M_tronc,20,9.0);


% ---- Verification de l'absence de vp degenerees
nb_dv = find(diff(diag(EigVal)) == 0); 

if (isempty(nb_dv))
    disp(['Il n''y a pas de valeurs propres degenerees.']);
else
    disp(['Il y a ', num2str(nb_dv), ' valeurs propres degenerees.']);
end


% ---- Normalisation des vecteurs propres
[maxEigVec, maxInd] = max(abs(EigVec)); 
phi = zeros(mesh.nbs,20);

for k=1:20
    phi(inconnues,k) = (sign(EigVec(maxInd(k),k)).*EigVec(:,k))/abs(maxEigVec(k)); 
end


% ---- Resolution du systeme linaire
alpha = zeros(20,20);

A_bis = phi(inconnues,:)'*A_tronc*phi(inconnues,:);
F_bis =  phi(inconnues,:)'*F(inconnues);

alpha = A_bis\F_bis;

sol = zeros(mesh.nbs,1);

for k = 1:20
     sol(inconnues) = sol(inconnues) + alpha(k).*phi(inconnues,k); 
end


% ---- Affichage de la solution
figure(1);
tri = mesh.elm_som;
x = mesh.som_coo(:,1);
y = mesh.som_coo(:,2);
trimesh(tri,x,y,sol);


% ---- Affichage des 20 premiers vecteurs propres
figure(2);
for j=1:20
    subplot(4,5,21-j)
    trimesh(tri, x, y, phi(:,j));
end
