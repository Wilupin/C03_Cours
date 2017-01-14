clear all;

% ---- Lecture du maillage
mesh = lect_mesh('../Meshs/L0');
mesh = raf_mesh(mesh);


% ---- Choix de kappa
kappa = ones(mesh.nbt,1);

% ---- Assemblage des matrices
A = assemb_A(kappa, mesh);
M = assemb_M(mesh);
%spy(A);


% ---- Identification du bord
dir = find(mesh.som_zon == 1);
inconnues = setdiff(1:mesh.nbs, dir);


% ---- Diagonalisation du probleme 
A_tronc = A(inconnues, inconnues); 
M_tronc = M(inconnues, inconnues);

[EigVec, EigVal] = eigs(A_tronc,M_tronc,6,9.0);

[maxEigVec, indMax] = max(EigVec); 
phi = zeros(mesh.nbs,6);


% ---- Normalisation des vecteurs propres
for k=1:6
    phi(inconnues,k) = (sign(EigVec(indMax(k)))*EigVec(:,k))/maxEigVec(k); 
end


% ---- Affichage des 6 premi?res valeurs propres
tri = mesh.elm_som;
x = mesh.som_coo(:,1);
y = mesh.som_coo(:,2);

for j=1:6
    subplot(2,3,7-j)
    trimesh(tri, x, y, phi(:,j));
end
