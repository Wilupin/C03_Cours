clear all;

% Lecture du mesh
mesh = lect_mesh('L0');
mesh = raf_mesh(mesh);


% choix de kappa
kappa = ones(mesh.nbt,1);

% Assemblage de la matrice de rigidite
A = assemb_A(kappa, mesh);
spy(A);

% Assemblage de la matrice de masse 
M = assemb_M(mesh);

dir = find(mesh.som_zon == 1);
inconnues = setdiff(1:mesh.nbs, dir);

% Resolution du systeme lineaire
A_tronc = A(inconnues, inconnues); 
M_tronc = M(inconnues, inconnues);

[EigVec, EigVal] = eigs(A_tronc,M_tronc,6,9.0);

[maxEigVec, indMax] = max(EigVec); 
phi = zeros(mesh.nbs,6);


for k=1:6
    phi(inconnues,k) = (sign(EigVec(indMax(k)))*EigVec(:,k))/maxEigVec(k); 
end

tri = mesh.elm_som;
x = mesh.som_coo(:,1);
y = mesh.som_coo(:,2);

for j=1:6
    subplot(2,3,7-j)
    trimesh(tri, x, y, phi(:,j));
end
