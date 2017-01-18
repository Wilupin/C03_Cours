function snew_t = conv_sca( mesh, sold_t, deltat )


% ---- Vitesse de transport
c = @(x,y) ([-y, x]);

% ---- Initialisation de la difference
dsol_t = zeros(mesh.nbt,1); 

% ---- Position du centre des arretes
x_a = mesh.fac_gra(:,1);
y_a = mesh.fac_gra(:,2);

% ---- Normales aux arretes
N1_KL = mesh.fac_mes.*mesh.fac_nor(:,1);
N2_KL = mesh.fac_mes.*mesh.fac_nor(:,2);


% Boucle sur l'ensemble des arretes
for ia = 1:mesh.nba
    
    % Triangles voisins de l'arrete ia
    ie = mesh.fac_elm(ia,:);
    
    if(ie(2) ~= 0)
        % flux numerique avec une arrete interieure
        CN_KL = c(x_a(ia), y_a(ia))*[N1_KL(ia) ; N2_KL(ia)];
        phi =  max(0, CN_KL)*sold_t(ie(1)) + min(0, CN_KL)*sold_t(ie(2));
    else
        % flux numerique avec une arrete du bord
        phi = 0;
    end 
    
    % Contribution au triangle K
    dsol_t(ie(1)) =  dsol_t(ie(1)) + phi;
    
    % Contribution au triangle L (s'il existe) 
    if(ie(2) ~=0 )
        dsol_t(ie(2)) = dsol_t(ie(2)) - phi;
    end

end

snew_t = sold_t - deltat*dsol_t./mesh.elm_mes; 
