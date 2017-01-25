function sol = init( u0, mesh )

sol = u0(mesh.elm_gra(:,1), mesh.elm_gra(:,2));

end

