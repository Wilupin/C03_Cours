function res = sweby(a,b, beta)

res = sign(a).* max(0, sign(a.*b)).*...
    max(min(abs(a), beta.*abs(b)), min(beta.*abs(a), abs(b)));

end

