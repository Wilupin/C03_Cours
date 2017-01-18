function res = minmod(a,b)

res = sign(a).*max(0, sign(a.*b)).*min(abs(a), abs(b)); 

end

