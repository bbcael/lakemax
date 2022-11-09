function g = gfbm(q);
g = q.^4./6.*hypergeom([1 1],[5/2 3],q.^2./2)-3.*q.^2 + pi.*(1-q.^2).*erfi(q./sqrt(2)) + sqrt(2*pi)*exp(q.^2/2).*q + (q.^2-2).*(log(2.*q.^2) + 0.577215664901532);
end