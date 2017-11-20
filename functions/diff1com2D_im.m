% FCS model function
% -----------------------------------------------------------------
% Copyright MIT 2012% 
% Laboratory for Computational Biology & Biophysics
% Apr 06, 2012
% -----------------------------------------------------------------
function F = diff1com2D_im(a,t,L,sigma)
m = mf(L,abs(a(2)),t,sigma) ;
m0 = mf(L,0,0,sigma) ;
F = a(3)+ abs(a(1)).*((erf(m) + 1/sqrt(pi)./m.*(exp(-m.^2)-1))/...
    (erf(m0) + 1/sqrt(pi)./m0.*(exp(-m0.^2)-1))).^2;  

end

function m = mf(L,D,t,sigma)

m = L/2./sqrt(D.*t + sigma^2);  


end