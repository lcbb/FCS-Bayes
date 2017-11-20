% FCS model function
% -----------------------------------------------------------------
% Copyright MIT 2012
% Developed by Syuan-Ming Guo
% Laboratory for Computational Biology & Biophysics
% Apr 06, 2012
% -----------------------------------------------------------------
function F = diff3com2D_im(a,t,L,sigma)
m0 = mf(L,0,0,sigma) ;
m1 = mf(L,abs(a(4)),t,sigma) ;
m2 = mf(L,abs(a(5)),t,sigma) ;
m3 = mf(L,abs(a(6)),t,sigma) ;
F = a(7)+ abs(a(1)).*((erf(m1) + 1/sqrt(pi)./m1.*(exp(-m1.^2)-1))/...
    (erf(m0) + 1/sqrt(pi)./m0.*(exp(-m0.^2)-1))).^2 ...
    + abs(a(2)).*((erf(m2) + 1/sqrt(pi)./m2.*(exp(-m2.^2)-1))/...
    (erf(m0) + 1/sqrt(pi)./m0.*(exp(-m0.^2)-1))).^2 ...
    + abs(a(3)).*((erf(m3) + 1/sqrt(pi)./m3.*(exp(-m3.^2)-1))/...
    (erf(m0) + 1/sqrt(pi)./m0.*(exp(-m0.^2)-1))).^2;

end

function m = mf(L,D,t,sigma)

m = L/2./sqrt(D.*t + sigma^2);  


end