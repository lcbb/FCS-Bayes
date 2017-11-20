% FCS model function 
% k = wz/wx, structure factor
% -----------------------------------------------------------------
% Copyright MIT 2012
% Developed by Syuan-Ming Guo
% Laboratory for Computational Biology & Biophysics
% Apr 06, 2012
% -----------------------------------------------------------------
function F = diff_flow1com3D(a,t,k)

    F = a(4)+ abs(a(1))./(1 + t./abs(a(2))).*(1+ t./abs(a(2).*k^2)).^(-0.5)...
        .*exp(-(t./abs(a(3))).^2./(1 + t./abs(a(2))));  
if any(~isreal(F))
error('MODELFUN has returned complex values.');
end

end