% FCS model function 
% k = wz/wx, structure factor
% -----------------------------------------------------------------
% Copyright MIT 2012
% Developed by Syuan-Ming Guo
% Laboratory for Computational Biology & Biophysics
% Apr 06, 2012
% -----------------------------------------------------------------
function F = trip_diff1com3D_free_k(a,t)

    F = a(3)+ a(1).* (1+ abs(a(5))./(1-abs(a(5)))*exp(-t./abs(a(6)))).*...
        (abs(a(2))./abs((a(2)) + t)).*(1+ t./abs((a(2)).*a(4)^2)).^(-0.5);  
end
