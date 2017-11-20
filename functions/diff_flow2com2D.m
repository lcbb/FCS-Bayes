% FCS model function 
% -----------------------------------------------------------------
% Copyright MIT 2012
% Developed by Syuan-Ming Guo
% Laboratory for Computational Biology & Biophysics
% Apr 06, 2012
% -----------------------------------------------------------------
function F = diff_flow2com2D(a,t)

    F = a(6)+ abs(a(1))./(1 + t./abs(a(3)))...
        .*exp(-(t./abs(a(5))).^2./(1 + t./abs(a(3))))+...
        abs(a(2))./(1 + t./abs(a(4)))...
        .*exp(-(t./abs(a(5))).^2./(1 + t./abs(a(4))));  
if any(~isreal(F))
error('MODELFUN has returned complex values.');
end

end