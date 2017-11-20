% -----------------------------------------------------------------
% Copyright MIT 2012
% Developed by Syuan-Ming Guo
% Laboratory for Computational Biology & Biophysics
% Apr 06, 2012
% -----------------------------------------------------------------
function F = diff3com_1_2D_2_3D(a,t,k)

% Used by the curve fitter to calculate values for the diffusion equation
% k=3.6 ;
    F = a(7) + a(1)./(1 + t./abs(a(4)))...
        + a(2)./(1 + t./abs(a(5))).*(1+ t./(abs(a(5)).*k^2)).^(-0.5)...
        + a(3)./(1 + t./abs(a(6))).*(1+ t./(abs(a(6)).*k^2)).^(-0.5); 
    if any(~isreal(F))
error('MODELFUN has returned complex values.');
    end
end