function F = diff2com3D_free_k(a,t)

% Used by the curve fitter to calculate values for the diffusion equation
% k=3.6 ;
    F = a(5) + a(1)./(1 + t./abs(a(3))).*(1+ t./(abs(a(3)).*a(6)^2)).^(-0.5)...
        + a(2)./(1 + t./abs(a(4))).*(1+ t./(abs(a(4)).*a(6)^2)).^(-0.5); 
    if any(~isreal(F))
error('MODELFUN has returned complex values.');
    end