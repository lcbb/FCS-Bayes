function F = diff3com3D_free_k(a,t)

% Used by the curve fitter to calculate values for the diffusion equation
% k=3.6 ;
    F = a(7) + a(1)./(1 + t./abs(a(4))).*(1+ t./(abs(a(4)).*a(8)^2)).^(-0.5)...
        + a(2)./(1 + t./abs(a(5))).*(1+ t./(abs(a(5)).*a(8)^2)).^(-0.5)...
        + a(3)./(1 + t./abs(a(6))).*(1+ t./(abs(a(6)).*a(8)^2)).^(-0.5); 
    if any(~isreal(F))
error('MODELFUN has returned complex values.');
    end
end