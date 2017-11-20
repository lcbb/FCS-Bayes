function F = diff1com3D_free_k(a,t)

% Used by the curve fitter to calculate values for the diffusion equation
% k = wz/wx, structure factor
% k=3.6 ;
    F = a(3)+ abs(a(1))./(1 + t./abs(a(2))).*(1+ t./abs(a(2).*a(4)^2)).^(-0.5);  
if any(~isreal(F))
error('MODELFUN has returned complex values.');
end

% if ~isfinite(F)
%     error('MODELFUN has returned complex values.');
%        a
%        F
end