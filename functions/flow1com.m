% -----------------------------------------------------------------
% Copyright MIT 2012
% Developed by Syuan-Ming Guo
% Laboratory for Computational Biology & Biophysics
% Apr 06, 2012
% -----------------------------------------------------------------
function F = flow1com(a,t)
F = a(3)+ abs(a(1)).*exp(-(t./abs(a(2))).^2);  
if any(~isreal(F))
error('MODELFUN has returned complex values.');
end

end