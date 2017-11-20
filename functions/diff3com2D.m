% FCS model function 
% -----------------------------------------------------------------
% Copyright MIT 2012
% Developed by Syuan-Ming Guo
% Laboratory for Computational Biology & Biophysics
% Apr 06, 2012
% -----------------------------------------------------------------
function F = diff3com2D(a,t)

    F = a(7) + abs(a(1))./(1 + t./abs(a(4)))...
        +  abs(a(2))./(1 + t./abs(a(5)))...
        +  abs(a(3))./(1 + t./abs(a(6))); 
    if any(~isreal(F))
error('MODELFUN has returned complex values.');
    end
end