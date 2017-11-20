% This funtion performs fast binning on the input 1D vector with "bin_size"
% "opt" specifies the variable type
% "opt" can be either 'native' or 'double'
% -----------------------------------------------------------------
% Copyright MIT 2012
% Developed by Syuan-Ming Guo
% Laboratory for Computational Biology & Biophysics
% Apr 06, 2012
% -----------------------------------------------------------------
function x_binned = binning_1d(x,bin_size,varargin)

if ~isempty(varargin)
    opt = varargin{1} ;    
else
    opt = 'double' ;
end

res = mod(length(x),bin_size) ;
% num_new  = floor(length(x)/bin_size);
x = x(1:length(x)-res) ;
x = reshape(x,bin_size,[]);
x_binned = sum(x,1,opt) ;
end