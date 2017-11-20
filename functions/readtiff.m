% This function reads the multi-frame TIFF file
function [im5, im_inf] = readtiff(path_name, range)
im_inf = imfinfo(path_name);
if nargin<2
    range = [1 size(im_inf,1)] ;    
else
end
o.num_frames = range(2)-range(1)+1 ;
im5 = zeros(im_inf(1).Height, im_inf(1).Width, o.num_frames, 'uint16');
h=waitbar(0,'loading data') ;
for j = 1:o.num_frames
im5(:,:,j) = imread(path_name,'index',j-1+range(1),'info',im_inf) ;
    if mod(j,(o.num_frames/100))==0
        waitbar(j/o.num_frames)
    end
end
close(h)
end