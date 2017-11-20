function [imA imB] = read_PAT(filename, opt)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Read photon arrival time and covert it to 8-bits photon count trace%
% ----Required input------------------------------------------------------
% filename: path and filename of the PAT file without "A.dat" (or "B.dat")
% opt: option, a struture with following fields: 
% .bin: integer indicating the binning size for the photon count trace.
% Increase this value when the memory is not enough to process the raw
% trace.
% .bits: integer indicating the data type (8 or 16 bits) ;
% .hdr: integer indicating number of values to remove from the beginning
% (header).
% .n_trc: number of traces. Use 1 if the data contain a single PAT and 2 if
% the data contain 2 PATs.
% ----Output-------------------------------------------------------
% imA,imB: photon count traces. imB is only output when opt.n_trc = 2.
%-----------------------------------------------------------------
% Ref: S. Guo et. al., Anal Chem, 2012
% -----------------------------------------------------------------
% Syuan-Ming Guo
% sguo@mit.edu, Sep 13, 2012
% -----------------------------------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if opt.n_trc ==2 
    fid = fopen([filename 'A.dat']);
    patA = fread(fid, ['*uint' num2str(opt.bits)]);       % Output has the same class as input. 
    fclose(fid);
    
    fid = fopen([filename 'B.dat']);
    patB = fread(fid, ['*uint' num2str(opt.bits)]);   % Output has the same class as input. 
    fclose(fid);
    
    if opt.hdr > 0
        patA(1:opt.hdr) = [] ;   % remove the header
        patB(1:opt.hdr) = [] ;   % remove the header
    end
             
    patA_c = cumsum(double(patA));   % convert photon arrival intervals to photon arrival time 
    patB_c = cumsum(double(patB));

    patA_c = ceil(patA_c/opt.bin) ;  % bin the arrival time
    patB_c = ceil(patB_c/opt.bin) ;

    im_length = max(patA_c(end), patB_c(end)) ;
    
    imA = zeros(1,im_length,'uint8');
    imB = zeros(1,im_length,'uint8');
    for j = 1:numel(patA_c)
        imA(patA_c(j)) = imA(patA_c(j))+1 ;
    end

    for j = 1:numel(patB_c)
        imB(patB_c(j)) = imB(patB_c(j))+1 ;
    end

elseif opt.n_trc ==1 
    fid = fopen([filename '.dat']);
    patA = fread(fid, ['*uint' num2str(opt.bits)]);       % Output has the same class as input. 
    fclose(fid);
    if opt.hdr > 0
        patA(1:opt.hdr) = [] ;   % remove the header        
    end
    
    patA_c = cumsum(double(patA));

    patA_c = ceil(patA_c/opt.bin) ;
    im_length = patA_c(end);

    imA = zeros(1,im_length,'uint8');    
    
    for j = 1:numel(patA_c)
        imA(patA_c(j)) = imA(patA_c(j))+1 ;
    end
else
    error('wrong option')
end

end

