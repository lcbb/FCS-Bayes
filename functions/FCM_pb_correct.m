% This function corrects photobleaching for individual intensity traces
% by fitting a double exponential function to the trace  corrected
% individualy. correction is done by subtracting the fit from the trace
% opt == 1 is nlinfit, 2 is lsqcurvefit
% -----------------------------------------------------------------
% Copyright MIT 2013
% Developed by Syuan-Ming Guo
% Laboratory for Computational Biology & Biophysics
% Oct 31, 2013
% -----------------------------------------------------------------
function imcc = FCM_pb_correct(imc,opt)
frame_v = [1:size(imc,3)]' ;  % frame number vector
imcc = zeros(size(imc));
h = waitbar(0,'Correcting photobleaching...');

switch opt
    case 'pixelwise1' ;
    for j = 1:size(imc,1)
        for k = 1:size(imc,2)
            Ib1 = squeeze(imc(j,k,:)) ;
            beta0 = [0.1 0.1 0.5*max(Ib1) 0.5*max(Ib1)] ; 
            beta = nlinfit((frame_v-1),Ib1,@bi_exponential,beta0) ;
            Ib_fit = bi_exponential(beta,(frame_v-1)) ;        
            imcc(j,k,:) = Ib1-Ib_fit + mean(Ib_fit); 
            waitbar((size(imc,2)*(j-1)+k)/numel(imc(:,:,1)),h)
        end
    end

    case 'pixelwise2' ;
    for j = 1:size(imc,1)
        for k = 1:size(imc,2)
            Ib1 = squeeze(imc(j,k,:)) ;
            beta0 = [0.1 0.0001 0.1*max(Ib1) 0.1*max(Ib1) mean(Ib1(end-10:end))] ; 
            beta = lsqcurvefit(@bi_exponential,beta0,(frame_v-1),Ib1,[0 0 0 0 min(Ib1)],[100 100 max(Ib1) max(Ib1) max(Ib1)]) ;
            Ib_fit = bi_exponential(beta,(frame_v-1)) ;    
            imcc(j,k,:) = Ib1-Ib_fit + mean(Ib_fit);
            waitbar((size(imc,2)*(j-1)+k)/numel(imc(:,:,1)),h)
        end
    end        
end
close(h)
end
