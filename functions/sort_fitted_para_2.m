function [para_s std_s] = sort_fitted_para_2(para, std, model, Aeff, wf)


switch model
    case 'null'        
        para_s.a = abs(para(1,:)) ;
        para_s.tauD = NaN ;
        para_s.D = NaN ;
        para_s.Ginf = NaN ;
        
        std_s.a = abs(std(1,:)) ;
        std_s.tauD = NaN ;
        std_s.D = NaN ;
        std_s.Ginf = NaN ;
                
    case 'flow'
        para(2,:) = abs(para(2,:))*1e6 ;
        para_s.a = abs(para(1,:)) ;
        para_s.tauD = para(2,:) ;
        para_s.D = wf./para_s.tauD*10^6 ;
        para_s.Ginf = para(3,:) ;
        
        std(2,:) = abs(std(2,:))*1e6 ;
        std_s.a = abs(std(1,:)) ;
        std_s.tauD = std(2,:) ;
        std_s.Ginf = std(3,:) ;
        
    case {'1comp3Ddiff' '1comp2Ddiff'}           
        para(2,:) = abs(para(2,:))*1e6 ;
        para_s.a = abs(para(1,:)) ;
        para_s.tauD = para(2,:) ;
        para_s.D = wf^2/4./para_s.tauD*10^6 ;
        para_s.Ginf = para(3,:) ;
        
        std(2,:) = abs(std(2,:))*1e6 ;
        std_s.a = abs(std(1,:)) ;
        std_s.tauD = std(2,:) ;
        std_s.Ginf = std(3,:) ;
        
    case '1comp2Ddiff_im'
        para(2,:) = abs(para(2,:)) ;
        para_s.a = abs(para(1,:)) ;
        para_s.D = para(2,:) ;
%         para_s.tauD = wf^2/4./para_s.D*10^6 ;
        para_s.tauD = Aeff/4./para_s.D*10^6 ;
        para_s.Ginf = para(3,:) ;
        
        std(2,:) = abs(std(2,:));
        std_s.a = abs(std(1,:)) ;
        std_s.D = std(2,:) ;
        std_s.Ginf = std(3,:) ;
        
    case {'1comp2Ddiff_0Ginf' '1comp2Ddiff_0Ginf_im'}           
        para(2,:) = abs(para(2,:))*1e6 ;
        para_s.a = abs(para(1,:)) ;
        para_s.tauD = para(2,:) ;
        para_s.D = wf^2/4./para_s.tauD*10^6 ;
        
        std(2,:) = abs(std(2,:))*1e6 ;
        std_s.a = abs(std(1,:)) ;
        std_s.tauD = std(2,:) ;
        
    case {'2comp3Ddiff' '2comp2Ddiff'}                
        tauD = abs(para(3:4,:)) ;
        alpha = abs(para(1:2,:)) ;       
       [tauD_s,IX] = sort(tauD,1,'ascend') ;
       

        for j = 1:size(tauD,2)
        alpha_s(:,j) = alpha(IX(:,j),j);
        end
        para_s.a = alpha_s ;
        para_s.tauD = tauD_s*1e6  ;
        para_s.D = wf^2/4./para_s.tauD*10^6 ;
        para_s.Ginf = para(5,:) ;
        
        tauD_std = abs(std(3:4,:)) ;
        alpha_std = abs(std(1:2,:)) ;            
        for j = 1:size(tauD,2)
        alpha_std_s(:,j) = alpha_std(IX(:,j),j);
        tauD_std_s(:,j) = tauD_std(IX(:,j),j);
        end       
        std_s.a = alpha_std_s ;
        std_s.tauD = tauD_std_s*1e6  ;
        std_s.Ginf = std(5,:) ;
        
    case '1comp2Ddiff+flow'                
        tauD_s = abs(para(2:3,:)) ;
        alpha_s = abs(para(1,:)) ;       
       
        para_s.a = alpha_s ;
        para_s.tauD = tauD_s*1e6  ;
        para_s.D(1,:) = wf^2/4./para_s.tauD(1,:)*10^6 ;
        para_s.D(2,:) = wf./para_s.tauD(2,:)*10^6 ; % flow velocity
        para_s.Ginf = para(4,:) ;
        
        tauD_std_s = abs(std(2:3,:)) ;
        alpha_std_s = abs(std(1,:)) ;            
        
        std_s.a = alpha_std_s ;
        std_s.tauD = tauD_std_s*1e6  ;
        std_s.Ginf = std(4,:) ;
        
    case '2comp2Ddiff+flow'
        tauD = abs(para(3:5,:)) ;
        alpha = abs(para(1:2,:)) ;       
        [tauD_s,IX] = sort(tauD(1:2,:),1,'ascend') ;       

        for j = 1:size(tauD,2)
        alpha_s(:,j) = alpha(IX(:,j),j);
        end
        para_s.a = alpha_s ;
        para_s.tauD = [tauD_s; tauD(3,:)]*1e6  ;
        para_s.D(1:2,:) = wf^2/4./para_s.tauD(1:2,:)*10^6 ;
        para_s.D(3,:) = wf./para_s.tauD(3,:)*10^6 ; % flow velocity
        para_s.Ginf = para(6,:) ;
        
        tauD_std = abs(std(3:4,:)) ;
        alpha_std = abs(std(1:2,:)) ;            
        for j = 1:size(tauD,2)
        alpha_std_s(:,j) = alpha_std(IX(:,j),j);
        tauD_std_s(:,j) = tauD_std(IX(:,j),j);
        end       
        std_s.a = alpha_std_s ;
        std_s.tauD = [tauD_std_s; abs(std(5,:))]*1e6  ;
        std_s.Ginf = std(6,:) ;
        
    case '1comp2Ddiff+1com_flow'                
        tauD_s = abs(para(3:4,:)) ;
        alpha_s = abs(para(1:2,:)) ;       
       
        para_s.a = alpha_s ;
        para_s.tauD = tauD_s*1e6  ;
        para_s.D(1,:) = wf^2/4./para_s.tauD(1,:)*10^6 ;
        para_s.D(2,:) = wf./para_s.tauD(2,:)*10^6 ; % flow velocity
        para_s.Ginf = para(5,:) ;
        
        tauD_std_s = abs(std(3:4,:)) ;
        alpha_std_s = abs(std(1:2,:)) ;            
        
        std_s.a = alpha_std_s ;
        std_s.tauD = tauD_std_s*1e6  ;
        std_s.Ginf = std(5,:) ;    
        
    case '2comp2Ddiff_im'
        D = abs(para(3:4,:)) ;
        alpha = abs(para(1:2,:)) ;        
       [D_s,IX] = sort(D,1,'descend') ;        

        for j = 1:size(D,2)
        alpha_s(:,j) = alpha(IX(:,j),j);
        end
        para_s.a = alpha_s ;
        para_s.D = D_s ;
%         para_s.tauD = wf^2/4./para_s.D*10^6 ;
        para_s.tauD = Aeff/4./para_s.D*10^6 ;
        para_s.Ginf = para(5,:) ;
        
        D_std = abs(std(3:4,:)) ;
        alpha_std = abs(std(1:2,:)) ;            
        for j = 1:size(D,2)
        alpha_std_s(:,j) = alpha_std(IX(:,j),j);
        D_std_s(:,j) = D_std(IX(:,j),j);
        end       
        std_s.a = alpha_std_s ;
        std_s.D = D_std_s  ;
        std_s.Ginf = std(5,:) ;
        
     case '2comp_2Ddiff+exp_im'        
        para_s.a = abs(para([1;4],:)) ;
        para_s.D = abs(para(2,:)) ;
        para_s.tauD = wf^2/4./para_s.D*10^6 ;
        para_s.Ginf = para(3,:) ;
        para_s.taut = abs(para(end,:))*1e6 ;
                
        std_s.a = abs(std([1;4],:)) ;
        std_s.D = abs(std(2,:)) ;
        std_s.Ginf = std(3,:) ;
        std_s.taut = abs(std(end,:))*1e6 ;                        
        
    case {'2comp2Ddiff_0Ginf' '2comp2Ddiff_0Ginf_im'}                
        tauD = abs(para(3:4,:)) ;
        alpha = abs(para(1:2,:)) ;     
        if strcmp(model, '2comp2Ddiff_0Ginf_im')
            [tauD_s,IX] = sort(tauD,1,'descend') ;
        else
           [tauD_s,IX] = sort(tauD,1,'ascend') ;
        end

        for j = 1:size(tauD,2)
        alpha_s(:,j) = alpha(IX(:,j),j);
        end
        para_s.a = alpha_s ;
        para_s.tauD = tauD_s*1e6  ;
        para_s.D = wf^2/4./para_s.tauD*10^6 ;
        
        tauD_std = abs(std(3:4,:)) ;
        alpha_std = abs(std(1:2,:)) ;            
        for j = 1:size(tauD,2)
        alpha_std_s(:,j) = alpha_std(IX(:,j),j);
        tauD_std_s(:,j) = tauD_std(IX(:,j),j);
        end       
        std_s.a = alpha_std_s ;
        std_s.tauD = tauD_std_s*1e6  ;
        
    case {'2comp2D+3Ddiff'}
        tauD = abs(para(3:4,:)) ;
        alpha = abs(para(1:2,:)) ;    
        para_s.a = alpha ;
        para_s.tauD = tauD*1e6  ;
        para_s.Ginf = para(5,:) ;
        
    case {'3comp3Ddiff' '3comp2Ddiff'}
        tauD = abs(para(4:6,:)) ;
        alpha = abs(para(1:3,:)) ;             
        [tauD_s,IX] = sort(tauD,1,'ascend') ;      
        for j = 1:size(tauD,2)
        alpha_s(:,j) = alpha(IX(:,j),j);
        end
        
        para_s.a = alpha_s ;
        para_s.tauD = tauD_s*1e6  ;
        para_s.D = wf^2/4./para_s.tauD*10^6 ;
        para_s.Ginf = para(7,:) ;
        
        tauD_std = abs(std(4:6,:)) ;
        alpha_std = abs(std(1:3,:)) ;            
        for j = 1:size(tauD,2)
        alpha_std_s(:,j) = alpha_std(IX(:,j),j);
        tauD_std_s(:,j) = tauD_std(IX(:,j),j);
        end       
        std_s.a = alpha_std_s ;
        std_s.tauD = tauD_std_s*1e6  ;
        std_s.Ginf = std(7,:) ;
        
    case  '3comp2Ddiff_im'
         D = abs(para(4:6,:)) ;
        alpha = abs(para(1:3,:)) ;             
        [D_s,IX] = sort(D,1,'descend') ;

        for j = 1:size(D,2)
        alpha_s(:,j) = alpha(IX(:,j),j);
        end
        
        para_s.a = alpha_s ;
        para_s.D = D_s  ;
%         para_s.tauD = wf^2/4./para_s.D*10^6 ;
        para_s.tauD = Aeff/4./para_s.D*10^6 ;
        para_s.Ginf = para(7,:) ;
        
        D_std = abs(std(4:6,:)) ;
        alpha_std = abs(std(1:3,:)) ;            
        for j = 1:size(D,2)
        alpha_std_s(:,j) = alpha_std(IX(:,j),j);
        D_std_s(:,j) = D_std(IX(:,j),j);
        end       
        std_s.a = alpha_std_s ;
        std_s.D = D_std_s  ;
        std_s.Ginf = std(7,:) ;
        
    case {'3comp2Ddiff_0Ginf' '3comp2Ddiff_0Ginf_im'}
        tauD = abs(para(4:6,:)) ;
        alpha = abs(para(1:3,:)) ;     
        if strcmp(model, '3comp2Ddiff_0Ginf_im')
            [tauD_s,IX] = sort(tauD,1,'descend') ;
        else
           [tauD_s,IX] = sort(tauD,1,'ascend') ;
        end

        for j = 1:size(tauD,2)
        alpha_s(:,j) = alpha(IX(:,j),j);
        end
        
        para_s.a = alpha_s ;
        para_s.tauD = tauD_s*1e6  ;
        para_s.D = wf^2/4./para_s.tauD*10^6 ;
        
        tauD_std = abs(std(4:6,:)) ;
        alpha_std = abs(std(1:3,:)) ;            
        for j = 1:size(tauD,2)
        alpha_std_s(:,j) = alpha_std(IX(:,j),j);
        tauD_std_s(:,j) = tauD_std(IX(:,j),j);
        end       
        std_s.a = alpha_std_s ;
        std_s.tauD = tauD_std_s*1e6  ;
        
                
     case {'3comp2_2D+1_3Ddiff'}
        para(1:6,:) = abs(para(1:6,:)) ;
         tauD = abs(para(4:5,:)) ;
        alpha = abs(para(1:2,:)) ;     
       [tauD_s,IX] = sort(tauD,1,'ascend') ;

        for j = 1:size(tauD,2)
        alpha_s(:,j) = alpha(IX(:,j),j);
        end
        
        para_s.a = [alpha_s; para(3,:)];
        para_s.tauD = [tauD_s; para(6,:)]*1e6  ;
        para_s.Ginf = para(7,:) ;
        
    case {'3comp1_2D+2_3Ddiff'}
        para(1:6,:) = abs(para(1:6,:)) ;
        tauD = abs(para(5:6,:)) ;
        alpha = abs(para(2:3,:)) ;     
       [tauD_s,IX] = sort(tauD,1,'ascend') ;

        for j = 1:size(tauD,2)
        alpha_s(:,j) = alpha(IX(:,j),j);
        end
        
        para_s.a = [para(1,:); alpha_s];
        para_s.tauD = [para(4,:); tauD_s]*1e6  ;
        para_s.Ginf = para(7,:) ;
        
    case {'4comp3Ddiff' '4comp2Ddiff'}
        tauD = abs(para(5:8,:)) ;
        alpha = abs(para(1:4,:)) ;     
       [tauD_s,IX] = sort(tauD,1,'ascend') ;

        for j = 1:size(tauD,2)
        alpha_s(:,j) = alpha(IX(:,j),j);
        end
        
        para_s.a = alpha_s ;
        para_s.tauD = tauD_s*1e6  ;
        para_s.Ginf = para(9,:) ;
        
        tauD_std = abs(std(5:8,:)) ;
        alpha_std = abs(std(1:4,:)) ;            
        for j = 1:size(tauD,2)
        alpha_std_s(:,j) = alpha_std(IX(:,j),j);
        tauD_std_s(:,j) = tauD_std(IX(:,j),j);
        end       
        std_s.a = alpha_std_s ;
        std_s.tauD = tauD_std_s*1e6  ;
        std_s.Ginf = std(9,:) ;
        
     case {'anomal3Ddiff'}           
        para(2,:) = abs(para(2,:))*1e6 ;
        para_s.a = para(1,:) ;
        para_s.tauD = para(2,:) ;
        para_s.Ginf = para(4,:) ;
        para_s.exp = para(3,:) ;
        
    case {'trip1comp3Ddiff' 'trip1comp3Ddiff_1comfixed' 'trip1comp3Ddiff_tripfixed'}  
        para(2,:) = abs(para(2,:))*1e6 ;
        
        para_s.a = para(1,:) ;
        para_s.tauD = para(2,:) ;
        para_s.Ginf = para(3,:) ;
        para_s.ft = abs(para(end-1,:)) ;
        para_s.taut = abs(para(end,:))*1e6 ;
        
        std(2,:) = abs(std(2,:))*1e6 ;
        std_s.a = abs(std(1,:)) ;
        std_s.tauD = std(2,:) ;
        std_s.Ginf = std(3,:) ;
        std_s.ft = abs(std(end-1,:)) ;
        std_s.taut = abs(std(end,:))*1e6 ;
        
    case {'trip2comp3Ddiff' 'trip2comp3Ddiff_1comfixed' 'trip2comp3Ddiff_tripfixed'} 
        tauD = abs(para(3:4,:)) ;
        alpha = abs(para(1:2,:)) ;     
        [tauD_s,IX] = sort(tauD,1,'ascend') ;
        for j = 1:size(tauD,2)
        alpha_s(:,j) = alpha(IX(:,j),j);
        end        
        para_s.a = alpha_s ;
        para_s.tauD = tauD_s*1e6  ;
        para_s.Ginf = para(5,:) ;
        para_s.ft = abs(para(end-1,:)) ;
        para_s.taut = abs(para(end,:))*1e6 ;
        
        tauD_std = abs(std(3:4,:)) ;
        alpha_std = abs(std(1:2,:)) ;            
        for j = 1:size(tauD,2)
        alpha_std_s(:,j) = alpha_std(IX(:,j),j);
        tauD_std_s(:,j) = tauD_std(IX(:,j),j);
        end       
        std_s.a = alpha_std_s ;
        std_s.tauD = tauD_std_s*1e6  ;
        std_s.Ginf = std(5,:) ;
        std_s.ft = abs(std(end-1,:)) ;
        std_s.taut = abs(std(end,:))*1e6 ;
        
    case {'trip3comp3Ddiff' 'trip3comp3Ddiff_1comfixed' 'trip3comp3Ddiff_tripfixed'}
        tauD = abs(para(4:6,:)) ;
        alpha = abs(para(1:3,:)) ;     
       [tauD_s,IX] = sort(tauD,1,'ascend') ;

        for j = 1:size(tauD,2)
        alpha_s(:,j) = alpha(IX(:,j),j);
        end
        
        para_s.a = alpha_s ;
        para_s.tauD = tauD_s*1e6  ;
        para_s.Ginf = para(7,:) ;
        para_s.ft = abs(para(end-1,:)) ;
        para_s.taut = abs(para(end,:))*1e6 ;
        
        tauD_std = abs(std(4:6,:)) ;
        alpha_std = abs(std(1:3,:)) ;            
        for j = 1:size(tauD,2)
        alpha_std_s(:,j) = alpha_std(IX(:,j),j);
        tauD_std_s(:,j) = tauD_std(IX(:,j),j);
        end       
        std_s.a = alpha_std_s ;
        std_s.tauD = tauD_std_s*1e6  ;
        std_s.Ginf = std(7,:) ;
        std_s.ft = abs(std(end-1,:)) ;
        std_s.taut = abs(std(end,:))*1e6 ;
        
    case {'trip4comp3Ddiff'}
        tauD = abs(para(5:8,:)) ;
        alpha = abs(para(1:4,:)) ;     
       [tauD_s,IX] = sort(tauD,1,'ascend') ;

        for j = 1:size(tauD,2)
        alpha_s(:,j) = alpha(IX(:,j),j);
        end
        
        para_s.a = alpha_s ;
        para_s.tauD = tauD_s*1e6  ;
        para_s.Ginf = para(9,:) ;
        para_s.ft = abs(para(end-1,:)) ;
        para_s.taut = abs(para(end,:))*1e6 ;
        
        tauD_std = abs(std(5:8,:)) ;
        alpha_std = abs(std(1:4,:)) ;            
        for j = 1:size(tauD,2)
        alpha_std_s(:,j) = alpha_std(IX(:,j),j);
        tauD_std_s(:,j) = tauD_std(IX(:,j),j);
        end       
        std_s.a = alpha_std_s ;
        std_s.tauD = tauD_std_s*1e6  ;
        std_s.Ginf = std(9,:) ;
        std_s.ft = abs(std(end-1,:)) ;
        std_s.taut = abs(std(end,:))*1e6 ;
        
        case 'trip1comp3Drotdiff' % 1-component 3D diffusion + triplet blinking
        para(2,:) = abs(para(2,:))*1e6 ;
        
        para_s.a = para(1,:) ;
        para_s.tauD = para(2,:) ;
        para_s.Ginf = para(3,:) ;
        para_s.ft = abs(para(4,:)) ;
        para_s.taut = abs(para(5,:))*1e6 ;
        para_s.l =  abs(para(6,:));
        para_s.Dr = abs(para(7,:));
        
        std(2,:) = abs(std(2,:))*1e6 ;
        std_s.a = abs(std(1,:)) ;
        std_s.tauD = std(2,:) ;
        std_s.Ginf = std(3,:) ;
        std_s.ft = abs(std(4,:)) ;
        std_s.taut = abs(std(5,:))*1e6 ;
        para_s.l =  abs(std(6,:));
        para_s.tauR = abs(std(7,:));
    otherwise
        error('models are not on the list')
end

end        
