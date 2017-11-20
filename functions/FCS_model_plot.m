function h = FCS_model_plot(a,t_extra,k,model,L,sigma, varargin)
switch model
    case 'null'
        f = @(t)null_model(a,t) ;
        
    case '1comp3Ddiff'
        f = @(t)diff1com3D(a,t,k) ;
    
    case {'2comp3Ddiff' '2comp3Ddiff_1comfixed'}
        f = @(t)diff2com3D(a,t,k) ;
    
    case {'3comp3Ddiff' '3comp3Ddiff_1comfixed'}
        f = @(t)diff3com3D(a,t,k) ;
    
    case '4comp3Ddiff'
        f = @(t)diff4com3D(a,t,k) ;
        
    case 'flow'
        f = @(t)flow1com(a,t) ;
              
    case '1comp3Ddiff+flow'
        f = @(t)diff_flow1com3D(a,t,k) ;
            
    case '2comp3Ddiff+flow'    
        f = @(t)diff_flow2com3D(a,t,k) ;
        
    case 'anomal3Ddiff'
        f = @(t)diff_anomal3D(a,t,k) ;
        
    case '1comp2Ddiff'
        f = @(t)diff1com2D(a,t) ;
    
    case '2comp2Ddiff'
        f = @(t)diff2com2D(a,t) ;

    case '3comp2Ddiff'
        f = @(t)diff3com2D(a,t) ;
        
    case '1comp2Ddiff_0Ginf'
        f = @(t)diff1com2D_0Ginf(a,t) ;
    
    case '2comp2Ddiff_0Ginf'
        f = @(t)diff2com2D_0Ginf(a,t) ;

    case '3comp2Ddiff_0Ginf'
        f = @(t)diff3com2D_0Ginf(a,t) ;
        
    case '1comp2Ddiff+flow'
        f = @(t)diff_flow1com2D(a,t) ;
    
    case '2comp2Ddiff+flow'    
        f = @(t)diff_flow2com2D(a,t) ;
        
    case '1comp2Ddiff+1com_flow'    
        f = @(t)diff1com_flow1com2D(a,t) ;
    
    case '2comp2D+3Ddiff'
        f = @(t)diff2com2D_3D(a,t,k) ;
    
    case '3comp2_2D+1_3Ddiff'
        f = @(t)diff3com_2_2D_1_3D(a,t,k) ;
    
    case '3comp1_2D+2_3Ddiff'
        f = @(t)diff3com_1_2D_2_3D(a,t,k) ;
    
    case {'trip1comp3Ddiff' 'trip1comp3Ddiff_1comfixed' 'trip1comp3Ddiff_tripfixed'}
        f = @(t)trip_diff1com3D(a,t,k) ;
    
    case {'trip2comp3Ddiff' 'trip2comp3Ddiff_1comfixed' 'trip2comp3Ddiff_tripfixed'}
        f = @(t)trip_diff2com3D(a,t,k) ;
    
    case {'trip3comp3Ddiff' 'trip3comp3Ddiff_1comfixed' 'trip3comp3Ddiff_tripfixed'}
        f = @(t)trip_diff3com3D(a,t,k) ;
    
    case 'trip4comp3Ddiff'
        f = @(t)trip_diff4com3D(a,t,k) ;   
    case '1comp3Ddiff_fk'
        f = @(t)diff1com3D_free_k(a,t) ;
    case '2comp3Ddiff_fk'
        f = @(t)diff2com3D_free_k(a,t) ;
    case '3comp3Ddiff_fk'
        f = @(t)diff3com3D_free_k(a,t) ;
    case '1comp3Ddiff_log'
        f = @(t)diff1com3D_log(a,t,k) ;
    case '2comp3Ddiff_log'
        f = @(t)diff2com3D_log(a,t,k) ;
    case '3comp3Ddiff_log'
        f = @(t)diff3com3D_log(a,t,k) ;
    case '1comp2Ddiff_log'
        f = @(t)diff1com2D_log(a,t) ;
    case '2comp2Ddiff_log'
        f = @(t)diff2com2D_log(a,t) ;
    case '3comp2Ddiff_log'
        f = @(t)diff3com2D_log(a,t) ;
    case '1comp2Ddiff_im'    
        f = @(t)diff1com2D_im(a,t,L,sigma) ;    
    case '2comp2Ddiff_im'    
        f = @(t)diff2com2D_im(a,t,L,sigma) ;
    case '3comp2Ddiff_im'    
        f = @(t)diff3com2D_im(a,t,L,sigma) ;
    case '1comp2Ddiff_0Ginf_im'    
        f = @(t)diff1com2D_0Ginf_im(a,t,L,sigma) ;    
    case '2comp2Ddiff_0Ginf_im'    
        f = @(t)diff2com2D_0Ginf_im(a,t,L,sigma) ;
    case '3comp2Ddiff_0Ginf_im'    
        f = @(t)diff3com2D_0Ginf_im(a,t,L,sigma) ;         
    case '2comp_2Ddiff+exp_im'
        f = @(t)diff1com2D_exp1com_im(a,t,L,sigma) ;
    case '1comp_exp'    
        f = @(t)exp1com(a,t) ;    
    case '2comp_exp'    
        f = @(t)exp2com(a,t) ;
        
end
if ~isempty(varargin)
    L = chol(varargin{1}, 'lower') ;
    h = plot(t_extra ,L\f(t_extra), 'linewidth',1) ;   
else
   h = plot(t_extra ,f(t_extra), 'linewidth',1) ;
end   
end