% This function generates the model probability map, the map and histogram 
% of diffusion coefficients, the map of blocking results (0--pass 1--fail) 
% ----input------------------------------------------------------
% imm: image of mean intensity 
% fit_para: 2D, m-by-6 cell array of fitted parameters of m model fits
%           - column 1: parameter estimates of model fits (see model
%           functions for details)
%           - column 2: fitting uncertainties of parameters
%           - column 3: log model probabilities
%           - column 4: residuals of fits
%           - column 5: reduced chi-square values
%           - column 6: corresponding p-values of the chi-square value
% block_min_ind: 2D matrix of the index of the fixed points for the blocking
% curve at each pixel
% o: parameter structure
% ----output------------------------------------------------------
% mp: 3D array of the model probability map (pixel, pixel, model)
% D_map: 1D cell array with each cell recording the fit diffsion coefficient 
% maps from each model. The diffsion coefficient map is a 3D array where each 
% 2D array is the map of each diffusion coefficient. The diffuison coefficient
% is set to NaN is the corresponding model probability <= 0.5
% figures:
%   931-model probability map: the probability of each model is represented
%   by the fraction of the color.
%   932-blocking map: binary map of the blocking results (0--pass 1--fail)  
%   36-mean intensity image with dots indicating pixels where two or three
%   diffusing components are detected (model probability >0.5) 
%   102-histogram of diffusion coefficients from preferred models
%   
% -----------------------------------------------------------------
% Syuan-Ming Guo
% sguo@mit.edu, Oct. 31, 2013
% -----------------------------------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [mp D_map] = model_likelyhood_map(imm, fit_para, block_min_ind, o) 
close all

%% pre-allocation and parameters
n_condi = size(fit_para,4) ;
n_model = size(fit_para,1) ;
n_cur = size(fit_para,3) ;
mp = zeros(n_condi, n_cur, n_model) ;
log_mp = zeros(n_condi, n_cur, n_model) ;
beta_c = (2*2.311*o.psf_sigma_um(1)/(o.binning*o.um_per_px)+0.06765)^-2+1 ; % correction of finite pixel size effect
wf = 2* o.psf_sigma_um(1)* beta_c ;
L = o.um_per_px*o.binning ;
m0 = L/(2*o.psf_sigma_um(1));  
Aeff = L^2/(erf(m0) + 1/sqrt(pi)/m0*(exp(-m0^2)-1))^2 ; % effective area
w1 = L/erf(m0/sqrt(2)) ;
gammaf = w1^2/Aeff ;
W0 = 2*o.psf_sigma_um(1)*sqrt(pi/2)*erf(m0/sqrt(2)) ; 
mp_thd = 0.5 ; % threshold for model probabilities
block_bi = block_min_ind==0 ;
model = o.model ;
%% Set NaN probablities to be zero  
for j = 1:n_condi      
    for k=1:n_cur
          logML = cell2mat(fit_para(:,3,k,j));
             for g = 1:numel(logML)
                if isnan(logML(g))||~isreal(logML(g))
                    logML(g)= -Inf ;
                end
             end
        PM = exp(logML-max(logML)) ;
        PM = PM./sum(PM) ;
        
        for m = 1:n_model
        mp(j,k,m) = PM(m) ;
        log_mp (j,k,m)= logML(m) ;
        end                
    end
end   

%% parameter map
tauD_map = cell(n_model, 1) ;
D_map = cell(n_model, 1) ;
% for j = 1:n_model
%     tauD_map{j} = zeros(n_condi, n_cur, numel(fit_para{j,1,1,1})) ;
%     D_map{j} = zeros(n_condi, n_cur, numel(fit_para{j,1,1,1})) ;
% end

B_map = cell(n_model, 1) ;
B_map{1} = zeros(n_condi, n_cur, 1) ;
B_map{2} = zeros(n_condi, n_cur, 1) ;
N_map = cell(n_model, 1) ;
N_map{1} = zeros(n_condi, n_cur, 1) ;
N_map{2} = zeros(n_condi, n_cur, 1) ;

for k = 1:n_model
        for j = 1:n_condi    
        clear para_N
        para = squeeze(cell2mat(fit_para(k,1,:,j)));
        para_std = squeeze(cell2mat(fit_para(k,2,:,j)));
        if size(para,1) == 1
            para = para';
            para_std = para_std';
        end        
        [para_s para_std_s]= sort_fitted_para_2(para, para_std, model{k}, Aeff, wf) ;
        
        % pre-allocate the map matrices
        if j==1
            tauD_map{k} = zeros(n_condi, n_cur, size(para_s.tauD,1)) ;
            D_map{k} = zeros(n_condi, n_cur, size(para_s.tauD,1)) ;
        else
        end
        
        tauD_map{k}(j,:,:) = permute(para_s.tauD, [3 2 1]) ;
        D_map{k}(j,:,:) = permute(para_s.D, [3 2 1]) ;        
        N_map{k}(j,:) = 1./sum(para_s.a,1) ;                    
        end
end
%% MP map
figure(931)
imshow(mp, 'InitialMagnification', 'fit')
axis on
colormap(jet)        
set(gca, 'fontsize',15)
caxis auto
hold on 
plot(0,0, 'rs', 'markerfacecolor', 'r')
plot(0,0, 'gs', 'markerfacecolor', 'g')
plot(0,0, 'bs', 'markerfacecolor', 'b')
legend('N_D=1','N_D=2','N_D=3')
legend('location', 'Eastoutside')
hold off
title('model probabilities','FontSize',25)
format_fig2(2)

%% set D with MP <= 0.5 to NaN
for k = 1:n_model
mp_f = mp(:,:,k)< mp_thd ;
mp_f_r = repmat(mp_f, [1 1 size(tauD_map{k},3)]) ;
tauD_map{k}(mp_f_r) = NaN ;
D_map{k}(mp_f_r) = NaN ;
end
%% blocking filter
figure(932)
imshow(block_bi, 'InitialMagnification', 'fit')
axis on
colormap(gray)    
set(gca, 'fontsize',15)
caxis auto
hold on 
title('blocking binary map','FontSize',25)
format_fig2(2)
  
%% intensity map with pb correction
figure(36)    
imshow(imm, 'InitialMagnification', 'fit')    
axis on
h = colorbar ;
caxis auto
title('mean intensity','FontSize',25)
format_fig2(2)     
%% overlay 2-component MP on the I map
color_1 = [86 180 233]./255 ;
color_2 = [213 94 0]./255 ;
color_3 = [0 158 115]./255 ;
mrkr_size = 200/max(size(mp,1), size(mp,2));
[mask_y mask_x] = find((mp(:,:,2)>mp_thd).*(~block_bi)) ;
hold on
plot(mask_x, mask_y,'o','markersize',mrkr_size, 'markerfacecolor', color_2,...
    'markeredgecolor', color_2,'LineWidth', 2.5); % note the inversed order of row & column index    
[mask_y mask_x] = find((mp(:,:,3)>mp_thd).*(~block_bi)) ;
plot(mask_x, mask_y,'s','markersize',mrkr_size, 'markerfacecolor', color_3,...
    'markeredgecolor', color_3, 'LineWidth', 2.5); % note the inversed order of row & column index
[mask_y mask_x] = find((mp(:,:,2)>mp_thd).*block_bi) ;
plot(mask_x, mask_y,'o','markersize',mrkr_size, 'markeredgecolor', color_2,'LineWidth', 2.5); % note the inversed order of row & column index    
[mask_y mask_x] = find((mp(:,:,3)>mp_thd).*block_bi) ;
plot(mask_x, mask_y,'s','markersize',mrkr_size,'markeredgecolor', color_3, 'LineWidth', 2.5); % note the inversed order of row & column index
legend(sprintf('N_D=2,\nBlocking passed'),sprintf('N_D=2,\nBlocking failed'),sprintf('N_D=3,\nBlocking passed'),sprintf('N_D=3,\nBlocking failed'))
legend('location', 'Eastoutside')
%% histogram of D
bin_l = linspace(-4,4,60) ;
[count1] = hist(log10(D_map{1}(:)),bin_l) ;
[count2] = hist(log10(D_map{2}(:)),bin_l) ;
[count3] = hist(log10(D_map{3}(:)),bin_l) ;

count1_n = count1./numel(D_map{1}) ;
count2_n = count2./numel(D_map{1});
count3_n = count3./numel(D_map{1}) ;

figure(102)
p0 = bar(bin_l, count1_n,'hist') ;
hold on
p1 = bar(bin_l, count2_n,'hist') ;
p2 = bar(bin_l, count3_n,'hist') ;

sh=findall(gcf,'marker','*');
     delete(sh);
set(p0, 'facecolor', color_1)
set(p1, 'facecolor', color_2)
set(p2, 'facecolor', color_3)
alpha(p0,0.5)
alpha(p1,0.5)
alpha(p2,0.5)
set(gca, 'xlim', [-2 2])
% set(gca, 'ylim', [0 0.2])
set(gca,'XTick',[-2:1:2])
a_pos =  get(gca, 'position');
a_pos(1) = a_pos(1)+0.05;
a_pos(4) = a_pos(4)*0.9;
a_pos(2) = a_pos(2)+ a_pos(4)*0.1;
set(gca, 'position', a_pos);
tx = format_ticks(gca,{'10^{-2}' '10^{-1}' '10^0','10^1','10^2'},[],[],[],[],[],[]) ;
set(tx,'units', 'normalized'); 
tx_pos = get(tx, 'position') ;
tx_pos_mid = tx_pos{3} ;
set(tx,'fontsize', 24); 
xl = xlabel('D (\mum^2/s)') ; 
xl_pos = get(xl, 'position') ;
xl_pos(2) = -0.08;

set(xl, 'position', xl_pos) ;
ylabel('Normalized count')
l_h = legend('N_D=1','N_D=2', 'N_D=3') ;
k = findobj(get(l_h, 'children'), 'type', 'patch') ;
set(k, 'FaceAlpha', 0.5)
format_fig2(5)
hold off

end

