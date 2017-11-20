% This script runs the Bayesian analysis on the simulated FCS multiple curves data 
% sets of two-component with vaying D2 for each data set. Each data set has
% 64 individual curves. See function "FCS_bayes_analysis_multi_curves" for 
% more info about input and output. This file was tested to run properly on MATLAB
% R2011a with Statistics Toolbox and Image Processing Toolbox installed.
% -----------------------------------------------------------------
% Copyright MIT 2012
% Developed by Syuan-Ming Guo
% Laboratory for Computational Biology & Biophysics
% Apr 23, 2012
% -----------------------------------------------------------------

addpath('.\functions')
load('.\examples\FCS_64cur_vary_Dr.mat') ;
model = {'1comp3Ddiff' '2comp3Ddiff' '3comp3Ddiff'};
D2 = 10.^linspace(0.8,2.8,15) ;  % simulation parameter for D2
D2 = D2(end:-1:1) ;
%% Compute model probabilities
[mp fit_para] = FCS_bayes_analysis_multi_curves...
    (t, corrFCS_cell, model, ks, 0, 8, 'GLS', '2way') ;

[mp_m err_l err_h] = mp_stac(mp) ;  % get statistics of model probabilities

%% plot model probabilities
Dr = D2./63.1; % simulation parameter for D-ratio
x1 = Dr ;
l_style = {'rs-','bo-', 'k^-', 'g*-', 'cv-'} ;
l_color = {'r','b', 'k', 'g', 'c'};
n_model = numel(model) ;
figure(11)
for k =1:n_model
errorbar(x1, mp_m(k,:), err_l(k,:),err_h(k,:),...
    l_style{k},'MarkerFaceColor',l_color{k},'markersize', 10,'linewidth', 2)
hold on
end
set(gca,'xscale','log');
tick = 10.^[-7:3] ;
set(gca, 'xtick',tick)
axis tight
set(gca, 'ylim', [0 1]) ;
legend('N_D=1','N_D=2','N_D=3')
legend('location', 'east')
xlabel('D_2/D_1')
ylabel('P(M)','FontSize',15)
format_fig2(5)
hold off
