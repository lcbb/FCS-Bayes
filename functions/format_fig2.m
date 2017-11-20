% figure export settings (default is opt=2)
% specify the font size (), and the figure size (cm)
% -----------------------------------------------------------------
% Copyright MIT 2012
% Developed by Syuan-Ming Guo
% Laboratory for Computational Biology & Biophysics
% Apr 06, 2012
% -----------------------------------------------------------------
function format_fig2(varargin)

figsizeX = 25 ;
figsizeY = 15 ;

if isempty(varargin) 
    opt = 1;
else
    opt=varargin{1};
end
if ~any(opt==[1 2 3 4 5 6 7 8 9])
    opt = 2;
end  

if opt == 1;
    fontsize = 18;
    fontsize_leg = 14;
    figsizeX = 20 ;
    figsizeY = 15 ;  
elseif opt == 2
    fontsize = 24;
    fontsize_leg = 22;
elseif opt == 3
    fontsize = 30;
    fontsize_leg = 26;
elseif opt == 4
    fontsize = 24;
    fontsize_leg = 22;
    figsizeX = 20 ;
    figsizeY = 20 ;
elseif opt == 5
    fontsize = 24;
    fontsize_leg = 22;
    figsizeX = 20 ;
    figsizeY = 15 ;  
elseif opt == 6
    fontsize = 22;
    fontsize_leg = 18;
    figsizeX = 30 ;
    figsizeY = 15 ;
elseif opt == 7
    fontsize = 18;
    fontsize_leg = 14;
    figsizeX = 25 ;
    figsizeY = 20 ;
elseif opt == 8
    fontsize = 24;
    fontsize_leg = 22;
    figsizeX = 20 ;
    figsizeY = 25 ;
elseif opt == 9
    fontsize = 22;
    fontsize_leg = 18;
    figsizeX = 20 ;
    figsizeY = 15 ;  
end


% legend font size
legend_handle =legend();
set(legend_handle,'FontSize', fontsize_leg)% legend font size

% axis label size
set(get(gca,'XLabel'),'FontSize',fontsize)
set(get(gca,'YLabel'),'FontSize',fontsize)

% axis font size
set(gca,'FontSize', fontsize)



pos = [10 5 figsizeX figsizeY];% figure position and size
% % export paper properties
% set(gcf,'PaperType','<custom>','unit','centimeters','PaperSize',[figsizeX figsizeY],...
%     'InvertHardcopy','off',...
%     'Color',[1 1 1]);
set(gcf,'unit','centimeters','position',pos) % figure size
set(gcf,'PaperPositionMode','auto')
% 
% % strech axis to fill window
% % TI = get(gca,'TightInset');  OP = get(gca,'OuterPosition');  Pos = OP + [ TI(1:2), -TI(1:2)-TI(3:4) ];  set( gca,'Position',Pos);
% 
% % set legend position to 'Best'
% set(legend_handle,'Location','Best')% legend font size
box on

