function [legend_h,object_h,plot_h,text_strings] = legalpha(varargin)
%LEGALPHA creates a legend whose patch objects match the transparency 
%(alpha values) of plotted patch objects. Syntax follows legend exactly. 
%In fact, if you so desire, you could always use legalpha instead of legend. 
% 
%
%% Example: 
% 
% plot(1:10,1:10,'b','linewidth',2); 
% hold on
% plot(7:.01:9,(1:.01:3).^2,'k','linewidth',2)
% plot(9,8,'rp','markersize',12)
% fill([6 7 7 6],[7 7 8 8],'blue')
% fill([7 9 9 7],[2 2 3 3],'blue','facealpha',0.4)
% plot(2:2:6,[4 4 4],'ks','markersize',12)
% fill([4 8 8 4],[1.5 1.5 5 5],'g','facealpha',0.2)
% box off
% 
% legalpha('blue line','black curve','red star','opaque blue box',...
%     'transparent blue box','black squares','transparent green','location','northwest')
% legend boxoff
% 
%%
% Written by Chad A. Greene of the University of Texas Instititue for
% Geophysics in August of 2014.  
% 
% See also legend. 


% Get handles of all patch objects in current axis: 
hla = findobj(gca,'type','patch'); 

% Get alpha properties of the current plotted patches: 
hlaFA = get(hla,'FaceAlpha');

% Create a legend: 
[legend_h,object_h,plot_h,text_strings] = legend(varargin{:});

% Get handles of patch objects in the legend: 
hlp = findobj(legend_h,'type','patch'); 

% Set alpha values: 
for k = 1:length(hlp)
set(hlp(k),'FaceAlpha',hlaFA{k}); 

end

% Clean up: 
if nargout == 0
    clear legend_h object_h plot_h text_strings
end