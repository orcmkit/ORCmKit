function figuresize( w , h , u )
%FIGURESIZE Set a figure to a specific size
%
% When saving a figure as a PDF, it is necessary to set the
% figure size appropriately. This function sets the "paper size"
% sufficient that the figure is saved with a tight bounding box.
% It will also set the figure size on screen correspondingly.
%
% figuresize(width,height) - sets the figure size in centimeters
% figuresize(width,height,units) - sets the figure size in <units>
%
% <units> can be any of the standard Matlab lengths.
%
% Will Robertson
% 28 April 2010
%
% Copyright and licence information appended.

p = inputParser;
p.addRequired('width', @(x) isnumeric(x) && all(size(x)==1) );
p.addRequired('height',@(x) isnumeric(x) && all(size(x)==1) );
p.addOptional('units','centimeters',...
  @(x) any(strcmpi(x,{'normalized','centimeters','inches','points'})) );

p.parse( w, h, u );
w = p.Results.width;
h = p.Results.height;
u = p.Results.units;

p = 0.01;

set(gcf,'Units',u);
screenpos = get(gcf,'Position');

set(gcf,...
  'Position',[screenpos(1:2) w h],...
  'PaperUnits',u,...
  'PaperPosition',[p*w p*h w h],...
  'PaperSize',[w*(1+2*p) h*(1+2*p)]);

end