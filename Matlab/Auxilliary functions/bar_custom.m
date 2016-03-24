function h = bar_custom(Y)
h=bar(Y);
yb = cat(1, h.YData); 
xb = bsxfun(@plus, h(1).XData, [h.XOffset]');
hold on;
for i = 1:size(yb,2)
    for j = 1:length(yb(:,1))
        
        text(xb(j, i),yb(j, i), cellstr(num2str(Y(i, j),3)), 'rotation', 90);
        
    end
end
hold off
end