function bar_minmax(minmax,location,color)

a = minmax(:,1);
b = minmax(:,2);

m = mean([a b],2);
errorbar(m,location,m-a,b-m,'.','horizontal','LineWidth',2,'CapSize',4, 'Color', color)

end