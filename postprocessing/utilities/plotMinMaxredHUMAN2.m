function plotMinMaxredHUMAN2(rxns,minmax1,minmax2,title_name)

n = length(rxns);
bar_minmax(minmax1,1:3:3*n,[0, 0, 0])
hold on
if contains(title_name, 'Recon 2')
bar_minmax(minmax2,1:3:3*n,[0.93, 0.57, 0.13])
% bar_minmax(minmax3,1:4:4*n,[0, 0, 0])
else
    bar_minmax(minmax2,1:3:3*n,[0.29, 0.59, 0.82])
end

hold off

yticks(1:3:3*n);
yticklabels(rxns);
xlabel('rate [mmol/gDW/h]')
set(gca,'FontSize',14)
%xlim([-100 100])
title(title_name,'FontSize',16)
end

%orange[0.93, 0.57, 0.13]
%yellow[0.98, 0.85, 0.5]
%darkgreen[0.09, 0.45, 0.27])
%pink[0.98, 0.38, 0.5]
%blue[0.29, 0.59, 0.82]
%purple[0.6, 0.4, 0.8]
%darkblue[0.0, 0.13, 0.28]
%red[0.93, 0.11, 0.14]
%grey[0.5, 0.5, 0.5]
