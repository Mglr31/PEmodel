function f=plot_audible_range(results)
%plot_audible_range(results)

colors=turbo;
m=length(results);
f=figure;
for i=1:m
    col=colors(round(i/m*256),:);
    plot(results(i).xr,results(i).yr,'LineWidth',2,'Color',col,'DisplayName',results(i).Name)
    hold on
    plot(-results(i).xl,results(i).yl,'LineWidth',2,'Color',col,'HandleVisibility','off')
    hold on
    plot(0,results(i).Zsource,'*','Color',col,'MarkerSize',30,'LineWidth',2,'HandleVisibility','off')


end
xlabel("Range (m)")
ylabel("Altitude (m)")
grid on
set(gca,'FontSize',22)
legend show

end