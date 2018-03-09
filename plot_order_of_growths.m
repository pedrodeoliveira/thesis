
% close all

% Problem size
n=0:0.1:9;
N=2.^n;
T = 1e-2;


color = [1 1 1];
fontsize = 14;
h = figure;
set(h,'Color',color)

axesh = axes('Parent',h);
set(axesh,'FontSize',fontsize);

plot(N,T*N.^2,'LineWidth',1.5)
hold on
plot(N,T*N.*log2(N),'r--','LineWidth',1.5)

h_legend = legend(axesh,'quadratic','linearithmetic');
set(h_legend,'FontSize',fontsize)
xlabel('Problem Size N','FontSize',fontsize)
ylabel('Time','FontSize',fontsize);
axis([0 500 0 500])
set(axesh,'YTick',0:100:500)
set(axesh,'YTickLabel',{'0','100T','200T','300T','400T','500T',})

% saveas(h,'figures/order_of_growth.eps')
