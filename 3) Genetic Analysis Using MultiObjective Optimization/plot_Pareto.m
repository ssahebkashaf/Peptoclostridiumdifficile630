
function plot_Pareto(num_file)

load(strcat('solution',num2str(num_file),'.mat'));
num_col = size(chromosome,2);


x = -chromosome(:,num_col-2);
y = -chromosome(:,num_col-3);

x = sort(x, 'ascend');
y = sort(y, 'descend');

plot(x,y,'*:')
axis([0 100 0 20])
grid on

%set(gca,'XTick',-pi:pi/2:pi)
%set(gca,'XTickLabel',{'-pi','-pi/2','0','pi/2','pi'})

title('');
xlabel('Biomass [h^{-1}]');
ylabel('Acetate [mmolh^{-1}gDW^{-1}]');