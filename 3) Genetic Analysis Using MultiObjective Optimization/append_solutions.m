function matrice = append_solutions(num_files, step)

primo_file=load('solution1.mat');
dim = size(primo_file.chromosome);
num_col=dim(2);
primo=num_col-2;
secondo=num_col-3;

%matrice=[];
%matrice=zeros(dim(1)*num_files, dim(2));

ultimo_file = load(strcat('solution',num2str(num_files),'.mat'));
chromosome = ultimo_file.chromosome;
    
x = -chromosome(:,num_col-2);
y = -chromosome(:,num_col-3);


    A = [x y];
    A_sort = sortrows(A,1);
    x = A_sort(:,1);
    y = A_sort(:,2);

   

fig1=figure(1);
plot(x,y,'*:', 'Color',[1 0 0])
grid on

title('');
xlabel('Biomass [h^{-1}]');
ylabel('Succinate [mmolh^{-1}gDW^{-1}]');



%subplot(1,1,1)
%axes('position',[0.5 0.5 0.4 0.4]);


fig2=figure(2);
colore=0;


for i=1:step:num_files
    disp(strcat(num2str(i),' out of ',num2str(num_files)));
    current_file = load(strcat('solution',num2str(i),'.mat'));
    chromosome = current_file.chromosome;
    %matrice=[matrice;chromosome];
    
    x = -chromosome(:,num_col-2);
    y = -chromosome(:,num_col-3);


    A = [x y];
    A_sort = sortrows(A,1);
    x = A_sort(:,1);
    y = A_sort(:,2);
    %colore = (step*1/num_files + colore)/colore;

    plot(x,y,'*:','Color',[0 0 colore])
    grid on

    %%%((i-1)*dim(1)-1,:)= current_file.chromosome;
    %%%plot(-chromosome(:,primo),-chromosome(:,secondo),'.');
    hold on
end

grid on

hold off

%axes('position',[0 0 1 1]);
title('');
xlabel('Biomass [h^{-1}]');
ylabel('Succinate [mmolh^{-1}gDW^{-1}]');

[h_m h_i]=inset(fig1,fig2);