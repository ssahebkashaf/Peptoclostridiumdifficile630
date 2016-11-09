num_populations = 0;
while (exist(['./solution' num2str(num_populations + 1) '.mat'], 'file') == 2)  % check the existance and counts how many solution files are in the folder
    num_populations = num_populations +1;
end 
disp([num2str(num_populations) ' populations detected.'])


primo_file=load('solution1.mat');
dim = size(primo_file.chromosome);
num_cols = dim(2);
num_rows = dim(1);
primo_obj = num_cols-3;
secondo_obj = num_cols-2;

disp([num2str(num_rows) ' individuals in each population.'])


%matrice=[];
matrice=zeros(num_rows*num_populations, num_cols);



for i=1:num_populations
    current_file = load(strcat('solution',num2str(i),'.mat'));
    chromosome = current_file.chromosome;
    matrice((i-1)*num_rows+1:i*num_rows,:) =  chromosome;
    iterations_left=num_populations-i
end

    %matrice = chromosome(:,primo:secondo);

array_aux = zeros(num_populations * num_rows,1);
array_aux_2 = zeros(num_populations * num_rows,1);
for i = 1:num_populations 
    array_aux((i-1)*num_rows+1 : i*num_rows) = i;
end
for i = 1:num_rows 
    array_aux_2(i:num_rows:end) = i;
end
matrice_tracked = [matrice array_aux array_aux_2];  % add last columns to tell (1) from which population the row comes from, and (2) which position it occupied in the population. the number of the comlumns for the first and second objectives are the same.
clear matrice

matrice_sorted = sortrows(matrice_tracked,primo_obj); %sort matrice according to the first objective


obj_sorted = matrice_sorted(:,[end-5 end-4 end-1 end]);

obj_non_redund = zeros(size(obj_sorted,1),size(obj_sorted,2));   %initialised at its maximum possible size, but later we eill eliminate all the final rows that are left with only zeros, and therefore have not been used to save the non redundant solutions
j = 1;

for i=2:size(obj_non_redund,1)
    if abs(obj_sorted(i,1)-obj_sorted(i-1,1))~=0 || abs(obj_sorted(i,2)-obj_sorted(i-1,2))~=0
        obj_non_redund(j,:)=obj_sorted(i,:);
        j = j+1;
    end
    iterations_left=size(obj_non_redund,1)-i
end
       
obj_non_redund = obj_non_redund(sum(abs(obj_non_redund'))~=0,:); % selects only those rows such that the elements are not all zeros, i.e. the sum of the absolute value of the row is nonzero. THis eliminates rows where all the elements are 0, and therefore have not been used from the initialized obj_not_refund matrix



feasible=[obj_non_redund(:,1:2),zeros(size(obj_non_redund,1),1),obj_non_redund(:,3:4)];
non_dominated_app = non_domination(feasible,2,0);

indici_non_dom = find(non_dominated_app(:,3)==1);
non_dominated = - non_dominated_app(indici_non_dom,:);
indici_others = find(non_dominated_app(:,3)==0);
others = - non_dominated_app(indici_others,:);

     %treshold = 1;
%%     %non_dominated = non_dominated(find((non_dominated(:,2) > treshold)==1),:);    %selects only certain rows of the non_dominated and others
     %others = others(find((others(:,2) > treshold)==1),:);

x = non_dominated(:,1);
y = non_dominated(:,2);
 
fig1=figure(1);
plot(-x,y, '*:', 'Color',[0 0 0], 'markers',9)
axis([110,200,0,3.5])

grid on
title('');
%xlabel('1,2-Propanediol [mmolh^{-1}gDW^{-1}]');
%ylabel('Biomass [h^{-1}]');

hold on

%subplot(1,1,1)
%axes('position',[0.5 0.5 0.4 0.4]);
%%
fig2=figure(2);
num_others = size(others,1);

%x = others(:,1);
%y = others(:,2);
%plot(x,y,'*','Color',[0 0 0])

colore=240;     % [200 200 200] is grey, almost white ([255 255 255]). Instead, [0 0 0] is black
for i=1:num_others
    x = others(i,1);
    y = others(i,2);
    colore = 10 + 245*(1 - (-others(i,4))/num_populations);  % others(i,4) contains the number (with the minus sign, though) of the population to which others(i,:) belongs. The +10 avoidt the points are too black and allows for more variability of colours
    plot(-1*x,y,'*','Color',[colore/255 colore/255 colore/255])
    hold on
end

save non_dominated.mat non_dominated;
save others.mat others;

non_dominated(:,4:5)=-non_dominated(:,4:5);
others(:,4:5)=-others(:,4:5);
%xlswrite('others.xlsx', others);
%xlswrite('non_dominated.xlsx', non_dominated);


x = non_dominated(:,1);
y = non_dominated(:,2);
plot(-1*x,y,'*:', 'Color',[0 0 0])

grid on
hold off

%axes('position',[0 0 1 1]);
title('');
%xlabel('1,2-Propanediol [mmolh^{-1}gDW^{-1}]');
%ylabel('Biomass [h^{-1}]');

[h_m h_i]=inset(fig2,fig1);


%exportfig(gcf, 'figura.pdf', 'color', 'cmyk', 'Width', '15', 'Height', '15', 'FontMode', 'scaled', 'FontSize', '1' );
