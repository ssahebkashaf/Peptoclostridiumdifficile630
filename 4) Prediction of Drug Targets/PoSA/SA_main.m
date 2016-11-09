[mu,sigma]= combinatorial_SA(fbamodel,100,4,10,40);
save sol_SA
for i=1:length(mu)
    name=(['file' num2str(i)]);
    name_file=strcat(name,'.txt');
    fid = fopen(name_file,'wt');
    fprintf(fid,'%d\t', mu(i));
    fprintf(fid,'%d', sigma(i));
    fclose(fid);
end

mpdc10 = distinguishable_colors(20)
label=fbamodel.subSystem_cluster;
figure;
A=sigma;
B=mu;
L=label
lg=legend(L)
th = findobj( lg, 'type', 'text' ) 
set( lg, 'fontsize', 9.5 ) 
gscatter(A,B,L,mpdc10,'ox+*sdv^<>ph.' , 10)

%set(findobj('type','axes'),'fontsize',10)

xlabel('\sigma', 'FontSize', 30)
ylabel('\mu', 'FontSize', 30)


