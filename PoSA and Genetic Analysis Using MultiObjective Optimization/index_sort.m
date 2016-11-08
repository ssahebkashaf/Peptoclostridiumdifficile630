load sol_SA
%mu,sigma
N=length(mu);

%
ind=find(mu==max(mu));
punto_max=[mu(ind) sigma(ind)];

% d=size(N,1);
for i=1:N
    d(i,1)= sqrt( (punto_max(1)-mu(i)).^2 + (punto_max(2)-sigma(i)).^2 );
end
load geobacter_sulf

fid= fopen('plot_gnuplot.txt','wt');
%"file15.txt"  u 1:2       with points title "Tyrosine, Tryptophan, and
%Phenylalanine Metabolism"  lt -1 pt 7 pointsize 1.3,\
[val,ind]=sort(d);
for i=1:fbamodel.nSS
    fprintf(fid,'%s \n',['"file' num2str(ind(i)) '.txt" u 1:2 with points title "'...
        char(fbamodel.subSystem_cluster{ind(i)}(1)) '" lt ' num2str(randi(10,1,1)) ' pt ' num2str(randi(10,1,1)) ',\']);
    fbamodel.subSystem_cluster{ind(i)}(1)
    
end

fclose(fid);