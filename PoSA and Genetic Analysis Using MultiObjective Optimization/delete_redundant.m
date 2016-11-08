function chromosome = delete_redundant(chromosome,fbamodel)

[M,N]=size(chromosome);
V=fbamodel.nbin;

for i=1:M
    y(i)=sum(chromosome(i,1:V));
end

for i=1:M
    disp(i)
    % analizzo ogni elemento della popolazione
    ind_gene=find(chromosome(i,1:V)==1);
    j=1; %contabiltità dei geni
    while j<=length(ind_gene)
        yt=chromosome(i,1:V);
        yt(1,ind_gene(j))=0;
        % yt vettore colonna
        fbamodel.present = ~(fbamodel.G' * yt');
        ac_old=-chromosome(i,V+1);
        biomass_old=-chromosome(i,V+2);
        [v1, fmax1] = flux_balance(fbamodel,true);
        biomass_new=fbamodel.f' * v1;
        ac_new=fmax1;
        % verifico che le funz obiettivo non cambiano
        if abs(ac_new-ac_old)<=10^-10
            %la funz onj non cambia
            % assegno il nuovo valore a chromosome
            chromosome(i,1:V)=yt;
            %sum(chromosome(i,1:V))
        end%if
        j=j+1;
    end%while
end