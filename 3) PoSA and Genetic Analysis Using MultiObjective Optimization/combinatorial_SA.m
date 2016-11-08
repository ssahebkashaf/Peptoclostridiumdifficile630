function [mu1,sigma1]= combinatorial_SA(fbamodel,n_val,K,KK,n_pathway)
%%
% fbamodel: FBA model to analyze
% n_val : number of evaluate
% k: number of step (es k=10)


%% combinatorial sensitivity analysis performs a sensitivity analysis for
% large scale FBA model.
% parameters of the models are the subSystem of the model.

%% fbamodel contains
% 1) subSystem 
% 2) subSystem_cluster 
% 3) SS_reaction 
% 4) SS_genes 

% 0 -> geneset turned on
% 1 -> geneset turned off

%%*************+
%% initialization
mu1=zeros(fbamodel.nSS,1);
sigma1=zeros(fbamodel.nSS,1);

%%

EF1_SS=zeros(K,n_val,fbamodel.nSS);

jj=n_pathway+1;
for j=jj:fbamodel.nSS % for each parameter
    EF1=zeros(K*KK,n_val);
    
    ind_geneset=find(fbamodel.SS_genes(j,:)==1);
    
    if ~isempty(ind_geneset)
        M=length(ind_geneset);
        k=floor(M/K);
        if k==0
            k=1;
        end
        flag=true;
        Kk=1;
        ii=1;
        
        while flag==true
            for i=1:K % livelli
                clc
                            
                YT=zeros(1,fbamodel.nbin);
                fbamodel.present=ones(fbamodel.nrxn,1);
                H1=randi([1 round(fbamodel.nbin*(0.1))],1);
                H2=randi([1 fbamodel.nbin],H1,1);
                YT(H2)=not(YT(H2));
                fbamodel.present = ~(fbamodel.G' * YT');
                [v1, fmax1] = flux_balance(fbamodel,true);           
                v1_ref=v1;
                disp('*************************************************')
                disp(['Subsystem number: ' num2str(j)])
                disp(['LEVEL: ' num2str(i)])
                disp(['repeated ' num2str(Kk) ' volte'])
                geni_spenti=zeros(n_val,1);
                for n=1:n_val
                    %disp(['evaluation: ' num2str(n)])
                    yt=zeros(1,M);
                    h1=i*k;    % livello i
                    if i==K
                        h1=M;
                    end
                    h2=randi([1 h1], 1);
                    h3=randi([1 M], h2, 1);
                    yt(h3)=not(yt(h3));
                    geni_spenti(n)=sum(yt);
                    YT(1,ind_geneset)=yt;                   
                         
                    fbamodel.present = ~(fbamodel.G' * YT');
                    [v1, fmax1] = flux_balance(fbamodel,true);
                    EF1(ii,n)=distance(v1,v1_ref);                    
                end
                media_geni_spenti=mean(geni_spenti);
                Delta=media_geni_spenti/M;
                EF1_SS(ii,:,j)=EF1(ii,:)./Delta;
                ii=ii+1;
            end% for
            Kk=Kk+1;
            solution=['EE_' num2str(j)]; 
            var=EF1_SS(:,:,j);
            save(solution, 'var')
            if Kk>KK
                flag=false;
            end
        end %while   
    end % if
    mu1(j,1)=mean(mean(EF1_SS(:,:,j)));
    sigma1(j,1)=mean(std(EF1_SS(:,:,j)));
end%for