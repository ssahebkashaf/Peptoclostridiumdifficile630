%% before running, please load geni and nomi from gene_reaction_mappings.xls in the folder model

empty_indices=[];
for i=1:numel(nomi)
    indice=find(strcmp(fbamodel.rxns,nomi{i}));
    
    if isempty(indice)
        empty_indices(end+1) = i; %finds reactions for which there is annotation but they are not in the model (bug)
    else
        fbamodel.grRules{indice}=geni{i};
    end
end

% %for i = 1:numel(empty_indices)
%     fbamodel = addReaction(fbamodel,'Sec_13pd','13pd[e] -->');  %add the missing reaction in the model (849th reaction in the aray of gene annotations, TSV file of the model
%     fbamodel.grRules{end}=geni{849};
% %end

C=[];
for i=1:numel(fbamodel.grRules)
    str = fbamodel.grRules{i};
    to_replace = '((.[^(]*))'
    ciao = regexp(str, to_replace, 'split');
    if ~isempty(str)
        final_genes = textscan(str,'%s','delimiter',{'(',')',' and ',' or '});
        ciao = final_genes{1}(find(~strcmp(final_genes{1},'')));
        C = [C; ciao];
    end
%     if find(strcmp(    C,'r'));
%         pause
%     end
end

C = unique(C);
fbamodel.genes = C;