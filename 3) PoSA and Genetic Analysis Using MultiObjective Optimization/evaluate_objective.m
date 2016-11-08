function f = evaluate_objective(x, M, V, fbamodel, geni, reaction_expression)

%% function f = evaluate_objective(x, M, V)
% Function to evaluate the objective functions for the given input vector
% x. x is an array of decision variables and f(1), f(2), etc are the
% objective functions. The algorithm always minimizes the objective
% function hence if you would like to maximize the function then multiply
% the function by negative one. M is the numebr of objective functions and
% V is the number of decision variables. 
%
% This functions is basically written by the user who defines his/her own
% objective function. Make sure that the M and V matches your initial user
% input.

yt=x';      % x' is the transpose of x, that is "child(j,1:V)" in the function genetic_operator, i.e. the array of gene expressions
yt=yt(1:V);     % needed because sometimes, and especially with child_3 in genetic_operator, there is a bug and all the child is passed to this function here, including the final rank, crowding distance and objective functions (instead of passing only the V decision variables)
fbamodel.present=ones(fbamodel.nrxn,1);

%fbamodel.present = ~(fbamodel.G' * yt);   %istruzione vecchio approccio knockout

%NEW CODE

%load('geni.mat');
%load('reaction_expression.mat');

%geni( cellfun(@isempty, reaction_expression) ) = {'1'};  %replaces all the empty names of genes with the name '1'
%reaction_expression( cellfun(@isempty, reaction_expression) ) = {'1'};  %replaces all the empty cells of gene expression (e.g. exchange reactions) with '1', i.e. gene expressed nomally

%reaction_expression = cellfun(@(c) c{1},reaction_expression);
%geni = cellfun(@(c) c{1},geni);

for i=1:length(yt)   %loop over the array of the gene expressions, in order to replace the names of genes in geni_reazioni.mat with their values
    if (strcmp(geni{i},'')==0) 
        %sost_reaction_expression = strrep({reaction_expression},  char(geni{i}), char((yt(i))));  %replaces all the names of the i-th gene (that are present in reaction_expression) with its numerical value
        %sost_reaction_expression = regexprep(num2str(reaction_expression), num2str(geni{i}), num2str(yt(i))) ;
        %reaction_expression(strcmp(geni(i), reaction_expression))={num2str(yt(i))};
        
        %substrmatch = @(x,y) ~cellfun(@isempty,strrep(y,  x, num2str(yt(i))));
        %substrmatch(geni{i},reaction_expression);
        
        %cellfun(strrep,reaction_expression,  geni{i}, num2str(yt(i)));
        
         matches = strfind(reaction_expression,geni{i});     %this and the following instruction find the locations of the gene 'bXXXX' in the array reaction_expression
         posizioni_gene = find(~cellfun('isempty', matches));
         for j=1:length(posizioni_gene)      %for each string of reaction_expression, we replace the substring 'bXXXX' with the number representing its gene expression
             %tmp=yt(i)*CAI(i);
             reaction_expression{posizioni_gene(j)}=strrep(reaction_expression{posizioni_gene(j)}, ['/' geni{i} '/'] , num2str(yt(i)));
             reaction_expression{posizioni_gene(j)}=strrep(reaction_expression{posizioni_gene(j)}, ['/' geni{i} ' /'] , num2str(yt(i)));
             reaction_expression{posizioni_gene(j)}=strrep(reaction_expression{posizioni_gene(j)}, ['/ ' geni{i} '/'] , num2str(yt(i)));
         end
        
    end
end
reaction_expression( cellfun(@isempty, reaction_expression) ) = {'1'};  %replaces all the empty cells of gene expression (e.g. exchange reactions) with 1, i.e. gene expressed nomally

num_reaction_expression=zeros(1,length(reaction_expression));

for i=1:length(reaction_expression)
    i
    reaction_expression{i}
    num_reaction_expression(i)=eval(reaction_expression{i});   %evaluates the cells like they are numerical expressions (so as to compute min and max of gene expressions)
end

for i=1:length(num_reaction_expression)   %loop over the array of the geneset_exressions gene expressions
    
%    aux=find(fbamodel.G(i,:)==1);   %reactions that depend on the i-th geneset
    
%    for j=1:length(aux)
        %fbamodel.vmin(aux(j)) = fbamodel.vmin(aux(j))*yt(i);
        %fbamodel.vmax(aux(j)) = fbamodel.vmax(aux(j))*yt(i);
        if num_reaction_expression(i)>=1
            fbamodel.vmin(i) = fbamodel.vmin(i)*(1+log(num_reaction_expression(i)));
            fbamodel.vmax(i) = fbamodel.vmax(i)*(1+log(num_reaction_expression(i)));
        else
            fbamodel.vmin(i) = fbamodel.vmin(i)/(1+abs(log(num_reaction_expression(i))));
            fbamodel.vmax(i) = fbamodel.vmax(i)/(1+abs(log(num_reaction_expression(i))));
        end
%   end
end
%END NEW CODE

[v1, fmax] = flux_balance(fbamodel,false);

% objective functions number M is 2
f(1)=-fmax; % acetate
f(2)=-fbamodel.f' * v1; % Biomass

%% Check for error
if length(f) ~= M
    error('The number of decision variables does not match you previous input. Kindly check your objective function');
end