
         for i=1:length(reaction_expression)      %for each string of reaction_expression, we replace the substring 'bXXXX' with the number representing its gene expression
             reaction_expression{i}=strrep(reaction_expression{i},  char('(/'),char('('));
         end
        
         for i=1:length(reaction_expression)      %for each string of reaction_expression, we replace the substring 'bXXXX' with the number representing its gene expression
             reaction_expression{i}=strrep(reaction_expression{i},  char('/,'),char(','));
         end
         
         for i=1:length(reaction_expression)      %for each string of reaction_expression, we replace the substring 'bXXXX' with the number representing its gene expression
             reaction_expression{i}=strrep(reaction_expression{i},  char(',/'),char(','));
         end
         
          for i=1:length(reaction_expression)      %for each string of reaction_expression, we replace the substring 'bXXXX' with the number representing its gene expression
             reaction_expression{i}=strrep(reaction_expression{i},  char('/)'),char(')'));
          end
          for i=1:length(reaction_expression)      %for each string of reaction_expression, we replace the substring 'bXXXX' with the number representing its gene expression
             reaction_expression{i}=strrep(reaction_expression{i},  char('/'),char(''));
          end