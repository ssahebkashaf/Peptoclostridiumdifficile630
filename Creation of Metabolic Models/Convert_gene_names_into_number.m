A=fbamodel.grRules;

for i=1:1093
    tmp=A{i};
    for j=1:806
        tmp=strrep(tmp, char(geni{j}), num2str(fbamodel.genes(j)));
    end
    A{i}=tmp;
end

    