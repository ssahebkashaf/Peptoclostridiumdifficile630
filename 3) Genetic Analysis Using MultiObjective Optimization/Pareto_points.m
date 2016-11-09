function matrice = append_solutions(num_files)

current_file = load(strcat('solution1.mat'));
[dim1,dim2] = size(current_file.chromosome);
matrice = zeros(dim1*num_files,dim2);


for i=1:num_files
    current_file = load(strcat('solution',num2str(i),'.mat'));
    chromosome = current_file.chromosome;
    matrice((i-1)*dim1+1:i*dim1,1:dim2) = chromosome;
    i
end

M = 2; %number of objectives
V = dim2-2-M; %length of input variables without crowding distance, rank, objective functions

matrice = non_domination_sort_mod(matrice, M, V);