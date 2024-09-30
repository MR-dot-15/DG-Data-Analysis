function allBehavDat = collateAll(markers_list)
fnames = [];
for marker = markers_list
    files = dir(marker);
    fnames_cell = {files.name};
    fnames = horzcat(fnames, string(fnames_cell));
end

n = length(fnames);

fileID = fopen("subject.txt", "w");
allBehavDat = cell(1,n);

for i = 1:n
    % save the subject code to "subject.txt"
    % and save the subject-specific 
    % behav data to the "allBehavDat" cell
    to_print = char(fnames(i));
    behavDat = readmatrix(to_print);

    to_print = to_print(1:5) + "\n";

    fprintf(fileID, to_print);

    allBehavDat{i} = behavDat;
end
fclose(fileID);