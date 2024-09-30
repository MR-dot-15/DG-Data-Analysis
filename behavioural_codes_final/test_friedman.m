function sub_accpt_mat = test_friedman(alldat, sublist)

% sub-specific acceptance cell
disp('n_subj:');
[~, sub_accpt, n_subj] = prep_accpt_matrix(behavDataSummary...
(alldat, sublist, 1, 360, 7));

% prep sub accpt matrix
sub_accpt_mat = zeros(n_subj*3, 3);

for off = 1:3
    temp = [];
    for emot = 1:3
        temp = [temp sub_accpt{emot, off}];
    end
    sub_accpt_mat(:, off) = temp';
end

% just fry the man
[p, ~, stats] = friedman(sub_accpt_mat, n_subj);
fprintf("p-val: %.4f\nmultcomp table---\n", p);
disp(multcompare(stats));