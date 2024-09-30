function [sub_mean_rt, rt_dat, sub_idx] = prep_rt_matrix(source_cell)
% source_cell: output of behavDataSummary
% output --
% sub_mean_rt: subject-wise mean rt across conds
% rt_dat: pulled rt dat across conds
% sub_idx: sub-wise grouping data for rt_dat

% NB: include input param sub indices

n_subj = length(source_cell);
disp(n_subj);

sub_mean_rt = cell(3, 3);
rt_dat = cell(3, 3);
sub_idx = cell(3, 3);

for emot = 1:3
    for offer = 2:4
        cell_ind = 10*emot+offer;

        pulled_for_cond = [];
        sub_mean_cond = [];
        sub_idx_cond = [];

        for subj = 1:n_subj
            temp_dat = source_cell{subj}{cell_ind};
            pulled_for_cond = horzcat(pulled_for_cond,...
                                        temp_dat);
            sub_idx_cond = horzcat(sub_idx_cond,...
                    subj*ones(1, length(temp_dat)));
            sub_mean_cond(end + 1) = mean(temp_dat);
        end

       rt_dat{emot, offer - 1} = pulled_for_cond;
       sub_mean_rt{emot, offer - 1} = sub_mean_cond;
       sub_idx{emot, offer - 1} = sub_idx_cond;
    end
end