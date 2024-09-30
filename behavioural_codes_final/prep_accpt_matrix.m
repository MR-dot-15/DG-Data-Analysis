function [cond_acceptance, sub_cond_acceptance, not_to_skip] = ...
    prep_accpt_matrix(source_cell, to_skip)

% source_cell: output of behavDatSummary
% to_skip: index of subjects whom you don't wanna include

% init params
not_to_skip = setdiff(1:length(source_cell), to_skip);
n_subj = numel(not_to_skip);
disp(n_subj);

% mean acceptance across all sub
cond_acceptance = zeros(3,3);
% mean acceptance for individual sub
sub_cond_acceptance = cell(3, 3);

for offer = 2:4
    for emot = 1:3
        cell_ind = 10*emot+offer;
        accpt_mean = 0;

        sub_accpt = 1:n_subj;
        iter_subj = 1;
        for subj = not_to_skip
            temp_dat = source_cell{subj}{cell_ind};

            % sub specific acceptance (0-1)
            acceptance = sum(temp_dat == 1)/length(temp_dat);
            sub_accpt(iter_subj) = acceptance;
            iter_subj = iter_subj + 1;

            % cumulatively adding sub spec accpt
            accpt_mean = accpt_mean + acceptance;
        end

        % save sub specific accpt values
        sub_cond_acceptance{emot, offer - 1} = sub_accpt;
        
        % all-over mean
        accpt_mean = accpt_mean/n_subj;
        cond_acceptance(emot, offer-1) = accpt_mean;
    end
end