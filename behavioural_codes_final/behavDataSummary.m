function [data, to_skip, block_label] = ...
    behavDataSummary(behavDatCell, llim, ulim, col)
% BEHAVDATASUMMARY(behavDatCell, col) spits out a cell of cells containing 
% the target variable sliced by emotions (1,2,3) X offers (2,3,4)
%
% Inputs:
% subject_names: txt file containing the subject names
% behavDatCell: cell indexed by subject id containing corresp behav data
% llim and ulim: llim-th trial to ulim-th trial--both included
% col: target variable (e.g. 4 for rt, 7 for resp, etc)
% accpt_val: if only acceptance/rejection related values are needed (0:
% both, 1: only accpt, -1: only rejec)
%
% Output:
% a cell of cells
% data{subj_id}{condition} = target var slice for that condition
% block_label{subj_id}{cond} = corresponding block labels
%
% Note that, condition = 10*emotion + offer
%
% Example usage:
% resp_dat = behavDataSummary(subjwisedat, 7)

% some unnecessary bit
% if ulim > 18*5
%     to_skip = horzcat(to_skip, sublist_idx(...
%         contains(sublist, ["mis" "vis"])));
% end
% if ulim > 18*10
%     to_skip = horzcat(to_skip, sublist_idx(...
%         contains(sublist, ["gau" "swa"])));
% end
% if ulim > 18*15
%     to_skip = horzcat(to_skip, sublist_idx(...
%         contains(sublist, ["kum" "biv"])));
% end

% subjects
subj_n = length(behavDatCell);
to_skip = [];
% if you want acceptance or rejection-spec slice
accpt_val = 0;

% preparing data cells
data = cell(1, subj_n-length(to_skip));
block_label = cell(1, subj_n-length(to_skip));
cond_slice = cell(1, 3);
cond_blk_slice = cell(1, 3);

for subj = 1:subj_n-length(to_skip)
    % display names
    % name = string(subjectNames(subj));
    % fprintf("%s ...\n", name);

    % access the behav data in the cell
    temp = behavDatCell{subj}; 
    %ulim = size(temp, 1);
    temp = temp(llim:ulim, :); 
    emot_col = temp(:,1);
    off_col = findBaseOff(temp(:,2));
    blk_idx = reshape(repmat(1:20, 18, 1), [], 1);

    % special slicing for aviyank: not needed anymore
    % load('aviyank_off.mat', 'aviyank_off');
    % aviyank_off = aviyank_off(llim:ulim);
    % off2 = (aviyank_off == 2); 
    % off3 = (aviyank_off == 3); 
    % off4 = (aviyank_off == 4);
    %k = 1;
    for emot = 1:3
        for offr = 2:4
            % "cell_index" defined here
            cell_index = 10*emot + offr; 

            %   change this cell idx later
            %cell_index = k; k = k+1;
            
            if accpt_val == 0
                slice = and(emot_col == emot,...
                            off_col == offr);
            else
                slice = emot_col == emot &...
                        off_col == offr &...
                        temp(:, 7) == accpt_val;
            end
            sliced_dat = temp(slice, col);
            sliced_blk_idx = blk_idx(slice);

            % store
            cond_slice{cell_index} = sliced_dat';
            data{subj} = cond_slice;
            cond_blk_slice{cell_index} =...
                            sliced_blk_idx;
            block_label{subj} = cond_blk_slice;
        end
    end
end
