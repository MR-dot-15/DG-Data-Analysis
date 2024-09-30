function [dat_extended_list, grp_param, stat] = ...
                test_anova(source_cell, exclude_index, emot, offer)
% TEST_ANOVA(param) for different ANOVA tests
%
% *Input-*
% source_cell: output of behavDatSumary
%
% *Output-*
% p = probab vals returned by ANOVAN

% WITHOUT NORMALIZATION

n_subj = length(source_cell);
subj_list = 1:n_subj;
subj_list = setdiff(subj_list, exclude_index);

if emot == 0
    emotions = 1:3; 
    grp_emot = repmat(emotions, 1, length(subj_list));
    grp_param = grp_emot;
else
    offers = 2:4;
    grp_offer = repmat(offers, 1, length(subj_list));
    grp_param = grp_offer;
end
dat_extended_list = zeros(1, 3*length(subj_list));

% TWO-WAY ANOVA for unbalanced---------------------------------------------

% offers = 2:4;
% grp_subj = zeros(1, 9*n_subj);
% grp_emot = repmat([1 1 1 2 2 2 3 3 3], 1, n_subj);
% grp_off = repmat([2 3 4], 1, 3*n_subj);

% k = 1;
% for subj = 1:n_subj
%     grp_subj(subj:subj+9-1) = subj;
%     for emot = emotions
%         for off = offers
%             key = 10*emot + off;
%             temp = dat{subj}{key};
%             dat_extended_list(k) = sum(temp == 1)/length(temp);
%             k = k + 1;
%         end
%     end
% end
% 
% p = anovan(dat_extended_list', {grp_emot grp_off}, ...
%     'model', "interaction", 'varnames', {'grp_emot', ...
%     'grp_off'});

% ONE-WAY ANOVA for different offers---------------------------------------

j = 1;

if offer == 0
    for subj = subj_list
        for offer = offers
            key = 10*emot + offer;
            temp = source_cell{subj}{key};
            dat_extended_list(j) = sum(temp == 1)/length(temp);
            j = j + 1;
        end
    end
    
    % Kruskal Wallis
    fprintf("1w-KW for emot = %d\n", emot);
    
    [p1, ~, stat] = kruskalwallis(dat_extended_list, grp_offer, 'off');
    fprintf("p = %.3f for len(group) = %d\n",...
        p1, length(dat_extended_list));
    
    % Wilcoxon
    fprintf("Wilcoxon for emot = %d\n", emot);
    
else
    for subj = subj_list
        for emot = emotions
            key = 10*emot + offer;
            temp = source_cell{subj}{key};
            dat_extended_list(j) = sum(temp == 1)/length(temp);
            j = j + 1;
        end
    end
    
    % ANOVA
    % fprintf("1w-ANOVA for offer = %d\n", offer);
    % 
    % [p1, ~, stat] = anova1(dat_extended_list, grp_emot, 'off');
    % fprintf("p = %.3f for len(group) = %d\n",...
    %     p1, length(dat_extended_list));
    
    % Kruskal Wallis
    fprintf("1w-KW for offer = %d\n", offer);
    
    [p1, ~, stat] = kruskalwallis(dat_extended_list, grp_emot, 'off');
    fprintf("p = %.3f for len(group) = %d\n",...
        p1, length(dat_extended_list));
    
    % Wilcoxon
    fprintf("Wilcoxon for offer = %d\n", offer);
end