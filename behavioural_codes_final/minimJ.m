function [p_thresh, idx] = minimJ(label, score)
% returns threshold p-val by minimizing Youden's J

% def threshold vals
thresh = 0:.05:1;
% store J-vals for each iter
J_vals = zeros(size(thresh));

for iter_idx = 1:numel(thresh)
    % predict labels
    label_hat = zeros(size(label));
    label_hat(score > thresh(iter_idx)) = 1;

    % calc confusion metrics
    tp = sum(label == 1 & label_hat == 1);
    fp = sum(label == 0 & label_hat == 1);
    fn = sum(label == 1 & label_hat == 0);
    tn = sum(label == 0 & label_hat == 0);

    % calculate and store J statistic
    J_vals(iter_idx) = tp/(tp + fn) + tn/(tn + fp) - 1;
    %J_vals(iter_idx) = tp / (tp + mean(fp + fn));
end

[~, idx] = max(J_vals); p_thresh = thresh(idx);