function [mu, se] = rm_anova_mean_se(rm_tbl, within_des_param)
% e.g.
% rm_tbl = tbl
% within_des_param = within_des.emot
within_param_val = unique(within_des_param);
mu = 1:numel(within_param_val);
se = 1:numel(within_param_val);

rm_tbl = table2array(rm_tbl(:,2:end));

for idx = 1:numel(within_param_val)
    slice = rm_tbl(:, within_des_param...
                       == within_param_val(idx));
    mu(idx) = mean(slice, 'all');
    se(idx) = std(slice, [], 'all')/...
                  (numel(slice))^.5;
end