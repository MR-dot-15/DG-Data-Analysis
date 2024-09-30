function rm = repeat_anova(cond_mat, display_model)
% cond_mat: n_sub x n_meas
% display_model: shows the output of rmfit if true (default false)

if nargin == 1
    display_model = 0;
end

[~, n_meas] = size(cond_mat);

% table
if n_meas == 2
    tab = array2table(cond_mat,...
        "VariableNames", ["c1", "c3"]);
else
    tab = array2table(cond_mat,...
        "VariableNames", ["c1", "c2", "c3"]);
end

% create the within-subjects design
withinDesign = table((1:n_meas)','VariableNames',{'cond'});
withinDesign.Condition = categorical(withinDesign.cond);

% fit the repeat-measure model
rm = fitrm(tab, 'c1-c3 ~ 1', 'WithinDesign', withinDesign);
if display_model
    display(rm);
end

% anova
ranovtab = ranova(rm);
disp("display ranova table ---------------------------------------------");
display(ranovtab);

% mulcompare
multcomptab = multcompare(rm, 'cond');
disp("display multcomp table -------------------------------------------");
display(multcomptab);