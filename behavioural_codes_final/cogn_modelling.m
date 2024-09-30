%% modelling reation time
end_block_idx = []; % 1-5, 1-10, ... etc
R_sq_idx = []; % adj R-sq across all iterations (pulled)

for end_block = 20
    % some sub with incomplete data
    %to_skip = findToSkipSub(end_block, sublist); 
    to_skip = [];
    % remove those sub
    subject_iter = setdiff(1:length(alldat), to_skip); 
    % remaining no of sub
    n_subj = length(subject_iter);
    % as R-sq is pulled, to find the idx we use
    % R-sq(prev_len : end)
    prev_R_sq_len = length(R_sq_idx); 

    % store model fits
    lm_store = {};
    
    % fill end_block indices
    end_block_idx = [end_block_idx ...
        (end_block/5) * ones(1, n_subj)];

    for sub_idx = subject_iter
        dat = alldat{sub_idx};    
        
        % extract the variables: rt, emot, off
        rt = log(dat(1:end_block*18, 4));
        remove_idx = rt < .25 & rt == 2;
        rt = rt(~remove_idx);
    
        emot = dat(~remove_idx,1);
        % 0,1,2
        off = dat(~remove_idx,2); off = off ./ max(off);
        block_id = reshape(...
            repmat(1:end_block, 18, 1),...
            1, [])';
        block_id = (block_id - 1)./19;
        block_id = block_id(~remove_idx);
        accpt = dat(~remove_idx, 7) + 1;
        
        % create the table
        tbl = table(emot, off, block_id, accpt, rt);
        tbl.emot = categorical(tbl.emot);
        %tbl.off = categorical(tbl.off);
        tbl.block_id = categorical(tbl.block_id);
        tbl.accpt = categorical(tbl.accpt);
        
        % fit the model
        % lm = fitlme(tbl, ['rt ~ 1 + emot*off +' ...
        %                   '(emot | block_id) +'...
        %                   '(off | block_id)']);
        lm = fitlme(tbl, ['rt ~ 1 + emot*off']);
        % lm = fitlm(tbl,...
        %    'rt ~ emot + off + emot*off + block_id');
        
        % store R_sq values, models
        R_sq_idx = [R_sq_idx lm.Rsquared.Adjusted];
        lm_store{end+1} = lm;
    end
    
    % find the best R-sq in this iter
    [best_fit, best_fit_idx] = max(R_sq_idx...
                (prev_R_sq_len+1 : end));

    fprintf("best R-sq %.3f for sub %s--------------------\n", best_fit,...
        sublist(subject_iter(best_fit_idx)));
    %disp(lm_store{best_fit_idx})
    disp(lm_store{best_fit_idx}.anova);
end

figure;
boxchart(end_block_idx, R_sq_idx);

%% try to classify based on emot or off effect
% 1:15 blocks used in the model
% from the anova table, extract p-val

emot_p = zeros(1, n_subj);
off_p = zeros(1, n_subj);
emot_off_p = zeros(1, n_subj);

for sub_idx = 1:40
    anova_tab = lm_store{sub_idx}.anova.pValue;
    % emot_p(sub_idx) = anova_tab{"emot", "pValue"};
    % off_p(sub_idx) = anova_tab{"off", "pValue"};
    % emot_off_p(sub_idx) = anova_tab{"emot:off", "pValue"};
    emot_p(sub_idx) = anova_tab(2);
    off_p(sub_idx) = anova_tab(3);
    emot_off_p(sub_idx) = anova_tab(4);
end

fprintf("with effect of emot:\n___________________________\n");
for sub_idx = 1:40
    if emot_p(sub_idx) < .1
        fprintf("sub: %s\tp: %.3f\tR: %.4f\n",...
        sublist(sub_idx), emot_p(sub_idx), R_sq_idx(sub_idx));
    end
end

fprintf("\nwith effect of offer:\n___________________________\n");
for sub_idx = 1:40
    if off_p(sub_idx) < .1
        fprintf("sub: %s\tp: %.3f\tR: %.4f\n",...
        sublist(sub_idx), off_p(sub_idx), R_sq_idx(sub_idx));
    end
end

fprintf("\nwith effect of offer-emot interaction:" + ...
    "\n___________________________\n");
for sub_idx = 1:40
    if emot_off_p(sub_idx) < .1
        fprintf("sub: %s\tp: %.3f\tR: %.4f\n",...
        sublist(sub_idx),...
        emot_off_p(sub_idx),...
        R_sq_idx(sub_idx));
    end
end

%% modelling acceptance, rejection 
% resp ~ off + emot + rt + subID +
%        off*sub + emot*sub + emot*off*sub

% init params
end_block = 10;
accpt_val = [];
subID = [];
emot = [];
off = [];
rt = [];

% collate data
to_skip = findToSkipSub(end_block, sublist); 
subject_iter = setdiff(1:length(alldat), to_skip); 

for sub_idx = subject_iter
    dat = alldat{sub_idx};
    slice_idx = dat(:, 4) ~= 2 &...
                dat(:, 4) > .25;

    emot = [emot dat(slice_idx,1)'];
    off = [off dat(slice_idx,2)']; off = off./max(off);
    subID = [subID...
        sub_idx * ones(1, sum(slice_idx))];
    rt = [rt log(dat(slice_idx,4))'];
    accpt_val = [accpt_val dat(slice_idx,7)'];
end
accpt_val(accpt_val == -1) = 0;

% build the model
tbl = table(emot', off', subID', rt', accpt_val',...
        'VariableNames', ["emot" "off" "subID" "rt" "accpt"]);
tbl.emot = categorical(tbl.emot);
%tbl.off = categorical(tbl.off);
tbl.subID = categorical(tbl.subID);

% fit the model
exp = 'accpt ~ (rt|subID) + (off|subID) + (emot|subID)';
glm = fitglme(tbl, exp, 'Distribution', 'Binomial',...
                    'Link', 'logit');

%% Cross-validation with different models
figure; warning('off');
fprintf("%s %10s %10s %10s\n",...
        "rt R2", "resp R2", "AUC 1", "BIC");

% define models for accpt
exp_array = "accpt ~ 1 + rt_res + emot + off +"+...
                "(1|subID) + (1|bID) +"+...
                "(rt_res-1|subID) + (emot-1|subID) + (off-1|subID) +"+...
                "(rt_res-1|bID) + (emot-1|bID) + (off-1|bID)";
% exp_array = "accpt ~ 1 + rt_res + emot + off +"+...
%             "(1|subID) + (1|bID) +"+...
%             "(rt_res-1|subID) + (emot-1|subID) + (off-1|subID) +"+...
%             "(rt_res-1|bID) + (emot-1|bID) + (off-1|bID)";

k = 1;
for exp = exp_array
    % init storing variables and parameters
    cv_iter = 10;
    run_iter = 3;
    n_elements = run_iter*numel(exp_array);
    rt_mdl = cell(n_elements,1);
    resp_mdl = cell(n_elements,1);
    train_dataset = cell(n_elements,1);
    test_dataset = cell(n_elements, 1);
    auc = 1:numel(n_elements);
    
    % % collate the entire data into one matrix
    % collated_dat = zeros(1,7);
    % subID = [];
    % blockID = [];
    % for sub_idx = 1:numel(alldat)
    %     dat = alldat{sub_idx};
    %     temp_blockID = reshape(repmat(1:20, 18, 1), 1, []);
    % 
    %     % RT filtering
    %     slice_idx = dat(:,4) ~= 2 & dat(:,4) > .25;
    %     collated_dat = cat(1, collated_dat, dat(slice_idx,:));
    %     subID = [subID...
    %         sub_idx * ones(1, sum(slice_idx))];
    %     blockID = [blockID temp_blockID(slice_idx)];
    % end
    % collated_dat = collated_dat(2:end, :);
    % 
    % % make the table
    % tbl = table(collated_dat(:,1),...
    %             collated_dat(:,2)./max(collated_dat(:,2)),...
    %             subID', blockID', log(collated_dat(:,4)),...
    %             zeros(size(subID')), (collated_dat(:,7)+1)./2,...
    %             'VariableNames', ["emot" "off" "subID" "bID"...
    %                               "rt" "rt_res" "accpt"]);
    % tbl.emot = categorical(tbl.emot);
    % tbl.subID = categorical(tbl.subID);
    
    % make table
    tbl = makeModelTable(alldat);
    
    % CV partition
    n = size(tbl, 1);
    c = cvpartition(n, "KFold", cv_iter);
        
    % exp_array = "accpt ~ rt + emot + off +" + ...
    %            " (rt|subID) +" +...
    %            " (emot|subID) +" +...
    %            " (off|subID)" ;
    
    for iter_idx = 1:run_iter
        % partition train-test
        tblTrain = tbl(training(c, iter_idx), :);
        train_dataset{iter_idx} = tblTrain;
        tblTest = tbl(test(c, iter_idx), :);
        test_dataset{iter_idx} = tblTest;
    
        % model rt from offer, emot ---
        lm = fitlme(tblTrain, ['rt ~ 1 +'...
                                'off + emot +'...
                                '(1|subID) + (1|bID) +'...
                                '(off-1|bID) + (off-1|subID) +'...
                                '(emot-1|bID) + (emot-1|subID)']);
        rt_mdl{iter_idx} = lm;   
        rt_fitted = lm.fitted;
    
        % model acceptance/rejection ---
        tblTrain.rt_res = tblTrain.rt - rt_fitted;
        glm = fitglme(tblTrain, char(exp), 'Distribution', 'Binomial',...
                            'Link', 'logit', 'StartMethod', 'default');
        resp_mdl{iter_idx} = glm;
        [~, ~, ~, auc_1] = perfcurve(tblTrain.accpt,...
                                     glm.fitted, 1);
        fprintf("%.2f \t %.3f", auc_1, glm.ModelCriterion.BIC); 

        % test, first rt then resp
        rt_fitted = predict(lm, tblTest);
        tblTest.rt_res = tblTest.rt - rt_fitted;
        pred_score = predict(glm, tblTest);
        
        % plot and display the auc value
        [x1, y1, ~, auc_1] = perfcurve(tblTest.accpt, pred_score, 1);
        disp(auc_1);
        auc(iter_idx) = auc_1;
        % fprintf("%.3f %10.3f %10.3f %10.3f\n",...
        %     lm.Rsquared.Adjusted, glm.Rsquared.Adjusted,...
        %     auc_1, glm.ModelCriterion.BIC);
    end
    disp("---------------------");
    % plot for each model
    subplot(2, numel(exp_array), k);
    plot(x1, y1); xlabel("FPR"); ylabel("TPR"); title(k);
    subplot(2, numel(exp_array), k + numel(exp_array)); k = k + 1;
    histogram(pred_score(tblTest.accpt == 1), 0:.05:1,...
        'FaceAlpha', .5); hold on;
    histogram(pred_score(tblTest.accpt == 0), 0:.05:1,...
        'FaceAlpha', .5); hold off;
    xlabel('predicted probability'); 
end

% close all;

%% visualize the coeff
% figure;
% subplot(2, 1, 1);
% temp_mu = lm.Coefficients.Estimate(2:end);
% temp_se = lm.Coefficients.SE(2:end);
% errorbar(1:3, temp_mu, temp_se, 'b-',... 
%          'LineWidth', 2, 'CapSize', 8,...
%          'Marker', 'o',...
%          'MarkerFaceColor', [.8 .8 .8],...
%          'MarkerEdgeColor', 'w');
% xticks(1:3);
% xticklabels(["e_2:Happy" "e_3:Disgusted" "Offer"]);
% xlim([0 4]);
% yline(0, 'k--', 'LineWidth', 1);
% ylabel("Coefficients");
% ylim([-0.04 0.027]);
% 
% subplot(2, 1, 2);
% temp_mu = glm.Coefficients.Estimate(2:end);
% temp_se = glm.Coefficients.SE(2:end);
% errorbar(1:4, temp_mu, temp_se, 'k-',... 
%          'LineWidth', 2, 'CapSize', 8,...
%          'Marker', 'o',...
%          'MarkerFaceColor', [.8 .8 .8],...
%          'MarkerEdgeColor', 'w');
% errorbar(1:4, temp_mu, temp_se, 'k-',... 
%          'LineWidth', 2, 'CapSize', 5,...
%          'Marker', 'o',...
%          'MarkerFaceColor', [.8 .8 .8],...
%          'MarkerEdgeColor', 'w');
% xticks(1:4);
% xticklabels(["e_2:Happy" "e_3:Disgusted" "Offer" "rt_r"]);
% xlim([0 5]);
% ylim([-4.3 5]);
% yline(0, 'k--', 'LineWidth', 1);
% ylabel("Coefficients");

%% periodic oscillation of coefficients with block ID
% find the block-spec random effects
[~, ~, stat] = randomEffects(glm);
e_2_b = double(stat(241:2:280,4));
e_3_b = double(stat(242:2:280,4));
off_b = double(stat(281:300,4));

% plot
figure; subplot(3, 1, 1); title('coeff disgust | bID');
hold on;
xline(10, 'k--'); yline(0, 'k');
plot(1:20, e_2_b, 'ro-', 'MarkerFaceColor', 'r');
hold off;
subplot(3,1,2); title('coeff neutral | bID');
hold on;
xline(10, 'k--'); yline(0, 'k');
plot(1:20, e_3_b, 'bo-', 'MarkerFaceColor', 'b');
hold off;
subplot(3,1,3); title('coeff offer | bID');
hold on;
xline(10, 'k--'); yline(0, 'k');
plot(1:20, off_b, 'ko-', 'MarkerFaceColor', 'k');

% note: offer coefficient being -ve doesn't mean "actual" coeff is -ve
% the fixed effect for offer is quite large (~15)
% this fluctuation is minor

%% rt vs rt_hat
figure;
cv_iter = 1;

for iter_idx = 1:cv_iter
    temp_mdl = lm; %rt_mdl{iter_idx};
    %tblTrain = train_dataset{iter_idx};
    tblTrain = lm.Variables;

    % lvl-1 (original rt and fitted rt - histograms)
    subplot(cv_iter, 3, iter_idx); hold on;
    histogram(tblTrain.rt); 
    histogram(temp_mdl.fitted);
    histogram(tblTrain.rt - temp_mdl.fitted);
    legend(["rt", "rt hat", "rt res"]);
    xlabel("t"); 
    title(temp_mdl.Rsquared.Adjusted);
    hold off;

    % lvl-2 (linear rel-n b/w rt and rt hat)
    subplot(cv_iter, 3, cv_iter + iter_idx); hold on;
    scatter(tblTrain.rt, temp_mdl.fitted, 'filled',...
            'MarkerFaceAlpha', .2);
    %plot(log(.25):log(2), log(.25):log(2), 'k--');
    plot(-1.5:.1:1, -1.5:.1:1, 'k--');
    xlabel('rt'); ylabel('rt-hat'); hold off;

    % lvl-3 (residual)
    subplot(cv_iter, 3, 2*cv_iter + iter_idx); 
    [~, idx] = sort(tblTrain.rt, 'ascend');
    rt_res = tblTrain.rt - temp_mdl.fitted;
    stem(rt_res(idx), 'filled', 'k--.'); 
    ylabel('res');
end

sgtitle("with log");
clear idx rt_res temp_mdl;

%% model figure
figure;
iter = 2;
temp_mdl = rt_mdl{iter}; glm = resp_mdl{iter};
tblTrain = train_dataset{iter};
tblTest = test_dataset{iter};
pred_score = predict(glm, tblTest);
[x1, y1, ~, auc_1] = perfcurve(tblTest.accpt, pred_score, 1);

% plot-1 rt vs rt_hat
subplot(2, 2, 1); hold on;
scatter(tblTrain.rt, temp_mdl.fitted, 'filled',...
        'MarkerFaceAlpha', .2);
%plot(log(.25):log(2), log(.25):log(2), 'k--');
xlim([-1.5 1]); ylim([-1.5 1]); 
plot(-1.5:.1:1, -1.5:.1:1, 'k--', 'LineWidth', 2);
title(sprintf('Showcasing the Reaction Time Fit'), 'FontSize', 18); 
set(gca, 'box', 'off');
ax = gca; ax.XAxis.FontSize = 12; ax.YAxis.FontSize = 12;
xlabel('Observed log(RT)', 'FontSize', 15); 
ylabel('Predicted log(RT)', 'FontSize', 15);

% plot-2 residual
subplot(2, 2, 2); 
[~, idx] = sort(tblTrain.rt, 'ascend');
rt_res = tblTrain.rt - temp_mdl.fitted;
stem(rt_res(idx), 'filled', 'k--.'); 
xlim([0 12200]); ylim([-1.3 1.3]);
xticks(0:4000:12000); 
xticklabels(["0" "4 \cdot 10^3" "8 \cdot 10^3" "12 \cdot 10^3"]);
yline(0, 'w--', 'LineWidth', 2);
title(sprintf('Residual from the fitted model'), 'FontSize', 18); 
ax = gca; ax.XAxis.FontSize = 12; ax.YAxis.FontSize = 12;
xlabel('Observations in the train set', 'FontSize', 15);
ylabel('Residual of log(RT)', 'FontSize', 15);
set(gca, 'box', 'off');

% plot-3 AUC curve
subplot(2, 2, 3); hold on;
plot(x1, y1, 'r-', 'LineWidth', 4); 
plot(0:.1:1, 0:.1:1, 'k--', 'LineWidth', 2);
ylim([0 1.1]); yline(1, 'k:', 'LineWidth', 2);
yticks(0:.2:1);
title(sprintf('ROC for the Response Model'), 'FontSize', 18); 
ax = gca; ax.XAxis.FontSize = 12; ax.YAxis.FontSize = 12;
xlabel("FPR", 'FontSize', 15); ylabel("TPR", 'FontSize', 15); 
set(gca, 'box', 'off');

% plot-4 separation plot
subplot(2, 2, 4);
histogram(pred_score(tblTest.accpt == 1), 0:.04:1,...
    'FaceAlpha', .5, 'EdgeColor','none'); hold on;
histogram(pred_score(tblTest.accpt == 0), 0:.04:1,...
    'FaceAlpha', .5, 'EdgeColor','none'); hold off;
legend(["Accept", "Reject"], 'FontSize', 15);
xlim([0 1]);
title(sprintf('Separability Plot'), 'FontSize', 18); 
set(gca, 'box', 'off'); 
ax = gca; ax.XAxis.FontSize = 12; 
ax.YAxis.FontSize = 12; ax.YScale = 'log';
xlabel('P(Acceptance)', 'FontSize', 15); 
ylabel('log(Count)', 'FontSize', 15);

%% grouping subjects based on regression coefficients : resp mdl
% extract coeff
glm = resp_mdl{1};
[~, ~, stat] = randomEffects(glm);
e_2_s = double(stat(101:2:180,4));
e_3_s = double(stat(102:2:180,4));
off_s = double(stat(181:220,4));

% plot
figure;
scatter3(e_2_s, e_3_s, off_s, 'filled'); 
text(e_2_s+.005, e_3_s+.005, off_s+.005, string(1:40),...
     'FontSize', 8);
xlabel("e_2", 'FontSize', 13); 
ylabel("e_3", 'FontSize', 13); 
zlabel("off", 'FontSize', 13);

%%
% save a rotating animation
% Set the view angle
% view_angle = [150, 30];
% 
% % Number of frames for the animation
% num_frames = 90;
% 
% % Create and save each frame as a PNG image
% for i = 1:num_frames
%     % Update the view angle
%     view_angle = view_angle + [0, 2];
%     view(view_angle);
% 
    % % Capture the frame
    % frame = getframe(gcf);
    % im = frame2im(frame);
    % 
    % % Save the frame as a PNG file
    % filename = sprintf('frame_%03d.png', i);
    % imwrite(im, filename, 'png');
% end
% 
% n_process = 10;
% ddm_chains = cell(1, n_process);
% 
% %v = .5; 
% sd = 1; 
% X = [0];
% k = 1;
% temp_rand = rand(n_process);
% for i = 1:n_process
%     temp = temp_rand > .5;
%     if temp == 1
%         v = .2;
%     else
%         v = -.2;
%     end
%     for t = 1:60
%         X(end+1) = X(end) + v + sd * randn(1);
%         if X(end) >= 5 || X(end) <= -5 || t == 1
%             %gcf(); close; figure;
%             X = [0];
%             continue
%         else
%             figure; plot(X, 'k-', 'LineWidth', 1.5); hold on;
%             xlim([-5 60]);
%             yline([5 -5], 'r-', 'LineWidth', 2);
%             xline(0, 'b--'); xticks(-5:5:60);
%             xticklabels(["" string(num2str((0:5:60)'))']);
%             yticklabels([]); 
%             text(55, 4, "A"); text(55, -4, "R");
%             xlabel("Time (s)"); ylabel('Latent Decision Variable');
%             fontsize(13, 'points'); set(gca(), 'Color', [243 243 243]./256)
% 
%             % Capture the frame
%             frame = getframe(gcf);
%             im = frame2im(frame);
% 
%             % Save the frame as a PNG file
%             filename = sprintf('frame_%03d.png', k); k = k + 1;
%             imwrite(im, filename, 'png');
%         end
%         %input('');
%         close all;
%     end
% end

%% clustering

e_2 = e_2_s ./ (max(e_2_s) - min(e_2_s));
e_3 = e_3_s ./ (max(e_3_s) - min(e_3_s));
off = off_s ./ (max(off_s) - min(off_s));

figure;
scatter3(e_2, e_3, off, 'filled'); text(e_2, e_3, off, string(1:40));


C_ = linkage(pdist([e_2 e_3 off]), 'centroid');
C = cluster(C_, 'criterion', 'distance', 'cutoff', .35); %.33

figure; hold on;
for iter = unique(C)'
    scatter3(e_2(C == iter),...
             e_3(C == iter),...
             off(C == iter), 60, 'filled');
    
    text(e_2(C == iter) + .05,...
         e_3(C == iter) + .05,...
         off(C == iter) + .05,...
         sublist_idx(C' == iter),...
         'FontSize', 8);
end
grid('on');
xlabel("e_2 : Happy", 'FontSize', 13); 
ylabel("e_3 : Disgusted", 'FontSize', 13); 
zlabel("offer", 'FontSize', 13);

%% model clustering figure -- work only on model fitted on the entire data
figure;
grp_cmap = parula(7);

fixed = glm.Coefficients.Estimate(2:4);
[~, ~, stat] = randomEffects(glm);
e_2_s = double(stat(101:2:180,4));
e_3_s = double(stat(102:2:180,4));
off_s = double(stat(181:220,4));
e_2 = e_2_s ./ (max(e_2_s) - min(e_2_s));
e_3 = e_3_s ./ (max(e_3_s) - min(e_3_s));
off = off_s ./ (max(off_s) - min(off_s));

C_ = linkage(pdist([e_2 e_3 off]), 'centroid');
leafOrd = optimalleaforder(C_, pdist([e_2 e_3 off]));
C = cluster(C_, 'criterion', 'distance', 'cutoff', .3);

% dendogram
subplot(2, 3, 1);
h = dendrogram(C_, 40, 'Orientation', 'top',...
    'ColorThreshold', .3, 'Reorder', leafOrd); 
set(h, 'LineWidth', 1.5);
set(gca, 'FontSize', 8);

% now change the dend colours hehe
dend_cmap = zeros(size(h, 1), 3);
mod_dend_cmap = zeros(size(h, 1), 3);
for iter = 1:size(h, 1)
    dend_cmap(iter, :) = h(iter).Color;
end
uniq_dend_cmap = unique(dend_cmap, "rows");

mod_dend_cmap(all(dend_cmap == uniq_dend_cmap(2, :), 2), :) =...
    repmat(grp_cmap(7,:),...
    sum(all(dend_cmap == uniq_dend_cmap(2, :), 2)), 1); % offer

mod_dend_cmap(all(dend_cmap == uniq_dend_cmap(3, :), 2), :) =...
    repmat(grp_cmap(5,:),...
    sum(all(dend_cmap == uniq_dend_cmap(3, :), 2)), 1); % emot - strong

mod_dend_cmap(all(dend_cmap == uniq_dend_cmap(4, :), 2), :) =...
    repmat(grp_cmap(4,:),...
    sum(all(dend_cmap == uniq_dend_cmap(4, :), 2)), 1); % emot - mid

mod_dend_cmap(all(dend_cmap == uniq_dend_cmap(5, :), 2), :) =...
    repmat(grp_cmap(1,:),...
    sum(all(dend_cmap == uniq_dend_cmap(5, :), 2)), 1); % AA, AR

mod_dend_cmap(all(dend_cmap == uniq_dend_cmap(6, :), 2), :) =...
    repmat(grp_cmap(6,:),...
    sum(all(dend_cmap == uniq_dend_cmap(6, :), 2)), 1); % neutral

for iter = 1:size(h, 1)
    h(iter).Color = mod_dend_cmap(iter, :);
end

xtickangle(70);
yline(.3, 'k--', 'LineWidth', 1.5);
xlabel("Subject index", 'FontSize', 15);
ylabel("Distance metric (Euclidean)", 'FontSize', 15);
title("Dendogram (threshold at 0.3)", "FontSize", 18);

% scatter
subplot(2, 3, [2 3 5 6]); hold on;
for iter = [5 4 1 6 7 2 3]
    scatter3(fixed(1) + e_2_s(C == iter),...
             fixed(2) + e_3_s(C == iter),...
             fixed(3) + off_s(C == iter), 150, 'filled',...
             'MarkerFaceColor', grp_cmap(iter,:),...
             'MarkerFaceAlpha', .6,...
             'MarkerEdgeColor', 'k',...
             'MarkerEdgeAlpha', .6);
    
    temp_i = (-1 .* ones(1, sum(C == iter)))  .^ (1:sum(C == iter));
    text(fixed(1) + e_2_s(C == iter) + temp_i'.*.02,...
         fixed(2) + e_3_s(C == iter) + temp_i'.*.02,...
         fixed(3) + off_s(C == iter) + temp_i'.*.03,...
         sublist_idx(C' == iter),...
         'FontSize', 10, 'Color', [1 1 1].*.4);
end
% [y, z] = meshgrid(-2:.1:2); x = zeros(size(y, 1));
% surf(x, y, z, 'FaceColor', 'b', 'FaceAlpha', .1, 'EdgeColor', 'none');
% [x, z] = meshgrid(-2:.1:2); y = zeros(size(x, 1));
% surf(x, y, z, 'FaceColor', 'b', 'FaceAlpha', .1, 'EdgeColor', 'none');
xlim([-.2 .6]); ylim([-.4 .25]); zlim([-.1 1]);
ax = gca; ax.XAxis.FontSize = 12; 
ax.YAxis.FontSize = 12;
grid('on');
xlabel("Coefficient : Happy", 'FontSize', 15, 'Rotation', 6);
%ylabel("Coefficient : Disgust", 'FontSize', 15, 'Rotation', 340);
zlabel("Coefficient : Offer", 'FontSize', 15);

title("Clustered scatter in coefficient space", "FontSize", 18);

% rationality triangle
subplot(2, 3, 4); hold on;
temp_off_emot_coeff = [e_2_s...
                       off_s];
for iter = [5 4 1 6 7 2 3]
    temp_mean = mean(temp_off_emot_coeff(C == iter,:), 1);
    temp_radius = sum(C == iter)*200;

    scatter(fixed(1) + temp_mean(1),...
            fixed(3) + temp_mean(2),...
            temp_radius, "filled",...
            'MarkerFaceColor', grp_cmap(iter,:),...
            'MarkerFaceAlpha', .6,...
            'MarkerEdgeColor', 'k',...
            'MarkerEdgeAlpha', .6);
end
xline(0, 'k-', 'LineWidth', 2);
%grid('on');
xlim([-.25 .55]); ylim([0 .9]);
ax = gca; ax.XAxis.FontSize = 12; 
ax.YAxis.FontSize = 12;
legend(["Strictly emotion-driven"...
        "Moderately emotion-driven"...
        "Generic rational"...
        "Moderately offer-driven"...
        "Strictly offer-driven"], 'FontSize', 13, "Box", "off",...
        "TextColor", [1 1 1].*.2);
xlabel("Coefficient : Happy", 'FontSize', 15);
ylabel("Coefficient : Offer", 'FontSize', 15); 
title("The Rationality Triangle", 'FontSize', 18);

%% 

figure; hold on;
for iter = unique(C)'
    scatter(e_2(C == iter),...
         off(C == iter),...
         50, 'filled');
    
    text(e_2(C == iter) + .02,...
         off(C == iter),...
         sublist_idx(C' == iter),...
         'FontSize', 8);
end
grid('on');
xline(0, 'k-', 'LineWidth', 2);
yline(0, 'k-', 'LineWidth', 2);
xlabel("e_2 : Disgusted", 'FontSize', 13); 
ylabel("offer", 'FontSize', 13);

%% generate the rationality triangle
figure; hold on;
% temp = glm.Coefficients.Estimate;
% temp_hap_coeff = temp(2); temp_off_coeff = temp(4);
% temp_off_emot_coeff = [e_2_s + temp_hap_coeff...
%                        off_s + temp_off_coeff];
temp_off_emot_coeff = [e_2_s'...
                       off_s'];

for iter = unique(C)'
    temp_mean = mean(temp_off_emot_coeff(C == iter,:), 1);
    temp_radius = sum(C == iter)*300;

    scatter(temp_mean(1), temp_mean(2), temp_radius, "filled",...
            'MarkerEdgeColor', 'k',...
            'MarkerEdgeAlpha', .4);
end

xline(0, 'k-', 'LineWidth', 2); xlabel("Coefficient: Happy");
ylabel("Coefficient: Offer"); xlim([-.3 .55]);
fontsize(15, 'points');
grid('on');

%% okie another scatter: intercept, offer, emot-2
int_s = double(stat(1:40,4));
interc = int_s ./ (max(int_s) - min(int_s));
figure;
% scatter3(glm.Coefficients.Estimate(1) + int_s,...
%          glm.Coefficients.Estimate(2) + e_2_s,...
%          glm.Coefficients.Estimate(4) + off_s, 60, 'filled');
scatter3(interc, e_2, off, 60, 'filled');
text(interc, e_2, off, sublist_idx);

xlabel("intercept"); ylabel("e_2"); zlabel("off");

%% intercept - offer scatter, grouped by prev. cluster

figure; 

subplot(1, 2, 1);
for iter = unique(C)'
    scatter(glm.Coefficients.Estimate(2)+e_2_s(C == iter),...
            glm.Coefficients.Estimate(4)+off_s(C == iter),...
            80, 'filled', 'DisplayName', sprintf('Cluster %d', iter)); hold on;
end
legend('FontSize', 8);
xlabel("e_2: Happy"); ylabel("Offer");
grid('on');
yline(0, 'k--', 'LineWidth', 1.5);

subplot(1, 2, 2);
for iter = unique(C)'
    % scatter(glm.Coefficients.Estimate(2)+e_2_s(C == iter)+e_2_s(C == iter),...
    %         interc_fixed + int_s(C == iter),...
    %         80, "filled", 'DisplayName', sprintf('Cluster %d', iter)); hold on;
    scatter(interc_fixed + int_s(C == iter),...
            glm.Coefficients.Estimate(4)+off_s(C == iter),...
            80, "filled",...
            'DisplayName', sprintf('Cluster %d', iter)); hold on;
end
legend('FontSize', 8);
xlabel("Intercept: predisposition to acceptance"); ylabel("Offer");
xline(0, 'k--', 'LineWidth', 1.5);
yline(0, 'k--', 'LineWidth', 1.5);
grid('on');

fontsize(13, 'points');
%%

figure; hold on;
temp_int = glm.Coefficients.Estimate(1) + int_s;

%% grouping subjects based on regression coefficients : rt mdl
% extract coeff
[~, ~, stat] = randomEffects(lm);
e_2_s = double(stat(161:2:240,4));
e_3_s = double(stat(162:2:240,4));
off_s = double(stat(81:120,4));

e_2 = e_2_s ./ (max(e_2_s) - min(e_2_s));
e_3 = e_3_s ./ (max(e_3_s) - min(e_3_s));
off = off_s ./ (max(off_s) - min(off_s));

% plot
% figure;
% scatter3(e_2_s, e_3_s, off_s, 'filled'); 
% text(e_2_s+.001, e_3_s+.001, off_s+.001, sublist);
% xlabel("e_2"); ylabel("e_3"); zlabel("off");
figure; hold on;
for iter = unique(C)'
    scatter3(e_2(C == iter),...
             e_3(C == iter),...
             off(C == iter), 80, 'filled');
    
    text(e_2(C == iter) + .05,...
         e_3(C == iter) + .05,...
         off(C == iter) + .05,...
         sublist_idx(C' == iter),...
         'FontSize', 8);
end
grid('on');
xlabel("e_2 : Happy", 'FontSize', 13); 
ylabel("e_3 : Disgusted", 'FontSize', 13); 
zlabel("offer", 'FontSize', 13);

%% fit the multilevel model for each subject
% n_sub = length(alldat);
% r_sq = zeros(n_sub, 2);
% mdl_rt_sub = cell(n_sub,1);
% mdl_resp_sub = cell(n_sub, 1);
% 
% 
% for sub_idx = 1:n_sub
%     dat = alldat{sub_idx};
% 
%     % init parameters
%     bID = reshape(repmat(1:20, 18, 1), 1, []);
%     tbl = table(dat(:,1), dat(:,2)./max(dat(:,2)),...
%             (dat(:,2)./max(dat(:,2))).^2,...
%             bID', log(dat(:,4)),...
%             zeros(size(bID')), dat(:,7),...
%             'VariableNames', ["emot" "off" "off_2" "bID"...
%                               "rt" "rt_res" "accpt"]);
%     tbl.emot = categorical(tbl.emot);
%     tbl = tbl(tbl.rt > log(.25), :);
%     tbl = tbl(tbl.accpt ~= 0, :); tbl.accpt = (tbl.accpt + 1)/2;
% 
%     % hierarchial modelling
%     lm = fitlme(tbl, ['rt ~ 1 +'...
%                         'off + emot + off_2 +'...
%                         '(off-1|bID) +'...
%                         '(off_2-1|bID) +'...
%                         '(emot-1|bID)']);
%     mdl_rt_sub{sub_idx} = lm;
%     tbl.rt_res = tbl.rt - lm.fitted;
% 
%     glm = fitglme(tbl, ['accpt ~ rt_res + emot + off +' ...
%                                 '(rt_res-1|bID) +'...
%                                 '(emot-1|bID) +'...
%                                 '(off-1|bID)'],...
%                     'Distribution', 'Binomial',...
%                     'Link', 'logit', 'StartMethod', 'default');
%     mdl_resp_sub{sub_idx} = glm;
% 
%     r_sq(sub_idx, :) = [lm.Rsquared.Adjusted glm.Rsquared.Adjusted];
% end

%% calculate thresh after i-th trial
thresh = {[], [], []};
d_blk = 4;
for idx = 1:d_blk:20
    slice_idx = ((idx-1)*18 + 1):...
                ((idx+d_blk-1)*18);
    for emot_idx = 1:3
        temp = calcThreshold(alldat,emot_idx,slice_idx);
        temp = repmat(temp, 1, d_blk * 18);
        thresh{emot_idx} = horzcat(thresh{emot_idx},...
                                   temp);
    end
end

thresh_complied = zeros(40, 360);
for sub_idx = sublist_idx
    dat = alldat{sub_idx};
    for idx = 1:360
        thresh_complied(sub_idx, idx) = ...
            thresh{dat(idx,1)}(sub_idx, idx);
    end
end

%% check offer vs rt 
figure; 

k = 1;
for sub_idx = 1:40
    dat = tbl(double(tbl.subID) == sub_idx,:);
    subplot(5, 8, k); title(sublist(sub_idx)); k = k + 1;
    hold on;
    
    % temp = thresh_complied(sub_idx,:)' - dat(:,2);
    % temp = temp ./ (max(temp)-min(temp));

    for resp = [1 0]
        % temp_x = temp(dat(:,end) == resp);
        % temp_y = log(dat(dat(:,end)==resp,4));
        temp_x = dat.off(dat.accpt == resp);
        temp_y = dat.rt_res(dat.accpt == resp);
        scatter(temp_x, temp_y, '.');
    end
end

%% speed accuracy trade-off
figure; 

rt_sample = [];
accpt_sample = [];
rt_fit_sample = []; 
accpt_fit_sample = [];
rt_fitted = temp_lm.fitted;
accpt_fitted = temp_glm.fitted > .5;
d_blk = 4;

% iter = 1000; n_samp = size(tbl,1);
% for idx = 1:iter
%     slice_idx = randsample(n_samp, round(100));
% 
%     accpt_sample = cat(2, accpt_sample,...
%                 sum(tbl.accpt(slice_idx)/...
%                 numel(slice_idx)));
%     rt_sample = cat(2, rt_sample,...
%                 mean(tbl.rt(slice_idx))); 
% 
%     accpt_fit_sample = cat(2, accpt_fit_sample,...
%                 sum(accpt_fitted(slice_idx)/...
%                 numel(slice_idx)));
%     rt_fit_sample = cat(2, rt_fit_sample,...
%                 mean(rt_fitted(slice_idx))); 
% end


for sub_idx = setdiff(1:40, [])
    for blk_idx = 1:20
        slice_idx = double(tbl.subID) == sub_idx &...
                    double(tbl.bID) == blk_idx;
        rt_sample = cat(2, rt_sample,...
                        median(tbl.rt(slice_idx)));
        accpt_sample = cat(2, accpt_sample,...
                        sum(tbl.accpt(slice_idx)==1)...
                        /sum(slice_idx));

        rt_fit_sample = cat(2, rt_fit_sample,...
                        median(rt_fitted(slice_idx)));
        accpt_fit_sample = cat(2, accpt_fit_sample,...
                        sum(accpt_fitted(slice_idx)==1)...
                        /sum(slice_idx));
    end
end

scatter(accpt_sample, rt_sample, 'filled', 'MarkerFaceAlpha', .4);
hold on;
scatter(accpt_fit_sample, rt_fit_sample, 'filled',...
        'MarkerFaceAlpha', .4);

%% check rt against attractiveness rating
% load vars
rt = face_Specific_Accept(alldat, 4, 'photIDbehavDat.txt');
attr = importdata('attractiveness_rating.txt');
% normalize attr as (attr - min)/range
attr_norm = zeros(size(attr));
for sub_idx = 1:size(attr,1)
    t = attr(sub_idx, :);
    attr_norm(sub_idx,:) = (t - min(t))/(max(t) - min(t));
end
clear t;

figure;
for sub_idx = 1:size(attr,1)
    subplot(5, 8, sub_idx);
    scatter(attr(sub_idx,:),...
            rt(sub_idx,:), 'filled');
    
end

%% model for attractiveness vs rt or acceptance 
% load vars
rt = face_Specific_Accept(alldat, 4, 'photIDbehavDat.txt');
accpt = face_Specific_Accept(alldat, 7, 'photIDbehavDat.txt');
attr = importdata('attractiveness_rating.txt');
% normalize attr as (attr - min)/range
attr_norm = zeros(size(attr));
for sub_idx = 1:size(attr,1)
    t = attr(sub_idx, :);
    attr_norm(sub_idx,:) = (t - min(t))/(max(t) - min(t));
end
clear t;
face_details = importdata('photIDbehavDat.txt');
face_details = sortrows(face_details, 4);

% prep vars
attr_flat = reshape(attr_norm', [], 1);
subID = reshape(repmat(1:size(attr,1),...
                size(attr,2), 1), [], 1);
faceID = reshape(repmat(1:size(attr,2),...
                size(attr,1), 1), [], 1);
%emot = repmat(face_details(:,1), size(attr,1), 1);
rt = reshape(rt', [], 1);
accpt = reshape(accpt', [], 1);
tbl = table(attr_flat, rt, accpt, faceID, subID,...
            'VariableNames', ["attr" "rt" "accpt" "faceID" "subID"]);


mdl_accpt = fitlme(tbl, "accpt ~ attr + (attr|subID) + (attr|faceID)");
% note: acceptance ~ attr : mdl_accpt; rt ~ attr : mdl_rt














































































































%% LPP and N170 in emotion identification
lpp = load('lpp_trial_data.mat'); lpp = lpp.lppHeights_19235657_dat;
n170 = load('n170_trial_data.mat'); n170 = n170.NpeakHeights_19235657_dat;
frn = load('FRN_trial_data_peakjalli.mat'); frn = frn.NpeakHeights_3_dat;
alldat_mod = load('behavDat_collated_final.mat'); alldat_mod = alldat_mod.alldat;

% pull the emot indices
emot_pulled = [];
off_pulled = [];
subID = [];
for iter_sub = 1:40
    dat = alldat_mod{iter_sub};

    emot_pulled = [emot_pulled dat(:,1)'];
    off_pulled = [off_pulled dat(:,2)'];
    subID = [subID iter_sub .* ones(1, size(dat,1))];
end
% pull N170 & LPP vals
n170_pulled = []; lpp_pulled = [];
for iter_sub = 1:40
    for iter_cond = 1:9
        n170_pulled = [n170_pulled...
                n170{iter_sub}{iter_cond}];
        lpp_pulled = [lpp_pulled...
                lpp{iter_sub}{iter_cond}];
    end
end

clear iter_cond iter_sub dat lpp n170;

%% pull & prepare the FRN data
figure;
fn = "FRN_trial_data_peakjalli.mat";

k = 1;
for iter_fn = fn
    frn = load(iter_fn);
    frn = frn.NpeakHeights_3_dat;
    
    frn_pulled = []; 
    for iter_sub = 1:40
        temp = [];
        for iter_cond = 1:9
            temp = [temp...
                    frn{iter_sub}{iter_cond}];
        end
        frn_pulled = [frn_pulled...
                    normalize(temp)];
    end
    
    off_uniq = unique(off_pulled);
    
    frn_off = [];
    for iter_off = off_uniq
    frn_off(end + 1) = mean(frn_pulled(off_pulled == iter_off));
    end
    
    subplot(1, numel(fn), k); k = k+1;
    scatter(off_uniq, frn_off, 'filled');
    title(iter_fn);

end

clear iter_cond iter_sub iter_off iter_fn k;
%% normalization across subjects
subID_uniq = unique(subID);

n170_pulled_norm = zeros(size(n170_pulled));
lpp_pulled_norm = zeros(size(lpp_pulled));
frn_pulled_norm = zeros(size(frn_pulled));

for iter_sub = subID_uniq 
    n170_pulled_norm(subID == iter_sub) = normalize(...
                         n170_pulled(subID == iter_sub));
    lpp_pulled_norm(subID == iter_sub) = normalize(...
                         lpp_pulled(subID == iter_sub));
    frn_pulled_norm(subID == iter_sub) = normalize(...
                         frn_pulled(subID == iter_sub));
end

%% FRN - offer modelling
% center and standardize the data
off_uniq_std = off_uniq./max(off_uniq) -...
               prctile(off_uniq./max(off_uniq), 50);
iter_cv = 100; % CV to be repeated
fold_cv = 5; % 5-fold CV

poly_deg = ["poly1" "poly2" "poly3"];
pred_rmse = zeros(iter_cv, numel(poly_deg)); % to store MSE
fit_r2 = zeros(iter_cv, numel(poly_deg)); % to store adj R2 

for k = 1:iter_cv
    c = cvpartition(numel(off_uniq_std), "KFold", fold_cv);
    % train set
    train_idx = training(c, fold_cv);
    off_train = off_uniq_std(train_idx); frn_train = frn_off(train_idx);
    % test set
    test_idx = test(c, fold_cv);
    off_test = off_uniq_std(test_idx); frn_test = frn_off(test_idx);
    
    % fit
    for iter_deg = 1:numel(poly_deg)
        [p, gof] = fit(off_train', frn_train', poly_deg(iter_deg));
        fit_r2(k, iter_deg) = gof.adjrsquare;
        rmse_err = sum((frn_test' - ...
                  feval(p, off_test)).^2)/numel(frn_test);
        pred_rmse(k, iter_deg) = rmse_err;
    end
end

clear iter_cv test_idx train_idx off_train off_test...
    frn_train frn_test poly_deg rmse_err;

% box plot of the statistics
subplot(1, 2, 1);
boxchart(fit_r2, 'MarkerStyle', 'none'); 
title("adj R^2 for different pol-deg");
xlabel("deg of pol fit");

subplot(1, 2, 2);
boxchart(pred_rmse, 'MarkerStyle', 'none'); 
title("RMSE pred for different pol-deg");
xlabel("deg of pol fit");
























