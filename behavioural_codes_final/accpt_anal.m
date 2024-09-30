%% check offer distribution across emotion
figure;
emot_idx_array = [];
offer_array = [];

for sub_idx = sublist_idx
    dat = alldat{sub_idx};
    emot_idx_array = horzcat(emot_idx_array,...
                     dat(:,1)');
    offer_array = horzcat(offer_array,...
                     dat(:,2)');
end

boxchart(emot_idx_array, offer_array);

%% Acceptance - attractiveness regression
% assuming ratio scale data
% standardizing ordinal data
accpt = face_Specific_Accept(alldat, 7, 'photIDbehavDat.txt');
attr = importdata('attractiveness_rating.txt');
attr_std = zeros(size(attr));
for face_idx = 1:size(attr, 2)
    temp = attr(:,face_idx);
    % attr_std(:, face_idx) = (temp - median(temp))./...
    %                         (max(temp) - min(temp));
    attr_std(:, face_idx) = temp - median(temp); 
end

% face ID
face_details = importdata('photIDbehavDat.txt');
face_details = sortrows(face_details, 4);
faceID = reshape(repmat(face_details(:,end)',...
                                40, 1), [], 1);
emotID = reshape(repmat(face_details(:,1)',...
                                40, 1), [], 1);

% flatten "attr" and "accpt" matrixes (40 X 24)
% and construct the table
tbl = table(reshape(attr, [], 1),...
            reshape(accpt, [], 1),...
            faceID,...
            emotID,...
            'VariableNames', ["attr" "accpt" "faceID" "emotID"]);
tbl.attr = categorical(tbl.attr);
tbl.emotID = categorical(tbl.emotID);
tbl.faceID = categorical(tbl.faceID);

% construct the linear model
lm = fitlme(tbl, 'accpt ~ attr + (attr|faceID) + (attr|emotID)');

clear temp;

%% visualization
figure; 
subplot(1,4,1); hold on;
for face_idx = 1:size(attr,2)
    scatter(tbl.attr(faceID==face_idx),...
            tbl.accpt(faceID==face_idx), 'filled');
end
title('Scatter, faceID coded');
xlabel('attr'); ylabel('accpt');

subplot(1,4,2);
plotPartialDependence(lm, "attr");
ylim([0 1]);

subplot(1,4,3);
plotPartialDependence(lm, "emotID");
ylim([0 1]);

subplot(1,4,4);
plotPartialDependence(lm, "faceID");
ylim([0 1]);

sgtitle("model assum categ attr");

%% attractiveness - modelling wrt emotion
% load vars
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
attr_flat = reshape(attr', [], 1);
subID = reshape(repmat(1:size(attr,1), 24, 1), [], 1);
emot = repmat(face_details(:,1), size(attr,1), 1);
tbl = table(attr_flat, emot, subID,...
            'VariableNames', ["attr" "emot" "subID"]);
tbl.emot = categorical(tbl.emot);

mdl = fitglme(tbl, 'attr ~ emot + (emot|subID)',...
              'Distribution', 'Poisson');

%% visualization
disp(mdl.Rsquared.Adjusted);
figure; 

[~, ~, stat] = randomEffects(mdl);
temp = mdl.Coefficients.Estimate;
c_2 = temp(2) + double(stat(2:3:end,4));
c_3 = temp(2) + double(stat(3:3:end,4));

subplot(1,2,1);
scatter(c_2, c_3, 'filled'); hold on;
text(c_2, c_3, sublist);
xlabel("emot 2"); ylabel("emot 3");

subplot(1,2,2);
mdl.plotPartialDependence("emot");
title(sprintf("p_2: %.2f, p_3: %.2f",...
      mdl.Coefficients.pValue(2),...
      mdl.Coefficients.pValue(3)));


%% Cond specific Acceptance over time
% clc;
% iter = 2;
% init_trial = 0:iter*18:360-iter*18;
% time = 1:length(init_trial);
% 
% % hap/dis/neu: mean across subjects
% % _dev: standard error
% hap = zeros(4, length(init_trial)); hap_dev = zeros(4, length(init_trial));
% dis = zeros(4, length(init_trial)); dis_dev = zeros(4, length(init_trial));
% neu = zeros(4, length(init_trial)); neu_dev = zeros(4, length(init_trial));
% 
% 
% for i = time
%     trial = init_trial(i);
% 
%     source_cell = behavDataSummary(alldat, sublist, trial+1, trial+iter*18, 7);
%     source_cell(AAsubj) = [];
%     [temp_accpt, temp_accpt_sub] = prep_accpt_matrix(source_cell);
% 
%     for off = 1:3
%         hap(off,i) = temp_accpt(1, off);
%         hap_dev(off, i) = std(temp_accpt_sub{1, off});
% 
%         dis(off,i) = temp_accpt(2, off);
%         dis_dev(off, i) = std(temp_accpt_sub{2, off});
% 
%         neu(off,i) = temp_accpt(3, off);
%         neu_dev(off, i) = std(temp_accpt_sub{3, off});
%     end
%     hap(4,i) = mean(temp_accpt(1, :)); 
%     hap_dev(4, i) = std([temp_accpt_sub{1, 1} ...
%                            temp_accpt_sub{1, 2} ...
%                            temp_accpt_sub{1, 3}]);
% 
%     dis(4,i) = mean(temp_accpt(2, :));
%     dis_dev(4, i) = std([temp_accpt_sub{2, 1} ...
%                            temp_accpt_sub{2, 2} ...
%                            temp_accpt_sub{2, 3}]);
% 
%     neu(4,i) = mean(temp_accpt(3, :));
%     neu_dev(4, i) = std([temp_accpt_sub{3, 1} ...
%                            temp_accpt_sub{3, 2} ...
%                            temp_accpt_sub{3, 3}]);
% end
% 
% figure;
% tiledlayout(2, 2);
% 
% for off = 1:4
%     %subplot(1,4,off);
%     nexttile;
% 
%     % plot the mean trace
%     plot(time, hap(off,:),...
%         time, dis(off,:),...
%         time, neu(off,:)); hold on;
% 
%     % plot the se
%     % patch([time fliplr(time)],...
%     %     [hap(off,:)-hap_dev(off,:)  fliplr(hap(off,:)+hap_dev(off,:))],...
%     %     'b', "FaceAlpha", 0.2, "EdgeColor", "none");
%     % patch([time fliplr(time)],...
%     %     [dis(off,:)-dis_dev(off,:)  fliplr(dis(off,:)+dis_dev(off,:))],...
%     %     'r', "FaceAlpha", 0.2, "EdgeColor", "none");
%     % patch([time fliplr(time)],...
%     %     [neu(off,:)-neu_dev(off,:)  fliplr(neu(off,:)+neu_dev(off,:))],...
%     %     'y', "FaceAlpha", 0.2, "EdgeColor", "none");
%     hold off
% 
%     xticks(time); grid("on");
%     xlabel("Segment"); ylabel("Acceptance");
%     legend(["Hap" "Dis" "Neu"], "Location","southeast");
% 
%     if off ~= 4
%         title("off: "+num2str(off + 1));
%     else
%         title("mean across off-s");
%     end
% end
%% Distribution of Acceptance rate (3X3) + heatmap table
condSpecAccpt = behavDataSummary(alldat,1,360,7);
AA_idx = [12 14 15 21 29 37:39];
[accpt_mat, ~, ~] = prep_accpt_matrix(condSpecAccpt, []);

n_subj = length(alldat) - numel(AA_idx);
emots = 1:3; offers = 2:4;
plot_index = 1;
figure;

condMatrix = zeros(9, n_subj); index = 1;

for emot = emots
    for offer = offers
        subj_temp = zeros(1,n_subj);
        cell_index = emot*10 + offer;
        
        iter_subj = 1;
        for subj = setdiff(1:n_subj, AA_idx)
            temp = condSpecAccpt{subj}{cell_index}; %change for RT
            subj_temp(iter_subj) = sum(temp == 1)/length(temp);
            iter_subj = iter_subj + 1;
            %subj_temp(subj) = mean(temp); %for RT
        end

        % store the accpt rate values
        condMatrix(index, :) = subj_temp;
        index = index + 1;
        
        if offer == 2
            align = 'right';
        else
            align = 'left';
        end
        
        % last two cols reserved for the heatmap
        subplot(3,5,plot_index); 
        if plot_index == 3 || plot_index == 8
            plot_index = plot_index + 3;
        else
            plot_index = plot_index + 1;
        end


        histogram(subj_temp, 7); hold on; grid minor;
        xline(mean(subj_temp), 'k--',...
            'LabelHorizontalAlignment',align,...
            'LabelOrientation','horizontal',...
            'Label', {"\mu="+num2str(mean(subj_temp),2)});
        title(num2str(cell_index));
    end
end

% heatmap table now

figure;
heatmap({'2:8', '3:7', '4:6'},...
        {'H', 'D', 'N'},...
        round(accpt_mat, 2));
colormap sky;
xlabel('Offer'); ylabel('Emotion');
fontsize(20, "points");
grid off;

clear emots offers plot_index index...
    emot offer subj n_subj temp accpt_mat;

%% t-test after eliminating AA/AR
space = 360;
start_block = 1:space:360;
% AA_idx = [14 15 21 29 37:39]; -- not excluding these => better p
AA_idx = [];

emot_label = ["hap" "dis" "neu"];
emot_series = [3 2; 3 2; 1 3]; off = 2;


for emot_idx = 1:3
emot = emot_series(emot_idx,:);
figure; k = 1;
for iter = start_block
    condSpecAccpt = behavDataSummary(alldat,iter,iter+space-1,7);
    [~, subCondAccpt, ~] = prep_accpt_matrix(condSpecAccpt,...
                                              AA_idx);
    temp_hap = subCondAccpt{emot(1),off-1};
    temp_dis = subCondAccpt{emot(2),off-1};

    subplot(numel(start_block), 2, k); 
    % pulled hist
    histogram(temp_hap, 0:.2:1); hold on;
    histogram(temp_dis, 0:.2:1);
    legend([emot_label(emot(1)) emot_label(emot(2))],...
        "Location", "northwest");
    xlabel("accpt rate");
    title(sprintf("%d - %d", iter, iter+space-1));
    grid('on');

    subplot(numel(start_block), 2, k+1);
    % ttest
    histogram(temp_hap - temp_dis, -.3:.05:.3);
    [~, p] = ttest((temp_hap - temp_dis), 0, 'Tail', 'Right');
    title(sprintf("p = %.3f", p));
    xlabel("diff in accpt rate");
    xlim([-.3 .3]);
    grid('on');

    k = k + 2;

end
if numel(AA_idx) > 0
    sgtitle(sprintf("ttest for accpt rate against %s & %s (excld. AAR)",...
                    emot_label(emot(1)), emot_label(emot(2))));
else
    sgtitle(sprintf("ttest for accpt rate against %s & %s (incld. AAR)",...
                    emot_label(emot(1)), emot_label(emot(2))));
end
end

clear temp_dis temp_hap k p space;

%% friedman test
%NOT NEEDED
% sub_accpt_mat = test_friedman(alldat, sublist);
% % note: sub_accpt_mat = emot X off
% 
% % effect of offer across emotions
% n_subj = 34;
% disp("effect of offers for part. emotion--------");
% for emot = 1:3
%     sub_accpt_mat_emot = sub_accpt_mat(...
%         (emot-1)*n_subj + 1 : (emot-1)*n_subj + 34, :);
%     disp(emot);
%     [p, ~, stats] = friedman(sub_accpt_mat_emot, 1, "off");
%     disp(p);
%     disp(multcompare(stats));
%     disp("********");
% end
% 
% disp("effect of emotions for part. offer--------");
% for off = 1:3
%     sub_accpt_mat_off = sub_accpt_mat(:, off);
%     sub_accpt_mat_off = reshape(sub_accpt_mat_off,...
%                         n_subj, 3);
%     disp(1 + off);
%     [p, ~, stats] = friedman(sub_accpt_mat_off, 1, "off");
%     disp(p);
%     %disp(multcompare(stats));
%     %disp("********");
% end
% 
% clear sub_accpt_mat_off sub_accpt_mat_emot;

%% rm-anova - now we can check the effect of learning

% calculate sub-wise accept rate, across conds
AA_idx = []; n_subj = 40 - numel(AA_idx);
[condSpecAccpt, ~, blk_idx] = behavDataSummary(alldat,1,360,7);

% old bit
% [~, subCondAccpt, ~] = prep_accpt_matrix(condSpecAccpt,...
%                                           AA_idx);
% 
% % store the data
% emot_label = 1:9; off_label = 1:9;
% subCondAccpt_mat = zeros(n_subj, 9);
% k = 1;
% for emot_idx = 1:3
%     for off_idx = 2:4
%         subCondAccpt_mat(:,k) = subCondAccpt{emot_idx, off_idx-1}';
%         emot_label(k) = emot_idx; off_label(k) = off_idx;
%         k = k + 1;
%     end
% end

% first construct the repeat-mes data matrix and within factors
cluster_n_blocks = 5;
emot_label = [];
off_label = [];
blk_label = [];
rm_accp_mat = zeros(40, 9*(20/cluster_n_blocks));

for sub_idx = 1:40
    
    temp_accp_holder = [];

    for emot_idx = 1:3
        for off_idx = 2:4
            key = 10*emot_idx + off_idx;

            % slice accpt data
            temp_accpt_array = condSpecAccpt{sub_idx}{key};
            temp_blk_idx = blk_idx{sub_idx}{key};

            % store the rt, emot, off & blk idx data
            blk_iter = 1;
            for start_block = 1:cluster_n_blocks:...
                                (21-cluster_n_blocks)
                slice_idx = temp_blk_idx >= start_block &...
                            temp_blk_idx < start_block +...
                            cluster_n_blocks;

                accpt_rate = sum(temp_accpt_array(...
                            slice_idx) == 1)...
                            /...
                            numel(temp_accpt_array(...
                            slice_idx));
                temp_accp_holder = horzcat(temp_accp_holder,...
                                            accpt_rate);

                % emot and offer idx
                if sub_idx == 1
                    emot_label(end+1) = emot_idx;
                    off_label(end+1) = off_idx;
                    blk_label(end+1) = blk_iter;
                    blk_iter = blk_iter + 1;
                end
            end
        end
    end

    rm_accp_mat(sub_idx, :) = temp_accp_holder;
end

% construct the table
tbl_sub = table((1:n_subj)', 'VariableNames', "sub");
tbl_sub.sub = categorical(tbl_sub.sub);
y = rm_accp_mat;
tbl_dat = array2table(y);
tbl = [tbl_sub tbl_dat];

within_des = table(emot_label', off_label', blk_label',...
    10*emot_label' + blk_label', 10*off_label' + blk_label',... 
    10*emot_label' + off_label',...
    'VariableNames', ["emot", "off", "blk", "e_b", "o_b", "e_o"]);
within_des.emot = categorical(within_des.emot);
within_des.off = categorical(within_des.off);
within_des.blk = categorical(within_des.blk);
within_des.e_b = categorical(within_des.e_b);
within_des.o_b = categorical(within_des.o_b);

% rm anova
variable_str = sprintf("y1-y%d~1", size(y,2));
rm = fitrm(tbl, variable_str, WithinDesign = within_des);
rm_mdl = ranova(rm, 'WithinModel',...
                'emot+off+blk+emot*blk+off*blk+emot*off');
disp(rm_mdl);

clear accpt_rate tbl_sub tbl_dat emot_idx off_idx y;

%% set colour map
cmap = parula(21);
cmap_off = cmap([1 4 10],:);
cmap_emot = cmap([13 17 19], :);

%% generate the anova figure;
figure;

% emot
subplot(2, 3, 1); 
[temp_mu, temp_se] = rm_anova_mean_se(tbl,...
                    double(within_des.emot));
for idx = 1:3
    errorbar(idx, temp_mu(idx), temp_se(idx), '-',... 
         'Color', cmap_emot(idx,:),...
         'LineWidth', 3, 'CapSize', 15,...
         'Marker', 'o',...
         'MarkerSize', 8,...
         'MarkerFaceColor', [.8 .8 .8],...
         'MarkerEdgeColor', 'w');
    hold on;
end
% errorbar(1:3, temp_mu, temp_se, 'k-',... 
%          'LineWidth', 2, 'CapSize', 8,...
%          'Marker', 'o',...
%          'MarkerFaceColor', [.8 .8 .8],...
%          'MarkerEdgeColor', 'w');
xticks(1:4);
xticklabels(["Happy" "Disgust" "Neutral"])
xlim([0 4]); 
%ylim([0.61 0.68]);
ylim([-.38 -.3]);
ylabel('Acceptance');
set(gca, 'FontSize', 15);
title('Emotion', 'FontSize', 18); 
set(gca, 'box', 'off');

% off
subplot(2, 3, 4);
[temp_mu, temp_se] = rm_anova_mean_se(tbl,...
                    double(within_des.off));
for idx = 1:3
    errorbar(idx, temp_mu(idx), 1.5.*temp_se(idx), '-',... 
         'Color', cmap_off(idx,:),...
         'LineWidth', 3, 'CapSize', 15,...
         'Marker', 'o',...
         'MarkerSize', 8,...
         'MarkerFaceColor', [.8 .8 .8],...
         'MarkerEdgeColor', 'w');
    hold on;
end
xticklabels(["Low" "Int" "Max"])
xlim([0 4]); 
%ylim([.3 1]);
ylim([-.42 -.26]);
ylabel('Acceptance');
set(gca, 'FontSize', 15);
title('Offer', 'FontSize', 18); 
set(gca, 'box', 'off');

% emot * off
subplot(2, 3, 2);
[temp_mu, temp_se] = rm_anova_mean_se(tbl,...
                    double(within_des.e_o));
emot_idx = 1;
for idx = 1:3:9
    slice_idx = idx:(idx+2);
    plot(slice_idx, temp_mu(slice_idx),...
         "LineWidth", 3, "Color", cmap_emot(emot_idx,:));
    hold on;
    emot_idx = emot_idx  + 1;
    for off_idx = 1:3
        errorbar(slice_idx(off_idx),...
             temp_mu(slice_idx(off_idx)),...
             temp_se(slice_idx(off_idx)), '-',...
             'Color', cmap_off(off_idx,:),...
             'LineWidth', 2, 'CapSize', 8,... 
             'Marker', 'o',...
             'MarkerFaceColor', [.8 .8 .8],...
             'MarkerEdgeColor', 'w');
        hold on;
    end
end
xticks(1:9);
xticklabels(["Low" "Int" "Max"...
             "Low" "Int" "Max"...
             "Low" "Int" "Max"]);
xlim([0 10]); 
%ylim([.15 1]);
ylim([-.43 -.25]);
ylabel('Acceptance');
set(gca, 'FontSize', 15);
title('Emotion x Offer', 'FontSize', 18); 
set(gca, 'box', 'off');

% blk
subplot(2, 3, 5);
[temp_mu, temp_se] = rm_anova_mean_se(tbl,...
                    double(within_des.blk));
errorbar(1:4, temp_mu, temp_se, '--',... 
         'Color', .3.*[1 1 1],...
         'LineWidth', 3, 'CapSize', 15,...
         'Marker', 'o',...
         'MarkerFaceColor', [.8 .8 .8],...
         'MarkerEdgeColor', 'w');
xticklabels(["1-5", "6-10", "11-15", "16-20"]);
xlim([0 5]); 
%ylim([.56 .72]);
ylabel("Acceptance");
set(gca, 'FontSize', 15);
title("Block-quarter", 'FontSize', 18); 
set(gca, 'box', 'off');

% blk * emot
subplot(2, 3, 3);
[temp_mu, temp_se] = rm_anova_mean_se(tbl,...
                    double(within_des.e_b));
emot_idx = 1;
for idx = 1:4:12
    slice_idx = idx:(idx+3);
    errorbar(slice_idx, temp_mu(slice_idx),...
             temp_se(slice_idx), '--',...
             'LineWidth', 3, 'CapSize', 8,... 
             'Marker', 'o',...
             'Color', cmap_emot(emot_idx, :),...
             'MarkerFaceColor', [.8 .8 .8],...
             'MarkerEdgeColor', 'w');
    emot_idx = emot_idx + 1;
    hold on;
end
xticks(1:12);
xticklabels(["H1" "H2" "H3" "H4"...
             "D1" "D2" "D3" "D4"...
             "N1" "N2" "N3" "N4"]);
xlim([0 13]); 
ylim([-.46 -.19]);
ylabel('Acceptance');
set(gca, 'FontSize', 15);
title('Emotion x Block', 'FontSize', 18);
set(gca, 'box', 'off');

% blk * off
subplot(2, 3, 6);
[temp_mu, temp_se] = rm_anova_mean_se(tbl,...
                    double(within_des.o_b));
off_idx = 1;
for idx = 1:4:12
    slice_idx = idx:(idx+3);
    errorbar(slice_idx, temp_mu(slice_idx),...
             temp_se(slice_idx), '--',...
             'LineWidth', 3, 'CapSize', 8,... 
             'Marker', 'o',...
             'Color', cmap_off(off_idx, :),...
             'MarkerFaceColor', [.8 .8 .8],...
             'MarkerEdgeColor', 'w');
    off_idx = off_idx + 1;
    hold on;
end
xticks(1:12);
xticklabels(["L1" "L2" "L3" "L4"...
             "I1" "I2" "I3" "I4"...
             "M1" "M2" "M3" "M4"]);
xlim([0 13]); 
ylim([-.5 -.18]);
ylabel('Acceptance');
set(gca, 'FontSize', 15);
title('Offer x Block', 'FontSize', 18); 
set(gca, 'box', 'off');

%% segment-wise Acceptance rate vector across emotion
figure;

for sub_idx = 1:length(alldat)
    dat = alldat{sub_idx};
    subplot(5, 8, sub_idx); hold on;

    for emot_idx = 1:3
        % Acceptance rate broken in n segments
        n_segments = 4;
        slice_idx = reshape(1:size(dat,1),...
                        n_segments, []);
        segment_accpt_rate = zeros(n_segments, 1);

        % calculate rho for each segment
        for segment_idx = 1:n_segments
            thresh = calcThreshold()
            temp = dat(slice_idx(segment_idx,:),:);
            temp = sum(temp(temp(:,1)==emot_idx, end)==1)/...
                       size(temp(temp(:,1)==emot_idx,:), 1);
            %temp = mean(temp(temp(:,1)==emot_idx,2));
            segment_accpt_rate(segment_idx) = temp;
        end

        plot(1:n_segments, segment_accpt_rate, 'o--',...
                           'LineWidth', 1.5);
    end
    title(sublist(sub_idx));
    ylim([0 1]); xlim([1 n_segments])
    xlabel("segment ID"); ylabel("\rho");
end


%% saturation of Acceptance rate
figure;
del = 10; % threshold +- del
k = 1;
n_block = 20;

% pre-calculate threshold
% threshold_matrix_mod = zeros(40,3);
% for emot_idx = 1:3
%     threshold_matrix_mod(:,emot_idx) = calcThresholdmod(...
%                     alldat, emot_idx, 1:n_block*18)';
% end
threshold_matrix = threshold_matrix_mod;

sub_iter =  sublist_idx(~ isnan(rsq));
area_bw_temp = zeros(numel(sub_iter),1);

for sub_idx = sub_iter
    dat = alldat{sub_idx};
    dat = dat(1:18*n_block,:);

    if length(sub_iter) > 1
        subplot(5, 7, k); k = k+1; hold on;
    else
        subplot(1, 1, 1); hold on;
    end
    
    temp = [];
    for emot_idx = 1:2
        slice_idx = dat(:,1) == emot_idx &...
                    dat(:,2) > threshold_matrix(sub_idx,emot_idx)-del &...
                    dat(:,2) < threshold_matrix(sub_idx,emot_idx)+del; 
        
        % slice_idx = dat(:,1) == emot_idx &...
        %             findBaseOff(dat(:,2)) == 2; 
        % slice_idx = dat(:,1) == emot_idx;

        emot_dat = dat(slice_idx,:);
        accpt = emot_dat(:,7) == 1;

        %accpt_local = movmean(accpt, 10);
        accpt = cumsum(accpt) ./ (1:numel(accpt))';
        plot(accpt, 'LineWidth', 1.5); 
        
        %plot(accpt_local, 'r-');

        % calculate "area" b/w curves
        if emot_idx == 1
            temp = accpt;
            N = numel(accpt);
        elseif emot_idx == 2
            N = min(N, numel(accpt));
            temp = temp(1:N) - accpt(1:N);
            temp = sum(temp)/N; 
            area_bw_temp(sub_idx) = temp;
        end
    end
    
    xlim([1 N]); 
    ylim([0 1]);
    title(sprintf("%s | %.1f | %.1f",...
                        sublist(sub_idx),...
                        threshold_matrix(sub_idx,1),...
                        threshold_matrix(sub_idx,2)));
end
sgtitle("Acceptance rate across emot thresh");

clear threshold_temp temp;

%% calculate the threshold
figure;
k = 1;
permitted_error = .02;
store_thresh = zeros(40,3);
n_block = 20;
for sub_idx = 1:40
    dat = alldat{sub_idx};
    accpt_rate = sum(dat(:,end)==1)/size(dat,1);
    subplot(5, 8, k); k = k+1; hold on;

    for emot_idx = 1:3
        slice_idx = dat(1:n_block*18,1) == emot_idx; % <-- emot
        emot_dat = dat(slice_idx,:); n_trial = sum(slice_idx);
        %emot_dat = dat; 
    
        % define threshold as mean(75th per of rej off, 
        %                          25th per of accp off)
        thresh = zeros(n_trial,1);

        for iter = 1:n_trial
            % leave scope for mistakes, approx 7 in 360
            if accpt_rate < permitted_error
                thresh(iter) = prctile(emot_dat...
                                (emot_dat(1:iter,end)==-1,2),...
                                95);
                
            elseif accpt_rate > 1-permitted_error
                thresh(iter) = prctile(emot_dat...
                                (emot_dat(1:iter,end)==1,2),...
                                5);
                
            else
                off_accp_min = prctile(emot_dat...
                                (emot_dat(1:iter,end)==1,2),...
                                25);
                off_rej_max = prctile(emot_dat...
                                (emot_dat(1:iter,end)==-1,2),...
                                75);
                if isnan(off_accp_min)
                    off_accp_min = off_rej_max;
                elseif isnan(off_rej_max)
                    off_rej_max = off_accp_min;
                end
                thresh(iter) = mean([off_accp_min off_rej_max]);
            end
        end
        store_thresh(sub_idx,emot_idx) = mean(thresh(end-5:end));
        plot(thresh, 'LineWidth', 1.5);
    end
    %xlim([1 120]); 
    ylim([15 45]);
    title(sublist(sub_idx));
end
sgtitle(sprintf("threshold (%d)", n_block));

%% bar for "near-threshold" offer, segmented by emot
pulled_behav = [];
for sub_idx = setdiff(1:40, [15 21 29 37:39])
    dat = alldat{sub_idx}; dat = dat(1:n_block*18,:);
    for emot_idx = 1:3
        thresh = round(store_thresh(sub_idx,emot_idx));
        near_thresh = [thresh-5 thresh+5];
        slice_idx = dat(:,1)==emot_idx&...
                    dat(:,2)>=near_thresh(1)&...
                    dat(:,2)<=near_thresh(2);
        temp_dat = dat(slice_idx,:);
        temp_dat(:,2) = temp_dat(:,2) - thresh;
        pulled_behav = vertcat(pulled_behav, temp_dat);
    end
end

accpt_rate_mat = zeros(11, 3);

for off_idx = -5:5
    for emot_idx = 1:3
        slice_idx = pulled_behav(:,1)==emot_idx &...
                    pulled_behav(:,2)==off_idx;
        accpt_rate = pulled_behav(slice_idx, 7);
        accpt_rate = sum(accpt_rate == 1)/numel(accpt_rate);
        accpt_rate_mat(off_idx+6, emot_idx) = accpt_rate;
    end
end

figure; bar(accpt_rate_mat); 
xticklabels(-5:5);
clear temp_dat thresh near_thresh slice_idx

%% redo all the analysis
% Acceptance rate vs offer across emot
figure; hold on;
n_block = 5; sub_idx = 3;
dat = alldat{sub_idx};
dat = dat(1:n_block*18,:);

off_lab = unique(dat(:,2))';

for emot_idx = 1:2
    accp_rate = [];
    for off_val = off_lab
        slice_idx = dat(:,2)==off_val & dat(:,1) == emot_idx;
        accp_rate(end + 1) = sum(dat(slice_idx,end)==1)/...
                             numel(dat(slice_idx,end));
    end
    scatter(off_lab, accp_rate, 'filled');
end
%% start with plain t/rm-anova

% note: block = 4 -> N = 9
% block = 5 -> N = 14

figure;
i = 1;

for start_block = 1:360
    % define llim & ulim of trials
    llim = start_block;
    ulim = start_block+18*5-1;

    % init param
    n_trial = 50;
    sub_iter = setdiff(1:40, sublist_idx(area_bw_20>0));
    accpt_rate = zeros(numel(sub_iter),2);
    
    % pre-calculate threshold
    threshold = [];
    for emot_idx = 1:2
        threshold_temp = calcThresholdmod   (alldat,...
                    [emot_idx 14 45], llim, ulim);
        threshold = horzcat(threshold, threshold_temp);
    end
    
    % iter over sub and store Acceptance values
    k = 1;
    for sub_idx = sub_iter 
        dat = alldat{sub_idx};
        dat = dat(llim:ulim,:);
        
        for emot_idx = 1:2
            slice_idx = dat(:,1) == emot_idx &...
                        dat(:,2) > threshold(sub_idx,emot_idx)-10 &...
                        dat(:,2) < threshold(sub_idx,emot_idx)+10;
            % n_trial, again, is heuristic; subject to change
            temp = dat(slice_idx, 7); temp = temp(1:n_trial);
            accpt_rate(k,emot_idx) = sum(temp==1)/n_trial; 
        end
        k = k+1;
    end
    
    subplot(1, 4, i); i = i+1;
    histogram(accpt_rate(:,1)-accpt_rate(:,2), 5);

    [~, p] = ttest(accpt_rate(:,1)-accpt_rate(:,2));
    title(sprintf("%d - %d | %.3f", llim, ulim, p));
end

clear temp threshold_temp;

%% thresh anal
% invalid_sub = [12 15 21 29 37:40]; % all acceptors

% llim_to_ulim = 1:5*18;
[thresh, ~, ~] = calcThresholdmod(alldat, 0, 1:360);

figure; hold on;
% scatter, sub label and x = y
scatter(dg_off, thresh, 60, 'filled');
text(dg_off + .5, thresh + .5, string(1:40), 'FontSize', 7);
plot(15:50, 15:50, 'k--', 'LineWidth', 2);
plot(15:50, (15:50)./2, 'k--', 'LineWidth', 2);
% quantile reg
[p, S] = quantreg(dg_off, thresh',  .5);
plot(15:50, polyval(p, 15:50), 'r--', 'LineWidth', 2);
% linear reg
% temp_lm = fitlm(dg_off, thresh', 'Intercept', false);
% plot(15:50, feval(temp_lm, 15:50), 'c--', 'LineWidth', 1.5)

xlabel('one-shot DG offer', 'FontSize', 13);
ylabel('repeat-game UG threshold', 'FontSize', 13);

%% new algo to calculate thresh
logistic = fittype('1/(1+exp(-b*(x-c)))',...
'dependent',{'y'},'independent',{'x'},...
'coefficients',{'b', 'c'});
warning('off');
n_block = 20;
x_off = 15:44; med = median(x_off);
x_off = (x_off - med)/numel(x_off);
emot_idx = 2;

for i = 3
    figure;
    
    ax = subplot(2, 2, 1); hold on;
    dat = alldat{i}; 
    
    dat = dat(dat(:,1)==emot_idx,:);

    p_accpt = bar_with_hist(dat(:,end), dat(:,2));
    title('Empirical distr of offers');
    
    subplot(2, 2, 2); hold on;
    plot(15:44, p_accpt, 'k.--', 'LineWidth', 1.5);
    plot(15:44, 1- p_accpt, 'kx--', 'LineWidth', 1.5);
    xlabel('offer'); ylabel('p(A) or p(R)');
    legend(["A" "R"]);
    title('p(accept|offer)');
    
    [f, gof] = fit(x_off', p_accpt', logistic,...
        'startPoint', [10 0],...
        'MaxIter', 800);

    subplot(2, 2, 3); 
    scatter(x_off, p_accpt, 'k', 'filled'); hold on;
    plot(f);
    xlabel('offer'); ylabel('p(Accept)');
    ylim([-.1 1.1]);
    legend(["p vs off" "fitted"]);
    title('Sigmoid fit');
    
    label = dat(:, end); 
    miss_idx = label == 0;
    label = label(~miss_idx); label(label==-1) = 0;
    off = (dat(:,2) - median(dat(:, 2)))/numel(x_off);
    score = feval(f, off(~miss_idx));
    
    [X, Y, T, AUC, O] = perfcurve(label, score, 1);
    subplot(2, 2, 4);
    plot(X, Y, 'LineWidth', 1.5); hold on; 
    %J = Y - X; [~, idx] = max(J); p_thresh = T(idx);
    [p_thresh, idx] = minimJ(label, score);
    
    plot(X, X, 'k--', 'LineWidth', 1); hold on;
    scatter(O(1), O(2), 'r', 'filled');
    text(.5, .3, sprintf("AUC = %.2f", AUC), "FontSize", 9);
    text(.5, .2, sprintf("p(thresh) = %.2f", p_thresh), "FontSize", 9);
    xlabel('FPR'); ylabel('TPR');
    title('ROC, threshold');
    
    find_thresh = @(b, c, y) ((-(log(1 - y)-log(y))/b + c)...
                    *numel(x_off) + med);

    h = yline(ax, find_thresh(f.b, f.c, p_thresh),...
                          'k--', 't',...
                          'LabelHorizontalAlignment','left');

    fprintf("thresh: %.2f\tR^2: %.2f\n",...
                find_thresh(f.b, f.c, p_thresh),...
                gof.adjrsquare);

    % input(sprintf("%d",i));
    % close;
end

% t = 0:.05:1;
% for i = t
%     pred_lab = score > i;
%     C = confusionmat(label, pred_lab);
% end
fontsize(11, 'points');
clear miss_idx off;

%% change of threshold over time 
gap_block = 5;
start_block = 0; 
thresh_diff_time = zeros(40, numel(start_block));

k = 1;
for iter_block = start_block
    thresh_hap = calcThresholdmod(alldat,...
                                  1,...
                                  iter_block*18+1:...
                                  (iter_block+gap_block)*18);
    thresh_dis = calcThresholdmod(alldat,...
                                  2,...
                                  iter_block*18+1:...
                                  (iter_block+gap_block)*18);
    thresh_diff_time(:,k) = thresh_dis - thresh_hap; k = k+1;
end

figure; boxchart(thresh_diff_time, "Notch", "on");

clear thresh_hap thresh_dis;

%% something crazy-ish
% b in the sigmoid vs reaction latency
temp = 1; % how many emots to consider
rt = zeros(40,temp);
[thresh_, params, rsquare] = calcThresholdmod(alldat, 0, 1:360);
% params =zeros(40,temp);
% rsquare = zeros(40,temp);
% thresh = zeros(40,temp);

% for emot_idx = 1:3
%     [temp1, temp2, temp3] = calcThresholdmod(alldat,...
%                                           emot_idx,...
%                                           1:360);
%     thresh(:,emot_idx) = temp1;
%     params(:,emot_idx) = temp2(:,2);
%     rsquare(:,emot_idx) = temp3;
% 
%     for sub_idx = 1:40
%         if ~isnan(rsquare(sub_idx))
%             dat = alldat{sub_idx};
%             temp = dat(dat(:,1)==emot_idx, 4);
%             temp = temp(temp > .2 &...
%                         temp < prctile(temp, 97));
%             rt(sub_idx, emot_idx) = mean(temp);
%         end
%     end
% end

for sub_idx = 1:40
    dat = alldat{sub_idx};
    rt(sub_idx) = median(dat(...
                        dat(:,4) > .2 &...
                        dat(:,4) < prctile(dat(:,4), 97),...
                        4));
end

figure;
temp1 = rt(~isnan(rsquare));
temp2 = params(~isnan(rsquare),1);
scatter(log(temp2(temp2>0)), temp1(temp2>0),...
                    'k', 'filled',...
                    'MarkerEdgeColor', 'w');
hold on;
[p, ~] = quantreg(log(temp2(temp2>0)), temp1(temp2>0), .5);
plot(log(temp2(temp2>0)), polyval(p, log(temp2(temp2>0))), 'k-');
ylabel('Mean RT'); 
xlabel(['$log(b) ; p_{accept} = ' ...
        '\frac{1}{1 + e^{-b \cdot (O - c)}}$'],...
        'Interpreter', 'latex');


%[r, p] = corr(log(params(slice_idx,1)), rt(slice_idx));

%% plot for the old thresholding method
sub_idx = 3;
figure;

% Define helper funcitons to normalize from axis coordinates to normalized position in figure.
ax = gca();
xnorm = @(x)((x-ax.XLim(1))./(ax.XLim(2)-ax.XLim(1)))...
                .*ax.InnerPosition(3)+ax.InnerPosition(1);
ynorm = @(y)((y-ax.YLim(1))./(ax.YLim(2)-ax.YLim(1)))...
                .*ax.InnerPosition(4)+ax.InnerPosition(2);

dat = alldat{sub_idx};
resp = dat(dat(:,1)==2,end);
off = dat(dat(:,1)==2,2);
bar_with_hist(resp, off);
hold on;
annotation('textarrow', xnorm([-.2 -.3]),...
    ynorm([prctile(off(resp==1),20)+1.5 prctile(off(resp==1),20)]),...
    'Color', 'b', 'String', '20%ile(A)');
annotation('textarrow', xnorm([-.2 -.3]),...
    ynorm([prctile(off(resp==-1),80)-1.5 prctile(off(resp==-1),80)]),...
    'Color', 'r', 'String', '80%ile(R)');

temp = (prctile(off(resp==1),20) + prctile(off(resp==-1),80))/2;
yline(temp,...
      'k--', 'Thresh',...
      'LineWidth', 1.5,...
      'LabelHorizontalAlignment','left');

% create custom legend
temp = {};
temp{1} = histogram(1, 'FaceColor', 'b', 'Visible', 'off');
temp{2} = histogram(1, 'FaceColor', 'r', 'Visible', 'off');
legend([temp{:}], {"Accept" "Reject"});

fontsize(11, 'points');

%% explain why thresholds are different, but accpts ain't
figure;
temp_start = 1;
temp_gap = 10;
k = 1;
%temp_start:temp_gap+temp_start-1 
for sub_idx = [2 14 7 5 10 25 1 4]
    dat = alldat{sub_idx};
    for emot_idx = 1:2
        temp_dat = dat(dat(:,1)==emot_idx, :);
        p_accpt = bar_with_hist(temp_dat(:,end),...
                                temp_dat(:,2), 0);
        
        subplot(temp_gap, 2, k); k = k+1;
        plot(15:44, p_accpt, 'k.-'); hold on;
        xline(thresh_mat(sub_idx, emot_idx), 'r-');
        title(sum(temp_dat(:,end)==1)/size(temp_dat,1));
    end
end
% okie, 17 looks like a good representative

%% now actually explain -- it's about variance babe
figure;

sub_idx = 16;
dat = alldat{sub_idx};
accpt_rate = 1:2;
col = ["b" "r"];

for emot_idx = 1:2
    temp_dat = dat(dat(:,1)==emot_idx,:);
    p_accpt = bar_with_hist(temp_dat(:,end),...
                            temp_dat(:,2), 0);
    plot(15:44, p_accpt, '.-',...
        'LineWidth', 2,...
        'Color', char(col(emot_idx))); hold on;
    histogram('BinEdges', 14.5:44.5,...
              'EdgeColor', 'none',...
              'BinCounts', p_accpt,...
              'FaceColor', char(col(emot_idx)),...
              'FaceAlpha', .1);

    accpt_rate(emot_idx) = sum(temp_dat(:,end)==1)/size(temp_dat,1);
    % xline(thresh_old(sub_idx, emot_idx), '-',...
    %       sprintf("%.1f", thresh_old(sub_idx, emot_idx)),...
    %       'Color', char(col(emot_idx)),...
    %       'LineWidth', 2);
end
xlabel("Offer"); ylabel("p(Acceptance|Offer)");
xlim([15 44]);
% create custom legend
temp = {};
temp{1} = histogram(1, 'FaceColor', 'b', 'Visible', 'off');
temp{2} = histogram(1, 'FaceColor', 'r', 'Visible', 'off');
legend([temp{:}], {...
        sprintf("%s", "Happy"),...
        sprintf("%s", "Disgusted")});
fontsize(12, 'points');

%% assume differential degree of noise in Acceptance
% take different percentile ranges and cram
sub_idx = setdiff(1:40, [14 15 21 29 37:39]);
error_prctile = 1:100;
recons_data = zeros(numel(error_prctile), numel(sub_idx), 3);
hap_dis_diff = 1:100;

for iter = error_prctile
    thresh_mat = zeros(numel(sub_idx), 3);
    for emot_idx = 1:3
        thresh_mat(:,emot_idx) = ...
            calcThreshold(alldat(sub_idx),emot_idx,1:360,iter);
    end
    % hap_dis_diff(int16(iter/error_prctile(1))) =...
    %     mean(diff(thresh_mat, [], 2));
    recons_data(iter,:,:) = thresh_mat;
end

% figure; plot(error_prctile, hap_dis_diff, '.--');

%%
figure;
% temp_dat = zeros(100, 2);
for emot_idx = 1:3
    temp = mean(recons_data(:,:,emot_idx),2);
    % temp_dat(:,emot_idx) = temp;
    p = polyfit(temp, error_prctile, 2);
    plot((temp-15)./30,...
         error_prctile./100, '.-',...
         'LineWidth', 2.5); hold on;
    %plot(temp, polyval(p, temp), '--'); hold on;
end

% patch([(temp_dat(:,1)'-15)./30 fliplr((temp_dat(:,2)'-15)./30)],...
%         [error_prctile./100 fliplr(error_prctile./100)], 'k',...
%         'FaceAlpha', .1);

legend(["CDF_1: Happy" "CDF_2: Disgusted" "CDF_3: Neutral"]);
ylabel('Percentile value');
xlabel('Offer (scaled)');
fontsize(13, 'points');

%% 
n_sub = size(recons_data, 2);
temp = 1:n_sub; 

k = 1;

for sub_idx = 1:n_sub
    temp_2 = squeeze(recons_data(:,sub_idx,2));
    temp_1 = squeeze(recons_data(:,sub_idx,1));
    
    temp(sub_idx) = sum(temp_2 - temp_1);
end

        

        












