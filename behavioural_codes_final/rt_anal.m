%% a flavour of bw sub variation
figure;
iter = 20;
[~, rt_pulled, pull_idx] = prep_rt_matrix(behavDataSummary(alldat,...
                                    1, 18*iter, 4));

tiledlayout(2, 3, "TileSpacing","compact");

% rt < certain hard threshold
below_thresh_count = 0;
threshold = 0.4;

for i = 1:2
    for j = 1:3
        nexttile;
        pull_idx_mod = pull_idx{i, j};
        %remove_AA_idx = ismember(pull_idx_mod, AAsubj);

        % remove all acceptors
        pull_idx_mod = pull_idx_mod(:);
        rt_pulled_mod = rt_pulled{i,j}(:);

        boxchart(pull_idx_mod, rt_pulled_mod, 'MarkerStyle', '.');
        title(sprintf("%s with %s",...
            emot_label(i), off_label(j)));
        xlim([0 length(unique(pull_idx{i,j})) + 1]); ylim([0 2]);
        xticks(1:length(sublist)); xticklabels(sublist); xtickangle(90);
        xlabel("subject ID"); ylabel("response latency (s)");
        yline(threshold);

        below_thresh_count = below_thresh_count +...
                             sum(rt_pulled_mod  <= threshold);
    end
end
sgtitle("subjectwise variability in rt");

clear i j iter;

%%
figure;

for off = [2 4]
    for emot = 1:2
        pull_idx_mod = pull_idx{emot,off-1}(:);
        
        % manage to put 2:8 and 4:6 together
        if off == 2
            pull_idx_mod = 3*pull_idx_mod-2;
        else
            pull_idx_mod = 3*pull_idx_mod-1;
        end

        rt_pulled_mod = rt_pulled{emot,off-1}(:);
        
        % colour code
        if emot == 1
            col = 'b';
        else
            col = 'r';
        end
        boxchart(pull_idx_mod, rt_pulled_mod,...
            'MarkerStyle', '.', 'BoxFaceColor', col, 'BoxFaceAlpha', 0.06);
        hold on;
    end
end

legend(["happy" "disgusted"]);

%% 1 14 21
sub_of_interest = [1 14 21];

tiledlayout(1, 3, "TileSpacing", "compact");
for sub = sub_of_interest
    nexttile;
    for off = [2 4]
        for emot = 1:2
            % colour code
            if emot == 1
                col = 'b';
            else
                col = 'r';
            end
     
            pull_idx_mod = pull_idx{emot,off-1} == sub;
            rt_pulled_mod = rt_pulled{emot,off-1}(pull_idx_mod);

            boxchart(off/2 .* ones(size(rt_pulled_mod)),...
                rt_pulled_mod,...
                'MarkerStyle', '.',...
                'BoxFaceColor', col,...
                'BoxFaceAlpha', 0.06);
            ylabel("rt (s)"); ylim([0 2]);
            xticks([]);
            text(.8, .2, "2:8");
            text(1.8, .2, "4:6");
            hold on;
        end
    end
    legend(["happy" "disgusted"]);
    title(sprintf("sub: %d", sub));
end


%% rt histograms: pulled, sliced by acceptance/rejection across cond

fig = figure;
tiledlayout(3, 3, "TileSpacing", "compact");

iter = 20;
threshold = .25; % max of min rt
AA_idx = [];
sublist_mod = sublist(~AA_idx);
pulled_rt_processed = behavDataSummary(alldat,...
                            1, 18*iter, 4);
pulled_resp_processed = behavDataSummary(alldat,...
                            1, 18*iter, 7);
% cond sliced rt for acc, rej
pulled_rt_cond_cell_a = cell(3, 3);
pulled_rt_cond_cell_r = cell(3, 3);

n_subj = length(pulled_rt_processed); nbin = 20;

for emot = 1:3
    for off = 2:4
        pulled_rt_cond = [];
        pulled_resp_cond = [];

        % iterate over subjects
        for subj = 1:n_subj
            index = 10*emot + off;
            pulled_rt_cond = vertcat(pulled_rt_cond,...
                pulled_rt_processed{subj}{index}');
            pulled_resp_cond = vertcat(pulled_resp_cond,...
                pulled_resp_processed{subj}{index}');
        end

        % KS test
        pulled_rt_cond_a = pulled_rt_cond(pulled_resp_cond == 1);
        pulled_rt_cond_a = pulled_rt_cond_a(pulled_rt_cond_a > threshold);

        pulled_rt_cond_r = pulled_rt_cond(pulled_resp_cond == -1);
        pulled_rt_cond_r = pulled_rt_cond_r(pulled_rt_cond_r > threshold);

        [~, p] = kstest2(pulled_rt_cond_a, pulled_rt_cond_r);

        % logn fit
        % p_a = lognfit(pulled_rt_cond_a);
        % %disp(lognstat(p_a(1), p_a(2)));
        % disp(mean(pulled_rt_cond_a));
        % 
        % p_r = lognfit(pulled_rt_cond_r);
        % %disp(lognstat(p_r(1), p_r(2)));
        % disp(mean(pulled_rt_cond_r));
        % 
        % disp('----');

        % save the sliced data
        pulled_rt_cond_cell_a{emot, off - 1} = pulled_rt_cond_a;
        pulled_rt_cond_cell_r{emot, off - 1} = pulled_rt_cond_r;
        
        % plot
        nexttile;
        [n1, m1] = hist(pulled_rt_cond_a, nbin);
        [n2, m2] = hist(pulled_rt_cond_r, nbin);
        bar(m1, n1, 'b',...
            'FaceAlpha', 0.6,...
            'EdgeColor','none',...
            'BarWidth', 1); hold on; 
        bar(m2, -n2, 'r',...
            'FaceAlpha', 0.6,...
            'EdgeColor','none',...
            'BarWidth', 1); hold off;
        xlim([0 2]); ylim([-170 230]);
        yticks([-150 0 200]); yticklabels(["150", "0", "200"]);
        % title(sprintf("%s with %s",...
        %     emot_label(emot), off_label(off-1)));
        %legend(["Accepted" "Rejected"], 'FontSize', 15);
        %legend boxoff;
        set(gca, 'FontSize', 15);
        set(gca, 'box', 'off');
    end
end
% sgtitle("Response Latency Across Conditions",...
%         'FontWeight', 'bold',...
%         'FontSize', 18);
ax = axes(fig,'visible','off'); 
ax.Title.Visible ='on';
ax.XLabel.Visible ='on';
ax.YLabel.Visible ='on';
ylabel(ax,'Frequency', 'FontSize', 20);
xlabel(ax,'Reaction time (s)', 'FontSize', 20);
% title(ax,"Response Latency Across Conditions",...
%          'FontWeight', 'bold',...
%          'FontSize', 20);

clear n_subj emot off index subj n1 n2 m1 m2 pulled_rt_cond...
    pulled_resp_cond p pulled_rt_cond_a pulled_rt_cond_r;

%% fit the data to lognorm/gamma
threshold = .25;

rt_clean_accpt = pulled_rt(pulled_resp == 1);
rt_clean_accpt = rt_clean_accpt(rt_clean_accpt > threshold);

rt_clean_rej = pulled_rt(pulled_resp == -1);
rt_clean_rej = rt_clean_rej(rt_clean_rej > threshold);

nbin = 30;
figure;

% ------------------- acceptance -------------------
subplot(2, 4, 1); histogram(rt_clean_accpt, nbin);
title("acceptance rt cleaned");
subplot(2, 4, 2); histogram(log(rt_clean_accpt), nbin);
title("log transformed");

% lognormal fit
p = lognfit(rt_clean_accpt);
subplot(2, 4, 3); 
histogram(rt_clean_accpt, nbin,...
        'Normalization', 'pdf',...
        'DisplayStyle', 'stairs'); hold on; 
xarray = threshold:.05:2;
plot(xarray, lognpdf(xarray, p(1), p(2))); hold off;
[m, v] = lognstat(p(1), p(2));
fprintf("mu: %.2f sig: %.2f\n", m, v);
title("log-normal fit");

% gamma fit
p = gamfit(rt_clean_accpt);
subplot(2, 4, 4); 
histogram(rt_clean_accpt, nbin,...
        'Normalization', 'pdf',...
        'DisplayStyle', 'stairs'); hold on; 
xarray = threshold:.05:2;
plot(xarray, gampdf(xarray, p(1), p(2))); hold off;
title("gamma fit");

% ------------------- rejection -------------------
subplot(2, 4, 5); histogram(rt_clean_rej, nbin);
title("rejection rt cleaned");
subplot(2, 4, 6); histogram(log(rt_clean_rej), nbin);
title("log transformed");

% lognormal fit
p = lognfit(rt_clean_rej);
subplot(2, 4, 7); 
histogram(rt_clean_rej, nbin,...
        'Normalization', 'pdf',...
        'DisplayStyle', 'stairs'); hold on; 
xarray = threshold:.05:2;
plot(xarray, lognpdf(xarray, p(1), p(2))); hold off;
[m, v] = lognstat(p(1), p(2));
fprintf("mu: %.2f sig: %.2f\n", m, v);
title("log-normal fit");

% gamma fit
p = gamfit(rt_clean_rej);
subplot(2, 4, 8); 
histogram(rt_clean_rej, nbin,...
        'Normalization', 'pdf',...
        'DisplayStyle', 'stairs'); hold on; 
xarray = threshold:.05:2;
plot(xarray, gampdf(xarray, p(1), p(2))); hold off;
title("gamma fit");

clear p m v;

%% rt-trend as per face
figure; 

for emot_idx = 1:3
    rt_pulled = [];
    face_idx_pulled = [];
    for face_idx = 1:8
        for accpt = [1 -1]
            for sub_idx = 1:40
                dat = alldat{sub_idx};
                face_spec_idx = (dat(1:180,1) == emot_idx &...
                     findBaseOff(dat(1:180,2)) == 2 &...
                                 dat(1:180,7) == accpt &...
                                 dat(1:180,5) == face_idx);
                n_face_spec_idx = sum(face_spec_idx);
                rt_pulled = [rt_pulled...
                    dat(face_spec_idx, 4)'];

                if accpt == 1
                    face_idx_pulled = [face_idx_pulled...
                        face_idx * ones(1, n_face_spec_idx)];
                else
                    face_idx_pulled = [face_idx_pulled...
                        face_idx * ones(1, n_face_spec_idx) + .5];
                end
            end
        end
        subplot(3, 1, emot_idx); 
        title(sprintf("%d emot | 2 off", emot_idx));
        boxchart(face_idx_pulled, log(rt_pulled));
        ylim([-2 log(2)]);
    end
end
%% two-way anova
p_val_array = zeros(40, 3);

for sub_idx = 1:40
    dat = alldat{sub_idx};
    slice_idx = dat(:, 4) > .25 &...
                dat(:, 7) == -1;
    rt_pulled = log(dat(slice_idx, 4));
    emot = dat(slice_idx, 1)';
    off = findBaseOff(dat(slice_idx, 2))';
    p_val_array(sub_idx,:) =...
        anovan(rt_pulled, {emot off},...
        'display', 'off',...
        'model', 'interaction')';
end

fprintf("with effect of emot:\n___________________________\n");
for sub_idx = 1:40
    if p_val_array(sub_idx, 1) < .05
        fprintf("sub: %s\tp: %.3f\n",...
        sublist(sub_idx), p_val_array(sub_idx, 1));
    end
end
fprintf("with effect of off:\n___________________________\n");
for sub_idx = 1:40
    if p_val_array(sub_idx, 2) < .05
        fprintf("sub: %s\tp: %.3f\n",...
        sublist(sub_idx), p_val_array(sub_idx, 2));
    end
end
fprintf("with effect of emot*off:\n___________________________\n");
for sub_idx = 1:40
    if p_val_array(sub_idx, 3) < .05
        fprintf("sub: %s\tp: %.3f\n",...
        sublist(sub_idx), p_val_array(sub_idx, 3));
    end
end
%% IT'S RM-ANOVA TIME BABYYY LESSSGOOOO!
% details:
% in the design, we have a 40 * 270 table
% for each sub, for each cond, we take first 30 rts
% first 30, cuz a portion of the data is unbalanced
% within design -- emot and off labels
% y1 + y2 + ... + y270 ~ 1, grouped by emot + off + emot*off

% cond-specific rt-s
[rt_behavDat, ~, blk_idx] = behavDataSummary(alldat, 1, 360, 4);
[accpt_behavDat, ~, ~] = behavDataSummary(alldat, 1, 360, 7);

% first construct the repeat-mes data matrix and within factors
n_rt_pts = 15;
emot_label = [];
off_label = [];
blk_label = [];
rm_rt_mat = zeros(40, 9*n_rt_pts); % n_sub = 40, we take 30 trials for each 9 cond

for sub_idx = 1:40
    dat = alldat{sub_idx};
    accpt_rate = sum(dat(:,7)==1)/360;
    
    temp_rt_holder = [];
    temp_blk_idx_holder = [];
    for emot_idx = 1:3
        for off_idx = 2:4
            key = 10*emot_idx + off_idx;

            % slice rt data
            temp_rt_array = rt_behavDat{sub_idx}{key};
            temp_accpt_array = accpt_behavDat{sub_idx}{key};

            okay_rt = temp_rt_array > .25 & temp_rt_array ~= 2;

            temp_rt_array = temp_rt_array(okay_rt);
            try
                temp_rt_array = temp_rt_array(1:n_rt_pts);
            catch
                temp_rt_array = NaN(1, n_rt_pts);
            end

            % store the rt & blk idx data
            temp_rt_holder = [temp_rt_holder...
                log(temp_rt_array)];
            
            if sub_idx == 3
                % store off, emot label
                emot_label = [emot_label...
                    emot_idx * ones(1, n_rt_pts)];
                off_label = [off_label...
                    off_idx * ones(1, n_rt_pts)];

                % store block label
                % slice block index
                temp_blk_idx = blk_idx{sub_idx}{key};
                temp_blk_idx = temp_blk_idx(okay_rt);

                try
                    temp_blk_idx = temp_blk_idx(1:n_rt_pts);
                catch 
                    temp_blk_idx = NaN(n_rt_pts, 1);
                end

                temp_blk_idx_holder = [temp_blk_idx_holder...
                temp_blk_idx'];
            end
        end
    end

    rm_rt_mat(sub_idx, :) = temp_rt_holder;
    if sub_idx == 3
        blk_label = [blk_label temp_blk_idx_holder];
    end
end

% construct the table
tbl_sub = table((1:40)', 'VariableNames', "sub");
tbl_sub.sub = categorical(tbl_sub.sub);
y = rm_rt_mat;
tbl_dat = array2table(y);
tbl = [tbl_sub tbl_dat];

within_des = table(emot_label', off_label', blk_label',...
    'VariableNames', ["emot", "off", "blk"]);
within_des.emot = categorical(within_des.emot);
within_des.off = categorical(within_des.off);
within_des.blk = categorical(within_des.blk);

% rm anova
variable_str = sprintf("y1-y%d ~ 1",...
                            size(y,2));
rm = fitrm(tbl, variable_str, WithinDesign = within_des);
disp(ranova(rm, 'WithinModel', 'blk+emot+off+blk*emot+blk*off+emot*off'));


clear temp* tbl_dat tbl_sub y variable_str;
%% another design for rm anova -- now with "learning effect"

% calculate sub-wise accept rate, across conds
AA_idx = []; n_subj = 40 - numel(AA_idx);
[condSpecrt, ~, blk_idx] = behavDataSummary(alldat,1,360,4);
[condSpecaccpt, ~, ~] = behavDataSummary(alldat,1,360,7);

% first construct the repeat-mes data matrix and within factors
cluster_n_blocks = 5;
emot_label = [];
off_label = [];
blk_label = [];
rm_rt_mat = zeros(40, 9*(20/cluster_n_blocks));

for sub_idx = 1:40
    
    temp_rt_holder = [];

    for emot_idx = 1:3
        for off_idx = 2:4
            key = 10*emot_idx + off_idx;

            % slice accpt data
            temp_rt_array = log(condSpecrt{sub_idx}{key});
            temp_accpt_array = condSpecaccpt{sub_idx}{key};
            temp_blk_idx = blk_idx{sub_idx}{key};

            % store the rt, emot, off & blk idx data
            blk_iter = 1;
            for start_block = 1:cluster_n_blocks:...
                                (21-cluster_n_blocks)
                slice_idx = temp_blk_idx >= start_block &...
                            temp_blk_idx < start_block +...
                            cluster_n_blocks;
                temp_rt_array_ = temp_rt_array(...
                            slice_idx' &...
                            temp_rt_array < prctile(temp_rt_array, 99));

                rt_val = mean(temp_rt_array_(...
                               temp_rt_array_ > log(.2)));

                temp_rt_holder = cat(2, temp_rt_holder,...
                                            rt_val);

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

    rm_rt_mat(sub_idx, :) = temp_rt_holder;
end

% construct the table
tbl_sub = table((1:n_subj)', 'VariableNames', "sub");
tbl_sub.sub = categorical(tbl_sub.sub);
y = rm_rt_mat;
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

clear temp_rt_array_ temp_rt_array rt_val...
        tbl_sub tbl_dat emot_idx off_idx y;

%% generate the anova figure;
figure;

% emot
subplot(2, 3, 1); 
[temp_mu, temp_se] = rm_anova_mean_se(tbl,...
                    double(within_des.emot));
errorbar(1:3, temp_mu, temp_se, 'k-',... 
         'LineWidth', 2, 'CapSize', 8,...
         'Marker', 'o',...
         'MarkerFaceColor', [.8 .8 .8],...
         'MarkerEdgeColor', 'w');
xticks(1:4);
xticklabels(["Hap" "Disg" "Neu"])
xlim([0 4]);
ylim([-.37 -.31]);
ylabel('log(rt)');
set(gca, 'FontSize', 15);
title('Emotion', 'FontSize', 18); 
set(gca, 'box', 'off');

% off
subplot(2, 3, 4);
[temp_mu, temp_se] = rm_anova_mean_se(tbl,...
                    double(within_des.off));
errorbar(1:3, temp_mu, temp_se, 'k-',...
         'LineWidth', 2, 'CapSize', 8,... 
         'Marker', 'o',...
         'MarkerFaceColor', [.8 .8 .8],...
         'MarkerEdgeColor', 'w');
xticklabels(["Low" "Interm" "Max"])
xlim([0 4]); ylim([-.4 -.25]);
ylabel('log(rt)');
set(gca, 'FontSize', 15);
title('Offer', 'FontSize', 18); 
set(gca, 'box', 'off');

% emot * off
subplot(2, 3, 2);
[temp_mu, temp_se] = rm_anova_mean_se(tbl,...
                    double(within_des.e_o));

for idx = 1:3:9
    slice_idx = idx:(idx+2);
    errorbar(slice_idx, temp_mu(slice_idx),...
             temp_se(slice_idx), '-',...
             'LineWidth', 2, 'CapSize', 8,... 
             'Marker', 'o',...
             'MarkerFaceColor', [.8 .8 .8],...
             'MarkerEdgeColor', 'w');
    hold on;
end
xticks(1:9);
xticklabels(["H/L" "H/I" "H/M"...
             "D/L" "D/I" "D/M"...
             "N/L" "N/I" "N/M"]);
xlim([0 10]); ylim([-.44 -.24]);
legend(["Hap", "Dis", "Neu"]);
legend boxoff;
ylabel('log(rt)');
set(gca, 'FontSize', 15);
title('Emotion x Offer', 'FontSize', 18); 
set(gca, 'box', 'off');

% blk
subplot(2, 3, 5);
[temp_mu, temp_se] = rm_anova_mean_se(tbl,...
                    double(within_des.blk));
errorbar(1:4, temp_mu, temp_se, 'k-',... 
         'LineWidth', 2, 'CapSize', 8,...
         'Marker', 'o',...
         'MarkerFaceColor', [.8 .8 .8],...
         'MarkerEdgeColor', 'w');
xticklabels(["1-5", "6-10", "11-15", "16-20"]);
xlim([0 5]); %ylim([.56 .72]);
ylabel("log(rt)");
set(gca, 'FontSize', 15);
title("Block-quarter", 'FontSize', 18); 
set(gca, 'box', 'off');

% blk * emot
subplot(2, 3, 3);
[temp_mu, temp_se] = rm_anova_mean_se(tbl,...
                    double(within_des.e_b));
for idx = 1:4:12
    slice_idx = idx:(idx+3);
    errorbar(slice_idx, temp_mu(slice_idx),...
             temp_se(slice_idx), '-',...
             'LineWidth', 2, 'CapSize', 8,... 
             'Marker', 'o',...
             'MarkerFaceColor', [.8 .8 .8],...
             'MarkerEdgeColor', 'w');
    hold on;
end
xticks(1:12);
xticklabels(["H1" "H2" "H3" "H4"...
             "D1" "D2" "D3" "D4"...
             "N1" "N2" "N3" "N4"]);
xlim([0 13]); 
ylim([-.47 -.18]);
legend(["Hap", "Dis", "Neu"]); 
legend boxoff;
ylabel('log(rt)');
set(gca, 'FontSize', 15);
title('Emotion x Block', 'FontSize', 18);
set(gca, 'box', 'off');

% blk * off
subplot(2, 3, 6);
[temp_mu, temp_se] = rm_anova_mean_se(tbl,...
                    double(within_des.o_b));
for idx = 1:4:12
    slice_idx = idx:(idx+3);
    errorbar(slice_idx, temp_mu(slice_idx),...
             temp_se(slice_idx), '-',...
             'LineWidth', 2, 'CapSize', 8,... 
             'Marker', 'o',...
             'MarkerFaceColor', [.8 .8 .8],...
             'MarkerEdgeColor', 'w');
    hold on;
end
xticks(1:12);
xticklabels(["L1" "L2" "L3" "L4"...
             "I1" "I2" "I3" "I4"...
             "M1" "M2" "M3" "M4"]);
legend(["2:8", "3:7", "4:6"]);
legend boxoff;
xlim([0 13]); 
ylim([-.5 -.17])
ylabel('log(rt)');
set(gca, 'FontSize', 15);
title('Offer x Block', 'FontSize', 18); 
set(gca, 'box', 'off');

%% rt vs (off, emot) scatter
% tbl needs the response glme model to be executed
figure; hold on;
for emot = 1:3
    off_uniq = unique(tbl.off(tbl.emot == emot));
    rt_uniq = zeros(size(off_uniq));

    for i = 1:length(off_uniq)
        rt_uniq(i) = mean(tbl.rt(...
                          tbl.emot == emot &...
                          tbl.off == off_uniq(i)));
    end
    plot(off_uniq, movmean(rt_uniq, 5), 'o-');
end

grid('on');
xline([25 35], 'k--', 'LineWidth', 2);
xlabel('offered money'); ylabel('log(rt)');
title('scatter: log(rt) against off, grouped by emot');
legend(["happy" "disgusted" "neutral"]);

%% find rt < .4 across subs (off-epoched ERP)
rt_index = cell(1, 40);

for sub_id = 1:40
    dat = alldat{sub_id};
    sub = sublist(sub_id);
    rt = dat(:, 4);

    if contains(sub, ["mis" "vis"])
        rt = rt(1:90);
    elseif contains(sub, ["gau" "swa"])
        rt = rt(1:180);
    elseif contains(sub, ["kum" "biv"])
        rt = rt(1:270);
    end

    rt_index{sub_id} = rt > .4;
end

%% quantile reg
figure; hold on;
dat = alldat{16};

scatter(dat(:,2), log(dat(:,4)), 'filled', 'MarkerFaceAlpha', .3);
emot = dat(:, 1); off = dat(:, 2); rt = log(dat(:, 4));

[p, stat] = quantreg(off(emot == 1), rt(emot == 1), .5);
disp(p);
plot(off(emot == 1), polyval(p, off(emot == 1)));
[p, stat] = quantreg(off(emot == 2), rt(emot == 2), .5);
disp(p);
plot(off(emot == 2), polyval(p, off(emot == 2)));
[p, stat] = quantreg(off(emot == 3), rt(emot == 3), .5);
disp(p);
plot(off(emot == 3), polyval(p, off(emot == 3)));

legend(["dat" "hap" "dis" "neu"]);

%% quantile reg param : emot specific
n_sub = 40;
quant = .5; % median
param = cell(2, 3); % to store the slope & intercept
for sub_idx = 1:6
    param{sub_idx} = [];
end

for sub_idx = 1:n_sub
    dat = alldat{sub_idx};
    
    for emot_idx = 1:3
        off_slice = dat(dat(:,1)  == emot_idx, 2);
        off_slice = off_slice ./ max(off_slice);
        rt_slice = log(dat(dat(:,1) == emot_idx, 4));

        [p, ~] = quantreg(off_slice, rt_slice, quant);
        param{1, emot_idx} = horzcat(param{1, emot_idx},...
                                        p(1)); % slope
        param{2, emot_idx} = horzcat(param{2, emot_idx},...
                                        p(2)); % intercept
    end
end

% KS-test performed (pair-wise)
% no result

%% param: accept-reject
param = cell(2, 2); % to store the slope & intercept
for sub_idx = 1:4
    param{sub_idx} = [];
end

for sub_idx = setdiff(1:n_sub, [21 29])
    dat = alldat{sub_idx};
    
    k = 1;
    for emot_idx = [1 -1]
        off_slice = dat(dat(:,7)  == emot_idx, 2);
        off_slice = off_slice ./ max(off_slice);
        rt_slice = log(dat(dat(:,7) == emot_idx, 4));

        [p, ~] = quantreg(off_slice, rt_slice, quant);
        param{1, k} = horzcat(param{1, k},...
                                        p(1)); % slope
        param{2, k} = horzcat(param{2, k},...
                                        p(2)); % intercept
        k = k+1;
    end
end

%% plot the param distribution : emot

figure;

lab_1 = ["slope" "interp"];
lab_2 = ["hap" "dis" "neu"];
for k = 1:2
    for emot_idx = 1:3
        subplot(2, 3, (k-1)*3 + emot_idx); hold on;

        % expected val
        xline(mean(param{k, emot_idx}), 'r--');
        % 0
        xline(0, 'k--');
        % hist
        histogram(param{k, emot_idx}, -1.5:.2:1);

        % legend
        legend("\mu");

        % accessories
        xlim([-1.5 1.5]); ylim([0 20]);
        title(lab_1(k) + " | " + lab_2(emot_idx));
        grid('on')
    end
end

%% plot param distribution : A/R
figure;


lab_1 = ["slope" "interp"];
lab_2 = ["acc" "rej"];
for k = 1:2
    for emot_idx = 1:2
        subplot(2, 2, (k-1)*2 + emot_idx); hold on;

        % expected val
        xline(mean(param{k, emot_idx}), 'r--');
        % 0
        xline(0, 'k--');
        % hist
        histogram(param{k, emot_idx}, -1.5:.2:1.5);

        % legend
        legend("\mu");

        % accessories
        ylim([0 20]); xlim([-1.5 1.5]);
        title(lab_1(k) + " | " + lab_2(emot_idx));
        grid('on')
    end
end

% kstest both significant

%% try eliminating neutral ka effect















