function prob_accpt = bar_with_hist(accpt_val_logical, offer, plot_fig)

bins = 15:44;
[c_off, ~] = hist(offer, bins);
target_resp = [-1 1];
col = ['r'; 'b'];


%figure; 
% subplot(1, 2, 1);
hold on;
k = 1;
for i = 1
    temp = offer(accpt_val_logical==i);
    [cnt, ~] = hist(temp, bins);
    cnt_norm = cnt ./ c_off;
    if sum(isnan(cnt_norm))>0
        cnt_norm(isnan(cnt_norm)) = 0;
    end

    if plot_fig
        histogram('BinEdges', (bins(1)-.5):1:(bins(end)+.5),...
                  'BinCounts', cnt_norm .* .2,...
                  'FaceAlpha', .2,...
                  'FaceColor', col(k), 'LineWidth', .8,...
                  'EdgeColor', col(k), 'EdgeAlpha', .5,...
                  'Orientation', 'horizontal');
        boxchart(-.3*ones(size(offer(accpt_val_logical==i))),...
                 offer(accpt_val_logical==i),...
                 "Notch", "on", "BoxFaceColor", col(k), "BoxWidth", .1);
    end
    k = k+1;
end

if plot_fig
    xlim([-.5 .4]);
    xticklabels([]);
    ylabel('offer');
end

prob_accpt = cnt_norm;