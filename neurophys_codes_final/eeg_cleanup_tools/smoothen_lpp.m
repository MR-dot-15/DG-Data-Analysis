function new_dat = smoothen_lpp(dat, electrode, channels, notch, smooth, show_plot)

fig = figure;
sgtitle(channels(electrode));

% extract relevant electrode details
new_dat = dat;
n_samp = size(dat, 2);
dat = squeeze(dat(electrode, :, :));

% plot the ERP
subplot(1,3,1);
plot(mean(dat, 2));
xlim([0 n_samp]); ylim([-20 20]);
xline(102, 'r--');
yline(0);

if notch
    % smooth the later part
    start_smoothing_after = round(n_samp * 550 / 1200);
    lfreq = 20;
    mov_win = round(512 / lfreq);

    untouched_part = dat(1:start_smoothing_after-1, :);
    smooth_part = dat(start_smoothing_after:end, :);
    smooth_part = movmean(smooth_part, mov_win, 1);
    dat = vertcat(untouched_part, smooth_part);
    
    % use notch-like removal
    % start_smoothing_after = 220;
    % fs = 512;
    % notch_f = 13;
    % bw = 5;
    % [b, a] = iirnotch(notch_f/(fs/2), bw/(fs/2));
    % 
    % untouched_part = dat(1:start_smoothing_after-1, :);
    % smooth_part = dat(start_smoothing_after:end, :);
    % smooth_part = filter(b, a, smooth_part);
    % dat = vertcat(untouched_part, smooth_part);
    
    % plot the smoothed ERP
    subplot(1,3,2);
    plot(mean(dat, 2));
    xlim([0 n_samp]); ylim([-35 35]);
    xline(102, 'r--');
    yline(0);
end
if smooth == 1
    % throughout smoothing
    lfreq = 45;
    mov_win = round(512 / lfreq);
    dat = movmean(dat, mov_win, 1);

    % plot
    subplot(1,3,3);
    plot(mean(dat, 2));
    xlim([0 n_samp]); ylim([-35 35]);
    xline(102, 'r--'); yline(0);
end

% update the data
new_dat(electrode, :, :) = dat;

if ~show_plot
    close(fig);
end