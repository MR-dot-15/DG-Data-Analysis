function new_dat = trend_correction(dat, elec, chan_label, smooth, show_plot)

new_dat = dat;
n_samp = size(dat, 2);
dat = squeeze(dat(elec, :, :));
f1 = figure;
ylim_arr = [-20 20];

ax1 = subplot(1,3,1);
plot(mean(dat, 2));
xlim([0 n_samp]); ylim(ylim_arr);
xline(102, 'r--'); yline(0);
title(num2str(elec));

bp = [];
dat = detrend(dat, 1, bp, "Continuous", true);
% T = getTrend(dat);
% dat = detrend(dat, T);
dat = dat - mean(dat(1:102, :), 1);
ax2 = subplot(1,3,2);
plot(mean(dat, 2));
xlim([0 n_samp]); ylim(ylim_arr);
xline(102, 'r--'); yline(0);
title(strcat(num2str(elec), " detr"));

if smooth == true
    lfreq = 60;
    mov_win = round(512 / lfreq);
    dat = movmean(dat, mov_win, 1);

    ax3 = subplot(1,3,3);
    plot(mean(dat, 2));
    xlim([0 n_samp]); ylim(ylim_arr);
    xline(102, 'r--'); yline(0);
    title(strcat(num2str(elec), " detr, smth"));
end

new_dat(elec, :, :) = dat;


if ~show_plot 
    close(f1);
else
    linkaxes([ax1 ax2], 'xy');
    sgtitle(chan_label(elec));
end