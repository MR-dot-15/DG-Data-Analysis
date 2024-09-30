function new_dat = remove_sin(dat, elec, chan_label,...
            preserve_portion, smooth, show_fig)

new_dat = dat;
n_samp = size(dat, 2);
dat = squeeze(dat(elec, :, :));
ylim_arr = [-7 7];

if smooth == true
    ax = [0 0 0 0];
    n_fig = 4;
else
    ax = [0 0 0];
    n_fig = 3;
end

f = figure;
ax(1) = subplot(1,n_fig,1);
plot(mean(dat, 2));
xlim([0 n_samp]); ylim(ylim_arr);
xline(102, 'r--'); yline(0);

% piece-wisree linear detrending
% get there noise, scale it a bit
% remove the noise from the data
amp = .6; %.2 .* min(mean(dat, 2));
bp = 1:15:n_samp;
noise = amp.*detrend(dat, 1, bp, "Continuous", true);

% preserve the initial peaks
if preserve_portion 
    time_pts = 120:n_samp;
    noise(time_pts, :) = zeros(length(time_pts),...
                                size(noise, 2));
end

ax(2) = subplot(1,n_fig,2);
plot(mean(noise, 2));
xlim([0 n_samp]); ylim(ylim_arr);
xline(102, 'r--'); yline(0);

ax(3) = subplot(1,n_fig,3);
dat = dat - noise;
plot(mean(dat, 2));
xlim([0 n_samp]); ylim(ylim_arr);
xline(102, 'r--'); yline(0);

if smooth == true
    ax(4) = subplot(1,n_fig,4);

    lfreq = 55;
    mov_win = round(512 / lfreq);
    dat = movmean(dat, mov_win, 1);

    plot(mean(dat, 2));
    xlim([0 n_samp]); ylim(ylim_arr);
    xline(102, 'r--'); yline(0);
end

new_dat(elec, :, :) = dat;

if show_fig == true
    linkaxes(ax, 'xy');
    sgtitle(chan_label(elec));
else
    close(f);
end