function all_dat = AR_correction(dat, electrode, p)
% AR-based noise removal

all_dat = dat;
dat = squeeze(dat(electrode, :, :));
dat_og = dat;

for idx = 1:size(dat, 2)
    col = dat(:, idx);
    [ARcoeff, ~] = aryule(col, p);
    estim_noise = filter(-ARcoeff, 1, col);
    
    cleaned_dat = col - estim_noise;

    dat(:, idx) = cleaned_dat;
end

figure;
subplot(1,2,1);
plot(mean(dat_og, 2));
xline(102, 'r--'); yline(0);
xlim([0 614]); ylim([-10 10]);

subplot(1,2,2);
plot(mean(dat, 2));
xline(102, 'r--'); yline(0);
xlim([0 614]); ylim([-10 10]);

all_dat(electrode, :, :) = dat;