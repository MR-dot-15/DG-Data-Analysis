function ind = detect_outlier(data)
% DETECT_OUTLIER(data) returns the indices of outlier vals
q1 = quantile(data, .25);
q3 = quantile(data, .75);
iqr = q3 - q1;

llim = q1 - 1.5*iqr;
ulim = q3 + 1.5*iqr;

ind = or(data<llim, data>ulim);
fprintf("%.2f points are outlier\nOutliers: ", sum(ind)/length(ind));
disp(data(ind));