function [threshold, params, rsquare] = calcThresholdmod(alldat, emot, llim_ulim)

% define the sigmoid
logistic = fittype('1/(1+exp(-b*(x-c)))',...
            'dependent',{'y'},'independent',{'x'},...
            'coefficients',{'b','c'});
% define some parameters/storage
off_uniq = 15:44; med = median(off_uniq);
threshold = 1:length(alldat);
params = zeros(length(alldat),2);
rsquare = 1:length(alldat);
permitted_error = .1; % rej/accpt by mistake

for sub_idx = 1:length(alldat)

    % slice data
    dat = alldat{sub_idx};

    % find all acceptors/ rejectors
    rej_rate = sum(dat(:,end)==-1)/size(dat,1);

    if rej_rate > permitted_error && rej_rate < 1 - permitted_error
        % centralise and scale data
        off_uniq_center = (off_uniq - med)/...
                        (max(off_uniq)-min(off_uniq));
        % emot-specific thresh
        if emot == 0
            dat = dat(llim_ulim, :);
        else
            dat = dat(llim_ulim, :);
            dat = dat(dat(:,1)==emot, :);
        end
        
        % step-1 calculate the PMF of acceptance
        label = dat(:, end); missed_idx = label == 0;
        label = label(~missed_idx); label(label == -1) = 0;
        off = dat(:,2); off = off(~missed_idx);
    
        % if acceptance, rejection both
        if numel(unique(label)) ~= 1
            [c_off, ~] = hist(off, off_uniq);
            accpt_off = off(label == 1);
            [c_accpt, ~] = hist(accpt_off, off_uniq);
            p_accpt = c_accpt ./ c_off;
            
            % if n_block < 20
            % some offers might not have appeared - NaN acceptance
            if sum(isnan(p_accpt))>0
                rem_idx = isnan(p_accpt);
                p_accpt = p_accpt(~rem_idx);
                off_uniq_center = off_uniq_center(~rem_idx);
            end
            
            % step-2 fit the sigmoid to PMF
            [f, gof] = fit(off_uniq_center', p_accpt', logistic,...
                'StartPoint', [10 0],...
                'MaxIter', 800);
            score = feval(f, (off - med)./...
                        (max(off) - min(off)));

            % store adj R^2 values
            rsquare(sub_idx) = gof.adjrsquare;

            % if the fit is extremely poor
            if gof.adjrsquare < .2
                threshold(sub_idx) = calcThreshold(alldat(sub_idx),...
                                            emot, llim_ulim, 20);
            else
                % store fitted param vals
                params(sub_idx,:) = [f.b f.c];
                
                % step-3 calculate ROC & optimal operating point
                % by maximizing the J-statistic
                p_thresh = minimJ(label, score);
            
                off_thresh = -(log(1-p_thresh) - log(p_thresh))/f.b + f.c;
                threshold(sub_idx) = off_thresh * (max(off_uniq)...
                                                  -min(off_uniq))...
                                                  + med;
            end

        % if everything's been accepted thus far
        elseif unique(label) == 1
            threshold(sub_idx) = min(off);
    
        % if everything's been rejected thus far
        elseif unique(label) == 0
            threshold(sub_idx) = max(off);
        end

    % for all acceptors
    elseif rej_rate < permitted_error
        threshold(sub_idx) = min(dat(:,2));
        rsquare(sub_idx) = nan;

    % for all rejectors
    elseif rej_rate > permitted_error
        threshold(sub_idx) = max(dat(:,2));
        rsquare(sub_idx) = nan;
    end
end