function [earned, max_earned, all_block_earning] = blockEarn(dat)
% BLOCKEARN(dat) returns the earning/learning pattern over one experiment
%
% Input:
% dat: behavData matrix
%
% Output:
% earned: earned in each block (dim- 1,20)
% max_earned: maximum possible earning in each block
% all_block_earning: earning in each trial in each block (dim- 20, 18)
% threshold: currently unavailable

% setting up variables
blockN = 20;
trialN = 18;
earned = 1:blockN;
all_block_earning = zeros(blockN, trialN);
max_earned = 1:blockN;
% threshold = 1:blockN;

% 7th col is resp
% note: R -1; A 1; miss 0
resp = dat(:,7);
resp_processed = (resp + 1)/2;
% 2nd col is offers
money = dat(:,2);
money_exceptMissed = money;
% set offer values = 0 when it's a miss
money_exceptMissed(resp == 0) = 0;

% if everything were accepted
resp_allAccepted = abs(resp);

% trial wise earning
perTrialEarn = money_exceptMissed...
    .* resp_processed;
perTrialMax = money_exceptMissed...
    .* resp_allAccepted;


for i = 1:blockN
    blkWiseSection = perTrialEarn( ...
        (i-1)*trialN+1:trialN*i);
    %disp(blkWiseSection(blkWiseSection>0));
    % store each trial earning into all_block_earning
    all_block_earning(i,:) = blkWiseSection;
    blkMaxEarned = perTrialMax( ...
        (i-1)*trialN+1:trialN*i);

    % thresholding 
    % set 10 percentile mark as the lower limit
%     llim = prctile(blkWiseSection(...
%         blkWiseSection>0), 10);
%     disp(llim)
%     disp(size(min(blkWiseSection(...
%         blkWiseSection>=llim))));
%     threshold(i) = min(blkWiseSection(...
%         blkWiseSection>=llim));

    earned(i) = sum(blkWiseSection);
    max_earned(i) = sum(blkMaxEarned);
end