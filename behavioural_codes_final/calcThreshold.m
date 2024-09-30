function store_thresh = calcThreshold(alldat, emot_idx, llim_to_ulim, perc)
% in //
% alldat: 40*1 cell containing behavdata
% condition: [emot off_lower off_upper] 
% (NOTE: emot = 0 means across all emot)
% llim_to_ulim: range of trials = llim:ulim format
% out //
% threshold: n_sub * 1 array

permitted_error = .02;
store_thresh = zeros(length(alldat),1);

for sub_idx = 1:length(alldat)
    dat = alldat{sub_idx}; dat = dat(llim_to_ulim,:);
    accpt_rate = sum(dat(:,end)==1)/size(dat,1);

    % find the slicing index based on the condition array
    if emot_idx ~= 0
        slice_idx = dat(:,1) == emot_idx;   
    else
        slice_idx = 1:size(dat,1);
    end
    dat = dat(slice_idx,:); 

    % define threshold as mean(75th per of rej off, 
    %                          25th per of accp off)
    
    % for all rejectors 
    if accpt_rate < permitted_error
        store_thresh(sub_idx) = prctile(dat...
                        (dat(:,end)==-1,2),...
                        99);
        
    % for all acceptors
    elseif accpt_rate > 1-permitted_error
        store_thresh(sub_idx) = prctile(dat...
                        (dat(:,end)==1,2),...
                        1);
        
    % for "normal" people
    else
        off_accp_min = prctile(dat...
                        (dat(:,end)==1,2),...
                        perc);
        off_rej_max = prctile(dat...
                        (dat(:,end)==-1,2),...
                        80);
        if isnan(off_accp_min)
            off_accp_min = off_rej_max;
            %store_thresh(sub_idx) = 44;
        elseif isnan(off_rej_max)
            off_rej_max = off_accp_min;
        % else
        %     store_thresh(sub_idx) = off_accp_min;
        end
        store_thresh(sub_idx) = mean([off_accp_min off_rej_max]);
    end
end