function current_chan= create_new_chanlist(EEG_chanlocs)
% Input:
% 1:57 channel labels
% EEG.chanlocs file
% Returns the current channel labels: with or without blank

len = length(EEG_chanlocs);
current_chan = cell(len, 1);

for idx = 1:len
    temp = EEG_chanlocs(idx);
    current_chan{idx} = string(temp.labels);
end

current_chan = string(current_chan);