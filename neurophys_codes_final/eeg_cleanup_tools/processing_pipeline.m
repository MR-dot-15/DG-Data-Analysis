%% standardized pipeline
%% first take backup
odat = EEG.icaact;
oeegdat = EEG.data;
chan = create_new_chanlist(EEG.chanlocs); chan_idx = 1:numel(chan);

%% manipulate the IC / ERPs
%dat = EEG.icaact;
dat = EEG.data;

% target_chan_set = chan_idx(contains(chan, ...
%                   ["T8" "P10" "FT8"]));
% 
target_chan_set = 1:numel(chan);

% trend correction
for elec_idx = target_chan_set
    dat = trend_correction(dat, elec_idx, chan, 1, 0);
end


% oscillation correction
% for elec_idx = target_chan_set
%     dat = remove_sin(dat,...
%         elec_idx, chan, 0, 1, 0);
% end

%%
EEG.data = dat;

%% ICA to EEG data
dat_eeg = EEG.icawinv * reshape(EEG.icaact,...
                                size(EEG.icaact,1), []);
dat_eeg = reshape(dat_eeg, size(EEG.data));
EEG.data = dat_eeg;