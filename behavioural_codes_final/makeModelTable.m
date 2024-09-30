function tbl = makeModelTable(alldat)
% collate the entire data into one matrix
% thresh_mat: threshold val over block 
% for 40 subs (40 x 460)
collated_dat = zeros(1,7);
subID = [];
blockID = [];
thresh = [];
for sub_idx = 1:numel(alldat)
    dat = alldat{sub_idx};
    temp_blockID = reshape(repmat(1:20, 18, 1), 1, []);

    % RT filtering
    slice_idx = dat(:,4) ~= 2 & dat(:,4) > .25...
                & dat(:,4) < prctile(dat(:,4), 97);
    collated_dat = cat(1, collated_dat, dat(slice_idx,:));
    subID = cat(2, subID,...
        sub_idx * ones(1, sum(slice_idx)));
    blockID = cat(2, blockID, temp_blockID(slice_idx));
end
collated_dat = collated_dat(2:end, :);

% make the table
% (collated_dat(:,2)-min(collated_dat(:,2)))...
%                             ./max(collated_dat(:,2))

tbl = table(collated_dat(:,1),...
            collated_dat(:,1),...
            collated_dat(:,2),...
            subID', blockID', log(collated_dat(:,4)),...
            zeros(size(subID')), (collated_dat(:,7)+1)./2,...
            'VariableNames', ["emot" "emot_og" "off" "subID"...
                              "bID" "rt" "rt_res" "accpt"]);
tbl.emot = rem(tbl.emot + 1, 4);
tbl.emot(tbl.emot == 0) = 1;
tbl.emot = categorical(tbl.emot);
tbl.off = (tbl.off - median(tbl.off))./(max(tbl.off) - min(tbl.off));
tbl.subID = categorical(tbl.subID);
tbl.bID = categorical(tbl.bID);