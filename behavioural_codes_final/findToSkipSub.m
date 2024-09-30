function to_skip = findToSkipSub(end_block, sublist)
% returns the sub idx to eliminate
% if first end_block many blocks are considered

to_skip = [];
sublist_idx = 1:length(sublist);
 
if end_block > 5
    to_skip = horzcat(to_skip, sublist_idx(...
        contains(sublist, ["mis" "vis"])));
end
if end_block > 10
    to_skip = horzcat(to_skip, sublist_idx(...
        contains(sublist, ["gau" "swa"])));
end
if end_block > 15
    to_skip = horzcat(to_skip, sublist_idx(...
        contains(sublist, ["kum" "biv"])));
end