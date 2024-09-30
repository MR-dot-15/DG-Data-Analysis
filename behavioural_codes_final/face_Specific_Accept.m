function out = face_Specific_Accept(behavDatCell, colID, photID)
% FACE_SPECIFIC_ACCEPT(behavDatCell, photID)
% Finds mean acceptance of a given subject against a given face
% **Input:**
% behavDatCell: output of prep_Cell_with_Subj_Data.m
% colID: 4 if RT, 7 if acceptance
% photID: prep_Cell_with_Subj_Data.txt
%   
% **Output:**
% accpt: matrix (# subjects X # face) containing the mean acceptance val

% id of photo files used
photID = importdata(photID);
% sort based on the order in the form
% so that the sequence is maintained
photID = sortrows(photID, 4);
phot_n = length(photID);

% subjects
subj_n = length(behavDatCell);

% data holder
out = zeros(subj_n, phot_n);

for subj_i = 1:subj_n
    ith_behavdat = behavDatCell{subj_i};
    for phot_i = 1:phot_n
        % emotion and index of the given photo
        emot = photID(phot_i,1);
        index = photID(phot_i, 2);
        % all emots and indices in the behav dat
        emot_col = ith_behavdat(:,1);
        index_col = ith_behavdat(:,5);
        % calculate output
        if colID == 7
            out_temp = ith_behavdat(and(emot_col==emot, ...
                index_col==index),colID);
            out_temp = sum(out_temp == 1)/length(out_temp);
        elseif colID == 4
            out_temp = median(log(ith_behavdat(and(emot_col==emot, ...
                index_col==index),colID)));
        else
            error("supported colID: 4, 7");
        end

        % update data holder
        out(subj_i, phot_i) = out_temp;
    end
end