function [keys, face_specific_off, face_off_array_mod] = check_Face_Specific_Off(behavDat)
face_specific_off = cell(1, 24);
keys = zeros(1, 24);

i = 1;
for emot = 1:3
    for ind = 1:8
        key = 10*emot + ind;
        keys(i) = key; i = i + 1;
        face_specific_off{key} = 0;
    end
end

trials = length(behavDat);
for i = 1:trials
    key = 10*behavDat(i, 1) + behavDat(i, 5);
    face_specific_off{key} = face_specific_off{key}...
        + behavDat(i, 2);
end

face_off_array = [face_specific_off{:}];

happ_mean = mean(face_off_array(1:8));
disg_mean = mean(face_off_array(9:16));
neut_mean = mean(face_off_array(17:24));
fprintf("Mean off-s:\nHapp: %.2f  Disg: %.2f  Neut: %.2f\n",...
    happ_mean, disg_mean, neut_mean);

face_off_array_mod = [face_off_array(1:8)./happ_mean...
    face_off_array(9:16)./disg_mean face_off_array(17:24)./neut_mean];

bar(keys, face_off_array_mod, 'cyan');
ylabel("Offer (normalized)"); grid("on");
title("Face specific offers");
text(13, 1.3, "Happ");
text(23, 1.3, "Disg");
text(33, 1.3, "Neut");
text(28, 1.4, "\downarrow");
ylim([0 1.5]);

% detect outlier
indx = detect_outlier(face_off_array_mod);