function attr_accpt_correlation(alldat)
% ATTR_ACCPT_CORRELATION() finds the correlation b/w-
% 1. average acceptance % against a face
% 2. attractive score assigned to that face
% alldat: behav data collated

% all-acceptors
%AA_idx = ismember(1:40, [15 21 29 37:39]);
AA_idx = ismember(1:40, []);

% calc face spec acceptance
face_spec_accpt = face_Specific_Accept(alldat, 'photIDbehavDat.txt');

% calc face spec attr rate
attr_rating = importdata('attr_rating.csv');
face_spec_attr_rating = attr_rating.data(2:end, :); % row 1 = heading

% calc pearson's r, cuz why not
[rho, p] = corr(face_spec_accpt(~AA_idx,:), face_spec_attr_rating(~AA_idx,:));
rho = diag(rho); p = diag(p);

% group using emot = "emotional breakdown" hehe
photid = importdata('photIDbehavDat.txt');
photid = sortrows(photid, 4);
emotid = photid(:,1);

% plot bitchhh!
figure;
emot_lab = ["hap" "dis" "neu"];
col_lab = ["b" "r" "y"];

% rho
for emot = unique(emotid)'
    subplot(2, 3, emot);
    bar(rho(emotid == emot), col_lab(emot));
    xticklabels(photid(emotid == emot, 2));
    xlabel("emot spec photID"); ylabel("\rho");
    title('\rho | ' + emot_lab(emot));
    ylim([-.6 .6]);
end

% p
for emot = unique(emotid)'
    subplot(2, 3, 3 + emot);
    bar(p(emotid == emot), col_lab(emot));
    xticklabels(photid(emotid == emot, 2));
    yline(0.05);
    xlabel("emot spec photID"); ylabel("p");
    title('p | ' + emot_lab(emot));
    ylim([0 1])
end

sgtitle("pearson's \rho b/w attractiveness rating & acceptance rate");