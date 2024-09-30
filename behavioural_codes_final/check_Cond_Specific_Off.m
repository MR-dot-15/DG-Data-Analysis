function [off_cond_spec_extendedList, mid_points, mean_offers] = check_Cond_Specific_Off()
tot = 360;
conds = [12 13 14 22 23 24 32 33 34];
off_cond_spec_extendedList = zeros(1, tot+16);
mean_offers = zeros(1, 9);
mid_points = zeros(1, 9);

i = 1; j = 1;
for cond = conds
    load('offer_condition_specific.mat', 'off_cond_spec');
    temp = off_cond_spec{cond};
    n = length(temp);
    off_cond_spec_extendedList(i:i+n-1) = temp;
    mean_offers(j) = mean(temp);
    mid_points(j) = i+n/2;

    text(i+n/3, 10*rem(cond,10)+7, num2str(cond)); hold on;

    i = i+n+2; j = j + 1;
end

bar(off_cond_spec_extendedList); 
plot(mid_points, mean_offers, 'LineWidth', 1.5);
ylim([0 50]); ylabel("Offers"); title("Condition specific offers");