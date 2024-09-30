function baseOff = findBaseOff(dat)
n = length(dat);
baseOff = zeros(n,1);

for i = 1:n
    baseOff(i) = getBaseOffer(dat(i));
end