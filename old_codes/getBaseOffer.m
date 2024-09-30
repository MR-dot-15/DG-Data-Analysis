function baseOffer = getBaseOffer(offer_val)
if and(15<= offer_val, offer_val < 25)
    baseOffer = 2;
elseif and(25<= offer_val, offer_val < 35)
    baseOffer = 3;
elseif and(35<= offer_val, offer_val < 45)
    baseOffer = 4;
end