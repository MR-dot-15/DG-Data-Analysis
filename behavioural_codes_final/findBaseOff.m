function off_col = findBaseOff(off_col)
% spits the "base offer" \in {2, 3, 4}

off_col(off_col >= 15 & off_col < 25) = 2;
off_col(off_col >= 25 & off_col < 35) = 3;
off_col(off_col >= 35 & off_col < 45) = 4;