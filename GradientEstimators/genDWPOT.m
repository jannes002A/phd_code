%Generate DWPOT.m
clear all
n=1500;
leftlim=-3.0;
rightlim=3.0;
diff=rightlim-leftlim;
x=(leftlim:diff/(n-1):rightlim)';
clear n
BETA=4;
POT='asym2wellPot';
V=Pot(x,POT);

save DWPOT_.mat