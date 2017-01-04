function fits = getRegionFits(msd)

ss = numel(msd);

t = 1:ss;

ln_msd = log10(msd);
ln_t = log10(t);

[P,s] = polyfit(ln_t(~isnan(ln_msd)),ln_msd(~isnan(ln_msd)),1);

fits = P(1);


end