function y = LogPrior_v2(para,para_det)

cv = 0.05*abs(para_det);
y = sum(log(normpdf(para,para_det,cv)));