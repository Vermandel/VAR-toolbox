function [res] = fn_OLS_min(x,theY,theX,theb)
%OBJ_COTREND Summary of this function goes here
%   Detailed explanation goes here
	
	% select only NaN to estimate
	idx = find(isnan(theb));
	% set new parameters
	theb(idx) = x;
	% compute residuals
	res = theY - theX*theb;
	% square them plz
	res = res'*res;
	
end