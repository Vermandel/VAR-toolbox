function [Yhat] = fn_simul(EstMdl,W)
%FN_SIMUL_SS Summary of this function goes here
%   Detailed explanation goes here
	
	if size(W,1) >  size(W,2)
		W = W';
	end

	N		= EstMdl.NumSeries;
	p		= EstMdl.P;

	[A,B] 	= fn_VAR_statespace(EstMdl);
	nB		= size(B,1);
	nT		= size(W,2);
	Yhat 	= repmat(EstMdl.Constant,p,nT+p);
	
	for t = (1+p):size(Yhat,2)
		
		%% add lag contribution
		Yhat(:,t) = B*Yhat(:,t-1) + [W(:,t-p);zeros((p-1)*N,1)];
		
	end
	Yhat=Yhat(1:N,(1+p):end);
end

