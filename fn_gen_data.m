function [ Yhat,What ] = fn_gen_data(EstMdl,W,Y)
%FN_GEN_DATA Summary of this function goes here
%   Detailed explanation goes here

	T		= size(W,1);
	N		= size(W,2);
	p		= EstMdl.P;
	nT		= size(Y,1);
	What    = nan(size(Y));
	gap     = nT-T;
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%%% REPLICATE THE DATA
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	Yhat	= zeros(T+p,N);
	for t = 1:nT

		% adding the constant
		if ~isempty(EstMdl.Constant)
			Yhat(t,:) = Yhat(t,:)+EstMdl.Constant';
		end
		
		% Not using yet all AR process
		if t < (1+p)
			
			% adding AR-contributions
			for o = 1:(t-1)
				Yhat(t,:) = Yhat(t,:)+Yhat(t-o,:)*EstMdl.AR{o}';
			end		
			What(t,:) = Y(t,:) - Yhat(t,:);
			Yhat(t,:) = Yhat(t,:) + What(t,:);
		else % using all AR 
			
			% adding AR-contributions
			for o = 1:p
				Yhat(t,:) = Yhat(t,:)+Yhat(t-o,:)*EstMdl.AR{o}';
			end
			% adding the error term
			Yhat(t,:) = Yhat(t,:)+W(t-p,:);
			What(t,:) = W(t-p,:);	
			
		end
	end


end

