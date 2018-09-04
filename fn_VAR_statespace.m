function [A,B] = fn_VAR_statespace(EstMdl)
%FN_VAR_STATESPACE
% Set the VAR(P) into a VAR(1)
% Y(t) = A + B*Y(t-1) + W;

	p = EstMdl.P;
	N = EstMdl.NumSeries;
	
	if isempty(EstMdl.Constant)
		EstMdl.Constant = zeros(N,1);
	end	
	
	%% STATE SPACE MODEL
	%A   = kron(EstMdl.Constant,[1;zeros(p-1,1)]);
	A	= repmat(EstMdl.Constant,p,1);
	% Constructing the B-matrix of the state-space model
	B1 = cell2mat(EstMdl.AR);
	B2 = eye(p*N);
	B  = [B1;B2(1:(end-N),:)];

end





