function [] = fn_cond_forecast_error(EstMdl,hori,y_names)
%FN_COND_FORECAST_ERROR Summary of this function goes here
%   Detailed explanation goes here
	%%% Business Cycles Moments
	%% func
	N = EstMdl.NumSeries;
	hori = sort(hori);

	if nargin < 3
		% if no y names
		% set one
		for i1 = 2:EstMdl.NumSeries
			y_names{i1} = ['var' num2str(i1)];
		end
	end

	%% first sampling
	nburns = 50;
	fce = nan(N,N,length(hori),nburns);
	for i1 = 1:nburns
		% fix the same seed across draws
		rng(i1);
		Z = normrnd(0,1,[max(hori) N]);%*EstMdl.Covariance;

		% for each contribution
		for i2=1:N
			
			% select the shock
			ZZ = zeros(size(Z));
			ZZ(:,i2) = Z(:,i2);
			% simulate with one shock
			Ysimul = filter(EstMdl,ZZ);
			
			% compute the conditional variance
			for i3=1:length(hori)
				if hori(i3) == 1 % then cannot compute the var
					thevar = zeros(N,1);
					thevar(i2) = 1;
					fce(:,i2,i3,i1) = thevar;
				else
					fce(:,i2,i3,i1) = var(Ysimul(1:hori(i3),:))';
				end
			end
		
		end
	end
	% mean from different draws
	c_mean = mean(fce,4);
	% compute % contrib
	tt_mean = 100*c_mean./sum(c_mean,2);
	disp('VARIANCE FORECAST ERROR')
	for i1=1:length(hori)
		disp(' ')
		disp('')
		disp(['horizon ' num2str(hori(i1)) ' period(s)'])
		disp(array2table(tt_mean(:,:,i1),'VariableNames',y_names,'RowNames',y_names))
	end

end

