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

	
	nburns = 500;		% number of paths
	Tpresample = 50;	% presample before the forecast
	
	%

	
	fce = nan(N,N,length(hori),nburns);
	for i1 = 1:nburns
		
		% draw future shocks after the forecast
		% to compute the contribution to the error
		rng(i1,'twister');		% fix the same seed across draws
		Z = normrnd(0,1,[max(hori) N]);%*EstMdl.Covariance;
		
		% simulate pre-sample before the forecast
		rng(i1+100,'twister');	% fix the same seed across draws
		Y0 	= filter(EstMdl,normrnd(0,1,Tpresample,N));
		
		% Run the forecast and feed the model with pre-samples
		%Y_fc 	= simulate(EstMdl,max(hori),'Y0',Y0);
		Y_fc 	= filter(EstMdl,0*Z,'Y0',Y0); % forecast using no shocks
		% Run the true (realized) value
		Ysimul = filter(EstMdl,Z,'Y0',Y0);
		
		% squared error forecast
		sq_err_fc = (Ysimul - Y_fc).^2;
		
		% for each contribution
		for i2=1:N
			
			% select the shock
			ZZ = zeros(size(Z));
			ZZ(:,i2) = Z(:,i2);
			% simulate with one shock
			Ysimul  = filter(EstMdl,ZZ,'Y0',Y0);
			sq_err_fc_i = (Ysimul - Y_fc).^2;
			
			% store the conditional variance
			fce(:,i2,:,i1) = sq_err_fc_i(hori,:)';
		
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

