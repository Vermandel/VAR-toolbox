function [ irfs ] = fn_compute_IRF(EstMdl,ICp,y_names,W)
%FN_COMPUTEIRF Summary of this function goes here
%   Detailed explanation goes here




if nargin < 2
	% if no uncertainty bounds
	% set if to 68%
	ICp		= .68;
end

if nargin < 3
	% if no y names
	% set one
	for i1 = 2:EstMdl.NumSeries
		y_names{i1} = ['var' num2str(i1)];
	end
end

%%% Generating the random variable
%%% options:
nb_simul 	= 200;
nb_simuly 	= 1000;
T_irf		= 50;

% compute the upper and lower bound
IC = norminv([ (1-ICp)/2 1-(1-ICp)/2],0,1);

irf_all_var = nan(T_irf,EstMdl.NumSeries,nb_simul,EstMdl.NumSeries);

for i2 = 1:nb_simul % sampling response to draw intervals
	
	rng(i2);
	
	% artificial draws
	simul_y_hat = filter(EstMdl,randn(nb_simuly,EstMdl.NumSeries));%fn_gen_data(VarObj,zeros(T_irf,VarObj.n));
	%simul_y_hat = filter(EstMdl,W(randi(size(W,1),nb_simuly,1),:));%fn_gen_data(VarObj,zeros(T_irf,VarObj.n));
	
	% re-estimate the model
	[EstMdl_i] = estimate(Mdl,simul_y_hat);
	
	
	% Randomize constants
	%EstMdl_i.Constant = EstMdl.Constant + EstStdErrors.Constant.*randn(EstMdl.NumSeries,1);
	%% Randomize AR
	%for p = 1:EstMdl.P
	%	EstMdl_i.AR{p} = EstMdl.AR{p} + EstStdErrors.AR{p}.*randn(EstMdl.P,EstMdl.P);
	%end
	
	% compute asymptotic mean
	simul_ss = filter(EstMdl_i,zeros(T_irf,EstMdl_i.NumSeries));

	for i1 = 1:EstMdl.NumSeries % for each shock	
		W_simul = zeros(T_irf,EstMdl.NumSeries);
		W_simul(1,i1) = 1;
		W_simul(1,:) = W_simul(1,:)*chol(EstMdl_i.Covariance)';
		irf_all_var(:,:,i2,i1)     = filter(EstMdl_i,W_simul)-simul_ss;
	end
end


% plot result
for i1 = 1:EstMdl.NumSeries
	figure;
	for i2 = 1:EstMdl.NumSeries
		subplot(2,round(EstMdl.NumSeries/2),i2)	
		IRFnow = squeeze(irf_all_var(:,i2,:,i1));
		plot(1:T_irf,zeros(1,T_irf),'r:')
		hold on;
		IRFmean = mean(IRFnow,2);
		IRFsd	= std(IRFnow,[],2);
		IRFlb 	= IRFmean+IC(1)*IRFsd;
		IRFub 	= IRFmean+IC(2)*IRFsd;
		X=[1:T_irf,fliplr(1:T_irf)];                %#create continuous x value array for plotting
		Y=[IRFlb',fliplr(IRFub')];              	%#create y values for out and then back
		h=fill(X,Y,'black','edgecolor','none');
		set(h,'facealpha',.15) % opacity
		plot(1:T_irf,IRFmean,'-b')
		plot(1:T_irf,IRFub,'Color',[0.7,0.7,0.7])
		plot(1:T_irf,IRFlb,'Color',[0.7,0.7,0.7])
		grid on;	
		hold off;
		xlim([0 T_irf])
		%plot(1:T_irf,IRFmean,1:T_irf,IRFlb,1:T_irf,IRFub)
		title([y_names{i2} ' response to ' y_names{i1}]); 
		eval(['irfs.mean.shock_' y_names{i1} '.var_' y_names{i2} ' = IRFmean;'])
		eval(['irfs.ub.shock_' y_names{i1} '.var_' y_names{i2} ' = IRFub;'])
		eval(['irfs.lb.shock_' y_names{i1} '.var_' y_names{i2} ' = IRFlb;'])	
	end
	save2pdf(['IRF_' y_names{i1}])
end

end

