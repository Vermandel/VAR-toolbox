function [ irfs ] = fn_compute_IRF(EstMdl,ICp,y_names)
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
nb_simul 	= 50;
T_irf		= 50;

% compute the upper and lower bound
IC = norminv([ (1-ICp)/2 1-(1-ICp)/2],0,1);

% compute the steady state
simul_ss = filter(EstMdl,zeros(T_irf,EstMdl.NumSeries));%fn_gen_data(VarObj,zeros(T_irf,VarObj.n));

for i1 = 1:EstMdl.NumSeries % for each shock

	Y_simul		= repmat({nan(T_irf,nb_simul)},1,EstMdl.NumSeries);
	for i2 = 1:nb_simul % sampling response to draw intervals
		W_simul = zeros(T_irf,EstMdl.NumSeries);
		%W_simul(1,i1) = abs(normrnd(0,1,1,1)); % random N(0,1)
		%W_simul(1,:) = W_simul(1,:)*VarObj.Q^.5;
		W_simul(1,i1) = abs(normrnd(0,1,1,1))*EstMdl.Covariance(i1,i1)^.5;
		simul_irf = filter(EstMdl,W_simul);%fn_gen_data(VarObj,W_simul);
		for i3 = 1:EstMdl.NumSeries
			irf_all_var = Y_simul{i3};
			irf_all_var(:,i2) = simul_irf(:,i3)-simul_ss(:,i3);%(simul_irf(:,i3)-simul_ss(:,i1))./simul_ss(:,i1);
			Y_simul{i3} = irf_all_var;
		end
	end


	% plot result
	figure;
	for i2 = 1:EstMdl.NumSeries
		subplot(2,round(EstMdl.NumSeries/2),i2)	
		plot(1:T_irf,zeros(1,T_irf),'r:')
		hold on;
		IRFmean = mean(Y_simul{i2}')';
		IRFsd	= std(Y_simul{i2}')';
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

