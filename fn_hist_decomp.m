function [outputArg1,outputArg2] = fn_hist_decomp(EstMdl,W,T,y_names)
%FN_HIST_DECOMP Summary of this function goes here
%   Detailed explanation goes here
	
	
	if size(W,1) <  size(W,2)
		W = W';
	end

	N		= EstMdl.NumSeries;
	p		= EstMdl.P;
	nT		= size(W,1);
	T		= T((1+end-nT):end);
	if nargin <4
		for i1 = 1:N
			y_names{i1} = ['var' num2str(i1)];
		end
	end
	
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%%% COMPUTE THE HISTORICAL CONTRIBUTION
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% Historical contribution by shock
	Ycontrib	= repmat({zeros(nT,N)},1,3);
	% get ss
	simul_ss = filter(EstMdl,zeros(1,EstMdl.NumSeries));
	% computing the contribution of each shock to 
	% business cycles
	for i1 = 1:N % each shock
		% create an empty shock matrix that includes lags
		W2 = zeros(nT+p,N);
		% Feeding with one shock
		W2(:,i1) = [zeros(p,1);W(:,i1)];
		% computing the contribution of that shock
		%shock_contrib = fn_simul(EstMdl,W2)';
		shock_contrib = filter(EstMdl,W2);
		% removing pre-sample induced by lags
		shock_contrib = shock_contrib((1+p):end,:)-simul_ss;
		
		% each cell must include the contribution of each shock
		% for one endogenous variable
		for i2 = 1:N % each variable
			tempY			= zeros(nT,N);
			tempY(:,i1)		= shock_contrib(:,i2);
			Ycontrib{i2}	= Ycontrib{i2} + tempY;
		end
	end

	Yfilt = filter(EstMdl,[zeros(p,N);W]);
	
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%%% PLOT THE HISTORICAL CONTRIBUTION
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	time_split	= mean(diff(T)); % 1=years; 1/4 = quarters 
	linecolor	= char('b','r','g','m','c','k');
	% plotting
	for i1 = 1:N % each variable
		figure;
		hold on;
		thecontribs = Ycontrib{i1};
		% find contribution above/below 0
		contrib_above = (thecontribs>=simul_ss(i1)).*thecontribs;
		contrib_below = (thecontribs<=simul_ss(i1)).*thecontribs;
		for i2 = 1:nT % each period
			for i3 = 1:N % each shock
				if thecontribs(i2,i3) >= simul_ss(i1)
					ymini = simul_ss(i1)+sum(contrib_above(i2,1:i3))-thecontribs(i2,i3);
					ymax = simul_ss(i1)+sum(contrib_above(i2,1:i3));
				else
					ymini = simul_ss(i1)+sum(contrib_below(i2,1:i3))-thecontribs(i2,i3);
					ymax = simul_ss(i1)+sum(contrib_below(i2,1:i3));
				end
				fill([T(i2)-time_split/2;T(i2)-time_split/2;T(i2)+time_split/2;T(i2)+time_split/2],[ymini;ymax;ymax;ymini],linecolor(i3),'EdgeColor','none')
			end
		end
		%plot(T,Yfilt((1+p):end,i1)-0*simul_ss(i1),'k','LineWidth',1.5)
		plot(T,Yfilt((1+p):end,i1)-0*simul_ss(i1),'Color',[47 141 231]/255,'LineWidth',1.5)
		hold off;
		xlim([min(T) max(T)])
		title(['Historical decomposition of ' y_names{i1}])
		legend(y_names)
		grid on;
		save2pdf(['hist_decomp_' y_names{i1}])
	end

end

