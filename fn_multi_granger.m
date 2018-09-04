function [Fstat,Fthresh] = fn_multi_granger(Mdl,Y,id_y,id_x,alpha,y_names,meth)
%FN_MULTI_GRANGER Summary of this function goes here
%   Detailed explanation goes here

	if nargin <5
		alpha = 0.05;
	end
	if nargin < 7
		meth = 'VAR';
	end
	%id_x = 1;
	%id_y = 2;
	N = Mdl.NumSeries;
	p = Mdl.P;
	nT = size(Y,1);
	
	i = length(id_x);
	
	%% ESTIMATING UNRESTRICTED MODEL
	
	if strcmp(meth,'VAR')
	
		% VAR ESTIMATION
		[~,~,~,WU] = estimate(Mdl,Y);
		RSS_U = WU(:,id_y)'*WU(:,id_y);
		
	else
		% LEAST SQUARES
		% set Y = b*X + e
		theY=Y(1:(nT-p),id_y);
		theX = [ones(nT-p,1)];
		for i1 = 1:N
			for i2 = 1:p
				theX = [theX Y((1+i2):(i2-p+nT),i1)];
			end
		end
		% b mat
		theb = Mdl.Constant(id_y);
		for i2 = 1:p
			theMAT = Mdl.AR{i2};
			theb = [theb ; theMAT(id_y,:)'];
		end
		idx 	= find(isnan(theb));
		b_0    	= [zeros(size(idx))];		% guess initial values
		options = optimoptions('fmincon','Display','off','StepTolerance',1e-20);
		[x,RSS_U] = fmincon('fn_OLS_min',b_0,[],[],[],[],...
												[],[],[],options,theY,theX,theb);
	end	
   
   
	%% ESTIMATING RESTRICTED MODEL
	
	Mdl_R = Mdl;
	for i1 = 1:p
		ARmat = Mdl_R.AR{i1};
		ARmat(id_y,id_x) = 0;
		Mdl_R.AR{i1} = ARmat;
	end
	
	if strcmp(meth,'VAR')
	
		% VAR ESTIMATION
		[~,~,~,WR] = estimate(Mdl_R,Y);
		RSS_R = WR(:,id_y)'*WR(:,id_y);
		
	else
	
		% LEAST SQUARES
		% b mat
		theb = Mdl_R.Constant(id_y);
		for i2 = 1:p
			theMAT 	= Mdl_R.AR{i2};
			theb 	= [theb ; theMAT(id_y,:)'];
		end
		idx 	= find(isnan(theb));
		b_0    	= [zeros(size(idx))];		% guess initial values
		options = optimoptions('fmincon','Display','off','StepTolerance',1e-20);
		[x,RSS_R] = fmincon('fn_OLS_min',b_0,[],[],[],[],...
												[],[],[],options,theY,theX,theb);

	end
	

	
	%The numerator of the F-statistic
	F_num = ((RSS_R - RSS_U)/(p*i));

	%The denominator of the F-statistic
	F_den = RSS_U/(nT-p*i-(p+2));

	%The F-Statistic
	Fstat = F_num/F_den;

	Fthresh = finv(1-alpha,p,(nT-(p*i+p+2)));
	PVAL 	= 1 -  fcdf(Fstat,p,(nT-(p*i+p+2)));
	
	%%% GET NAMES
	for i1 = 1:i
		if nargin < 6
			thename = ['var' num2str(id_x(i1))];
		else
			thename = [ y_names{id_x(i1)}];
		end
		if i1==1
			names_x = thename;
		else
			names_x = [names_x ', ' thename];
		end
	end
	if nargin < 6
		names_y = ['var' num2str(id_y)];
	else
		names_y = y_names{id_y};
	end

	if Fstat > Fthresh
		fprintf('%s does Granger cause %s: Fstats %0.3g > %0.3g, pval=%0.3g \n',names_x,names_y,Fstat,Fthresh,PVAL)
	else
		fprintf('%s does NOT Granger cause %s: Fstats %0.3g < %0.3g, pval=%0.3g \n',names_x,names_y,Fstat,Fthresh,PVAL)
	end
end

