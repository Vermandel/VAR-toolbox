function [Mdl] = fn_select_best_orderVAR(Mdl,Y,CritName)
%FN_SELECT_BEST_ORDER Summary of this function goes here
%   Detailed explanation goes here

	if nargin < 3
		% likelihood comparison test
		CritName = 'LL';
	end
	p		= 0;
	pValue	= 0;

	Mdl.AR  = {nan(Mdl.NumSeries)};
	
	while pValue < 0.05
		p 		   = p + 1;
		Mdl.AR{p}  = nan(Mdl.NumSeries);
		[EstSpec,~,logL] = estimate(Mdl,Y,'Display','off');
		n1  = (Mdl.P-1)*Mdl.NumSeries^2;
		n1p = Mdl.P*Mdl.NumSeries^2;
		if p > 1 
			switch CritName
			
				case 'LL'
				
					[h,pValue]  = lratiotest(logL,logC,n1p-nC,0.05);
					disp(['lag:' num2str(p) ' LL:' num2str(logL) ' pVal:' num2str(pValue) ' Reject: ' num2str(h)])

				case {'AIC','BIC'}

					[aic,bic]  =  aicbic([logL logC],[n1p n1],size(Y,1));
					vc  = eval([lower(CritName) '(1)']);
					vnc = eval([lower(CritName) '(2)']);
					disp(['lag:' num2str(p) ' ' CritName 'nc:' num2str(vc) ' ' CritName ':' num2str(vnc) ])
					if vc < vnc
						h = 1;% reject new model
					else
						h = 0;
						break;
					end
					
					
			end
		else
			disp(['lag:' num2str(p) ' LL:' num2str(logL) ' (initial model)'])
		end 
		logC	= logL;
		nC		= n1p;
	end
	Mdl.AR{p} = zeros(Mdl.NumSeries); % remove last lag	
end

