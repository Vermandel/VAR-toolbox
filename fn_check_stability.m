function [] = fn_check_stability(EstMdl)
%fn_check_stability Summary of this function goes here
%   Check whether VAR eigenvalues are below the unit circle

	% get the F-matrix of the state state VAR
	[~,F] 	= fn_VAR_statespace(EstMdl);
	disp(' ')
	disp('VAR Stability check')
	
	%get the absolute value of the roots
	myeig = abs(eig(F))';
	g=sprintf('%04.3g ',myeig);
	fprintf('Eigenvalues : %s\n', g)
	
	% display a warning in case of an unstable VAR
	if sum(myeig>1) > 0
		warning('This VAR is unstable')
	else
		disp('The VAR is stable')
	end

end

