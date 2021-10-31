function [Yfit] = fit_gpr(X,Y,cov,hyp,Ncg)
% function [Yfit] = fit_gp(X,Y,ard)
%
% Fits a Gaussian Process with RBF kernel to the pairs (X,Y)
%
% INPUT:  X    should be N*d matrix (N data points, d dimensions)
%         Y    should be N*1 matrix (N data points)
%         ard  if 1, use covSEard, otherwise covSEiso (default)
%
% OUTPUT: Yfit contains the fitted y values
%
% NOTE:   Uses GPML 3.0 code (which should be in the matlab path)
%
% Copyright (c) 2008-2010  Joris Mooij  [joris.mooij@tuebingen.mpg.de]
% All rights reserved.  See the file COPYING for license terms.

% check input arguments
if size(Y,2)~=1 | size(X,1)~=size(Y,1)
    error('X should be Nxd and Y should be Nx1');
end
if regexp(cov,'covSEard')
    hyp(1)= log(hyp(1) * ones(size(X,2),1));
else
    hyp(1)=log(hyp(1));
end
covfunc = {'covSum',{cov,'covNoise'}};
hyp = minimize(hyp, 'gpr', -Ncg, covfunc, X, Y);
[Yfit]=gpr(hyp, covfunc,X, Y, X);
return

