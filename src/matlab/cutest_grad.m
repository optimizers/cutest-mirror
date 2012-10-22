function [varargout] = cutest_grad( varargin )
% Return function gradient.
% Usage:  f = cutest_grad(x).
    varargout = cell(1,nargout);
    [varargout{:}] = mcutest('grad',varargin{:});
