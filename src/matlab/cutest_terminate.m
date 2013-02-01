function [varargout] = cutest_terminate( varargin )
% Return variable names.
% Usage: vnames = cutest_terminate()
    varargout = cell(1,nargout);
    [varargout{:}] = mcutest('terminate',varargin{:});
