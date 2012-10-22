function [varargout] = cutest_setup( varargin )
% Set up main problem structure
% Usage:  prob = cutest_setup()
    varargout = cell(1,nargout);
    [varargout{:}] = mcutest('setup',varargin{:});
