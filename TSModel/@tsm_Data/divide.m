%% Divide data for identification / validation
%
% Divides a tsdm_dataSet object into two reanges, 
% the size of which is determined by the ratio of s1 and s2
% A time vector will be corrected to start equally
%
% Inputs:
%  s1 proportion of identification data
%  s2 proportion of validation data
%
% Outputs:
%  d1 tsm_Database object with identification data
%  d2 tsm_Database object with validation data
%
% $Id$

function [ d1, d2 ] = divide( self, s1, s2, comment )

% sum of proportions
s = s1 + s2;

l1 = round( self.n * s1/s );
l2 = self.n - l1;
i2 = l1 + 1;

d1 = copy( self );
d1.n = l1;
d1.u = d1.u(1:l1,:);
d1.y = d1.y(1:l1);

d2 = copy( self );
d2.n = l2;
d2.u = d2.u(i2:end,:);
d2.y = d2.y(i2:end);

if isvector( self.t )
    d1.t=self.t(1:l1);
    d2.t=self.t(i2:end)-d2.t(i2);
end

if nargin > 3
    d1.addComment( comment );
    d2.addComment( comment );
end
d1.addComment( sprintf('Split source: ''%s'' / identification',self.Name) );
d2.addComment( sprintf('Split source: ''%s'' / validation',self.Name)  );

