%% Normzahl


function n = Normzahl( unten, oben, R )

if nargin < 2
    unten = 0
    oben = 21
end
    
delta = oben - unten;

if nargin < 3
    % R10
    R = [1,1.25,1.6,2,2.5,3.15,4,5,6.3,8,10 ];
end

n = fix( log(delta)/log(10) );
zn =10^n;
x = delta / zn;

u = R( find( unten/zn <= R, 1 ) ) * 10^n
o = R( find( oben/zn <= R, 1 ) ) * 10^n
s = R( find( x <= R, 1 ) ) * 10^n;


