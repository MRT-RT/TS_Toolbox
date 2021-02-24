
for nD=[3,5]
for n=[5,6]
    %N = 6 % x in [0,1]
    N=n^nD;

x0 = lhsdesign(N,nD,'criterion','maximin') ;
yx0 = Friedmann_Function(x0,nD);

save(sprintf('LHS-Data-nD%d-N%d.mat',nD,N),'x0','yx0')
end
end
