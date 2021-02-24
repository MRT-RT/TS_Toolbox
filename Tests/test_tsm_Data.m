%% Test class tsm_Base 

clear
clc

addpath( '../TSModel' )

o1 = tsm_Base( 'Test' )
o2 = tsm_Base( 'Test', 'Author: AxelD' )

o3 = tsm_Base( 'Test', { 'Author: AxelD', 'EMail: axeld@uni-kassel.de'} )
o3.addComment( 'URL: http://www.uni-kassel.de/mrt' )

o3.save( 'Test.mat' )
o3.setDate( datetime( 'tomorrow' ) )
o3

%% Test Signal

t = 0:0.1:10;
u = sin( t );

su = tsm_Signal( 'inputs', u )
su.addComment( 'synthetisches Signal!' )
su.addComment( { 'Author: Axel Dürrbaum', ... 
                 'EMail: axeld@uni-kassel.de', ... 
                 'URL: http://www.uni-kassel.de/mrt'} )
su.Unit = 'V';
su.Limit = [-2,2];
plot( su )

%% Test DatSet
su1 = tsm_Signal( 'input 1', sin( t ) )
su1.Unit = 'v';
su2 = tsm_Signal( 'input 2', exp( sin( t ) ) )
su2.Unit = 'v/s';
sy = tsm_Signal( 'output', cos( t ) )
sy.Unit = 'm';

ds = tsm_DataSet( 'ident', [su1,su2], sy, 'Ident dataset' )
figure(1),clf
plot( ds );