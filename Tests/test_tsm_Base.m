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

