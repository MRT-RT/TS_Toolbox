%% Tests of tsm_DataSet object

close all
addpath( '../TSModel' ) % Path to TSModel files

%% Create data t,u,y
t = transpose( 0: 0.1 : 10 );
u = [ sin(t) , cos(t) ];
y = u(:,1) .* u(:,2);

%% Create object
d = tsm_DataSet( 'Test',u,y,'comment', {'Test','Validation'} );
% Add a time vector
d.set_t( t )
% Set labels for the signals
d.set_Label('t','Zeit')
d.set_Label('u',{'u1:sin','u2:cos'})
d.set_Label('y','Ausgang')
% Set limits
d.set_Limit( 'u', [-1.1, 1.1 ; -1.1, 1.1 ] )
d.set_Limit( 'y', [-0.6, 0.6 ] )
% Set unit
d.set_Unit( 'y', 'm/s' )
%%
plot(d,1);

%% Create two objects as copy
[t1,t2] = d.divide( 1,1 );

plot(t1,2);
plot(t2,3);
