%% Generate PDF document of TS Toolbox example
% in temporary directory "/tmp"
%
% Works only with Linux!

function MakeHTML( mfile )

publish( sprintf( '%s.m',mfile) ,'format','html',...
    'imageFormat','png',...
    'outputDir', 'html')

% Trim plots
png = dir( sprintf('html/%s_*.png', mfile ) );
for i=1:length(png)
    system( sprintf( '/usr/bin/convert -trim html/%s /tmp/xxx && mv /tmp/xxx html/%s', ...
        png(i).name,png(i).name ) )
end
