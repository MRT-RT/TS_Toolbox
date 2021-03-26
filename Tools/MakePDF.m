%% Generate PDF document of TS Toolbox example
% in temporary directory "/tmp"
%
% Works only with Linux!

function makePDF( mfile )

publish( sprintf( '%s.m',mfile) ,'format','latex',...
    'imageFormat','png',...
    'stylesheet','../Tools/mxdom2latex.xsl',...
    'outputDir', '/tmp')

% Trim plots
png = dir( sprintf('/tmp/%s_*.png', mfile ) );
for i=1:length(png)
    system( sprintf( '/usr/bin/convert -trim /tmp/%s /tmp/xxx && mv /tmp/xxx /tmp/%s', ...
        png(i).name,png(i).name ) )
end

% Run twice for TOC
system( sprintf( 'cd /tmp && pdflatex %s && pdflatex %s', mfile, mfile ) )

% Move resulting PDF to actual directory
movefile( sprintf( '/tmp/%s.pdf', mfile ), '.' )

% % cleanup LaTeX aux files
% for ext = { 'aux','log','out','toc' }
%     delete( sprintf('PDF/%s.%s', mfile,ext{1} ) )
% end
% mfile.tex, mfile_*.eps/mfile_*-eps-converted-to.pdf    