%% MakeDoc: Erzeugen Dokumentation HTML->PDF

% $Id$

% Benötigt Linux/Python3/weasyprint

%feature('DefaultCharacterSet')
%feature('locale')
%char(unicode2native('ü','utf-8'))

cd( '../TSModel' )
publish( 'tsm_Manual.m','html')
if strcmp( computer, 'GLNXA64')
     system('weasyprint html/tsm_Manual.html ../TS_Toolbox-Manual.pdf')
end
