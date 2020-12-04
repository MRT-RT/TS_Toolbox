%% MakeDoc: Erzeugen Dokumentation HTML->PDF

% $Id$

% Benötigt Linux/Python3/weasyprint

cd( '../TSModel' )
publish( 'tsm_Manual.m','html')
if strcmp( computer, 'GLNXA64')
     system('weasyprint html/tsm_Manual.html ../tsm_Manual.pdf')
end
