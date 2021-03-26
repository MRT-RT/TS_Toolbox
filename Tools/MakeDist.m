%% Erstellen eines Archives

% $Id$

% Axel Dürrbaum, 10.5.2020 MRT / Uni Kassel <axeld@uni-kassel.de>

%% 
% Eingabe:
%  datei: Name des Zip-Archivs (Vorgabe: Name_Vorname-Datum-Zeit.zip)
%  verzeichnis: Zielverzeichnis für das Archiv (Vorgabe: ./Backups)
%%
% Ausgabe:
%  datei: vollständiger Pfad zum ZIP-Archiv

% Tip: git archive --format zip --output "output.zip" master -0

function datei = MakeDist( datei )

cd( '../..' )
verzeichnis = 'TS_Toolbox/Dists';

%% Backup als ZIP-Datei erstellen
if ~exist( verzeichnis, 'dir' )
    mkdir( verzeichnis )
end

%% Dateiname
if nargin < 1
    datei = sprintf('TS_Toolbox-%s-dist.zip',...
        datestr(now,'yyyy_mm_dd-HH_MM') );
end
if length(datei) < 4 || ~strcmp( datei(end-3:end), '.zip')
    datei = [ datei, '.zip' ];
end


%% Erstellen ZIP-Archiv

datei = sprintf('%s/%s', verzeichnis,datei );
zip( datei, ...
    { 'TS_Toolbox/TSModel/*.m' ; 
      'TS_Toolbox/TSModel/@tsm_Base/*.m' ; 
      'TS_Toolbox/TSModel/@tsm_Conc/*.m' ; 
      'TS_Toolbox/TSModel/@tsm_Prem/*.m' ; 
      'TS_Toolbox/TSModel/@tsm_DataSet/*.m' ; 
      'TS_Toolbox/Examples/*.m' ; 
      'TS_Toolbox/Examples/*.pdf' ; 
      'TS_Toolbox/Examples/Data/*' ;
      'TS_Toolbox/Functions/*.m' ; 
      'TS_Toolbox/*.pdf' } )

%%
dir( [ verzeichnis, '/*.zip' ] )
