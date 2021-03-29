%% TS-Toolbox: ToDo-Liste
%
% Ideen und Verbesserungen


%% Prem/Gauss: individuelles sigma für jedes lokale Modell
%%
% * sigma = Skalar -> sigma(1:nv)=sigma0
% * sonst =sigma(1:nv)
%%
% 29.3.2029/ad

%% Idee: Scheduling
% negativer Lag -> gradient
%%
% 
% * lag=-1 -> u = gradient(y)/dt
% * lag=-2 -> u = gradient( gradient(y)/dtITEM2
%%
% Alternative: neues u für Ableitungen aufnehmen statt der Lags
%%
% 29.3.2029/ad

%% Zeiger auf Benutzer-Bewertungsfunktion
%%
% * für Optimierung
% * für Modellauswahl (Wrapper-Ansatz)
%%
% 29.3.2029/ad

%% Local Models: Init ARX/OE parameters
%%
% * random 
% * from other TS model
%%
% 29.3.2029/ad

%% Offnen Punkte
% 
%% 
% * Zeitvektor $t$ immer äquidistant / t_s konstant?

%% Befehle
% git clone /mrt/Software/Repos/git/TS_Toolbox.git