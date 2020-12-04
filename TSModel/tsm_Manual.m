%% TS-Toolbox
%
% Matlab-Toolbox zur nichtlinearen Systemidentifikation mittels lokal affiner Tagaki-Sugeno-Modelle
%
% <<Titelbild-TSM.jpg>>
% 
% Version: 1.3 vom 14.9.2020
%
% Prof. Dr.-Ing. Andreas Kroll, FG Mess- und Regelungstechnik, FB 15 Maschinenbau, Universität Kassel
% 
% <<Logo-MRT.png>>
%
% URL: <http://www.uni-kassel.de/go/mrt>
%
% Author: Axel Dürrbaum (<mailto:axel.duerrbaum@mrt.uni-kassel.de>)


%% Aufgabe: Nichtlineare Systemidentifikation und Regression
%%
% * für statische MISO-Modelle
%%
% $$y(t) = f( u_1(t),\ldots,u_m(t) )$$
%%
% * oder dynamische MISO-Modelle
%%
% $$y(t) = f(u_1(t),\ldots,u_1(t-m_1),\ldots,u_m(t-1),\ldots,u_m(t-m_m),\ldots, y(t-1),\ldots,y(t-n)$$
%
%% Modellansatz: lokal affine Tagaki-Sugeno-Modelle (TS)
%
% Überlagerung der $c$ lokal affinen Teilmodelle $y_i(x)$ zu einem Gesamtmodell
%%
% $$\hat{y}(t) = \sum_{i=1}^c \mu_i(z) \cdot \hat{y}_i(x)$$
%% 
% * mit den Eingangssignalen $u(t)$ und dem Ausgangssignal $y(t)$,
% * der Scheduling-Variablen $z(u,y)$, 
% * der Zugehörigkeitsfunktionen $\mu_i(z)$,
% * der Regressor-Variablen $x(u,y)$
% * und den lokalen TS-Modellen $\hat{y}_i(x)$ 


%% Funktionsprinzip
%%
% * Datensatz $\{ u(t), (y(t) \}$, ggf. Normierung und Split in Identifikations- und
% Validierungsdaten
% * Ggf. Anpassung der Standard-Einstellungen der Hyperparameter
%   * Clusterung $\nu$, $c$, $\epsilon_{FCM}$
%   * Multistart
%   * NL-Optimierung
% * Vorgabe der Anzahl der lokalen Modelle $c$ und des Unschärfeparameters
% $\nu$ 
% * Clustering zur Emittlung der Partitionierung bzw. Lage der Teilmodelle im Scheduling-Raum,
% Multistartstrategie mit Auswahl des besten Ergebnisses auf Basis des
% Modellfehlers auf Identtifikationsdaten
% * Initiale Schätzung der lokalen Modelle mittels Least-Squares-Verfahren (lokal
% oder global)
% * Optionale Optimierung der Zugehörigkeitsfunktionen $\mu_i$ und/oder der lokalen Teilmodelle $\hat{y}_i$ 
% mittels nichtlinearer Optimierung der Simulation (Matlab-Funktion |lsqnonlin|)
% * Unterschiedliche Wahl der Scheduling- und Regeressor-Variablen möglich
% * Validierung auf neuen Daten

%% Clusterung
%
% Eingangs- ($u$) oder Produktraum $(u|y)$
%
% Implementierte Algorithmen:
%%
% * Abstandsnormen:  Euklid, Mahalanobis
% * Fuzzy C-Means (FCM)
% * Gustafson-Kessel (GK)

%% Verfügbare Zugehörigkeitsfunktionen
%%
% * FCM-Type-Funktionen
% * Gauss-Funktionen

%% Lokale TS-Modelle
%%
% * Linear: $y(t) = \sum_{i=0}^n a_i\cdot u_i(t) + a_0$
% * ARX: $y(t) =  A\cdot y(t-dt) + B\cdot u(t) + C$
% * OE:  $y(t) =  A\cdot y(t-dt) + B\cdot u(t) + C + e(t)$

%% Modellgütemaße
%
% auf Identifikations- und Validierungsdaten
%%
% * Maximum Absolute Error (MAE)
% * Sum of Squared Errors (SSE)
% * Mean Squared Error (MSE)
% * Root Mean Squared Error (RMSE)
% * Normalized Mean Squared Error (NMSE)
% * Best Fit Rate (BFR)
% * ??? Akaike Information Criterion (AIC)
% * ??? Bayesian Information Criterion (BIC)

%% Visualisierung
%%
% * Clusterung (2D, n-dimensional als mehrfache 2D)
% * ??? Zugehörigkeitsmaße / Regelaktivierung
% * Residuen
% * ??? Residualhistogramm
% * Simulation oder 1-Schritt-Prädiktion auf Identifikations- oder Validierungsdaten

%% Dokumentation
%
%% Verfügbare Objekte 
%%
% * Daten
% Parameter, Methoden
% * Modell
% Parameter, Methoden

%% Verfügbare Funktionen
%
% Funktionen, die nicht auf Objekten arbeiten

%% Installation
%
% ToDo
%
% Verzeichnisse
%% 
% * TS_Toolbox
% * TS_Toolbox/Functions
% * TS_Toolbox/Examples

%% Musterprojekte
%
% im Verzeichnis |Examples|
%%
% * statisch: Akademisches Beispiel    |Test_LS_Akad|
% * statisch: Friedmann-Funktion 2D/3D |Test_LS_Friedman|
% * statisch: Kompressor-Kennlinie 3D  |Test_LS_Kompressor|
% * dynamisch: Narendra (SISO)         |Test_ARX_Narendra.m|
% * dynamisch: Regelkappe (SISO)       |Test_ARX_Throttle.m|
% * dynamisch: Drosselkappe IAV (MISO) |Test_ARX_Ladedruck.m|

%% Implementierung
%
% Objektorientierte Realisierung: 
%%
% * Objekt |Daten|: $u(t),y(t)$
% * Objekt |Modell|: Daten, Premisse, Konklusion
% * Objekt |Premisse|: Scheduling / Zugehörigkeitsfunktion
% * Objekt |Konklusion|: Regresser / lokale Modelle (ARX/OE)

%% Benötigte Software
%% 
% * Matlab R2019a oder höher (Windows//Linux/MacOS)
% * Matlab Fuzzy Toolbox (FCM)
% * Matlab Optimzation Toolbox (lsqnonlin)

%% Geplante Erweiterungen
%
% <html>
% <table>
% <th>Aufgabe<th>Zeitraum<th>Prio<th>Status<tr>
% <!--->
% <td>Toolbox als OO-Klasse in Matlab<br>
% Identifikation+Optimierung TS-Modelle<br>
% statisch + ARX/OE (mg)
% <td>1<td> 5/20<td>Erledigt 100%<tr>
% <!--->
% <td>ARX MISO-Modelle<br>
% Optimierung MISO<br>
% Tests Matlab-Optimierungsverfahren/-parameter
% <td>2<td>5/20 – 6/20<td>Erledigt: 100%<tr>
% <!--->
% <td>Subklassen für Modelle <br>
% (LS/ARX/OE) und Datensätze
% <td><td>6/20<td>60%<tr>
% <!--->
% <td>Testsignalentwurf (mg)
% <td><td>6/20 - 7/20<td>0%<tr>
% <!--->
% <td>Regelung (as)
% <td><td>7/20<td>0%<tr>
% <!--->
% <td>Strukturselektion (mk)
% <td><td><td>10%<tr>
% <!--->
% <td>Maximum-Likelihood (jf)
% <td><td><td>0%<tr>
% <!--->
% <td>BETS (fw)
% <td><td><td>0%<tr>
% <!--->
% <td>Datascreening (da)
% <td><td><td>0%<tr>
% <!--->
% <td>Implementierung weitere Cluster-Verfahren
% <td><td><td>0%<tr>
% </table>
% </html>

%%
% $Id$