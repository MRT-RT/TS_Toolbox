%% TS-Toolbox Symbols
%
% Symbolx used in the TS Toolbox
%%
% <latex>
% \def\lag{\mathrm{lag}}
% \begin{tabular}{llll}
% Sign & Meaning & Matlab & Comment \\\hline
% $N$ & number of data points & \verb|N| & \\
% $t$ & time vector   & \verb|t| & $N\times 1$ \\
% $t_s$ & sampling time  & \verb|dt| & 0=undefined\\
% $n_u$ & number of inputs  & \verb|nu| & \\
% $u$ & input  & \verb|u| & $N\times n_u$ \\
% $n_y$ & number of outputs  & \verb|ny| & always 1\\
% $y$ & output  & \verb|y| & $N\times n_y$ \\
% \hline
% $n_v$ & number of local models  & \verb|nv| & \\
% $v$ & cluster centers & \verb|v|  & \\
% $z$ & scheduling matrix  & \verb|z| & \\
% $z_{\lag,u}$ & scheduling input lags & \verb|z_lag_u| &\\ 
% $z_{\lag,y}$ & scheduling output lags & \verb|z_lag_u| &\\ 
% $\mu$ & membership degree  & \verb|mu| & \\
% $\phi$ & membership degree  & \verb|phi| & $\phi=\sum_{i=1}^{n_v} \mu_i$\\
% $\nu$ & FBF exponent  & \verb|nue| & $\nu\in\{1,2\}$\\
% $\sigma$ & Gauss exponent  & \verb|sigma| & \\
% \hline
% $x$ & regression matrix & \verb|x| &\\ 
% $x_{\lag,u}$ & regresser input lags & \verb|x_lag_u| &\\ 
% $x_{\lag,y}$ & regresser output lags & \verb|x_lag_u| &\\ 
% \hline
% $A$ & local models y & \verb|A| &\\
% $B$ & local models u & \verb|B| &\\
% $c$ & local models affine term & \verb|c|  &
% \end{tabular}
% </latex>
