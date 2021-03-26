#!/bin/sh

cd ../../

zip -r TS_Toolbox/Dists/TS_Toolbox-$(date '+%Y-%m-%d').zip \
TS_Toolbox/TSModel/*.m \
TS_Toolbox/TSModel/@tsm_Base/*.m \
TS_Toolbox/TSModel/@tsm_Prem/*.m \
TS_Toolbox/TSModel/@tsm_Conc/*.m \
TS_Toolbox/TSModel/@tsm_DataSet/*.m \
TS_Toolbox/Examples/*.m \
TS_Toolbox/Examples/*.pdf \
TS_Toolbox/Examples/Data*.mat \
TS_Toolbox/*.txt \
TS_Toolbox/*.pdf