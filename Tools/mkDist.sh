#!/bin/sh

cd ../../
if test -d TS_Toolbox/Dists; then
 dest=TS_Toolbox/Dists
else
 dest=. 
fi
zip -r $dest/TS_Toolbox-$(date '+%Y-%m-%d').zip \
TS_Toolbox/TSModel/*.m \
TS_Toolbox/TSModel/@tsm_Base/*.m \
TS_Toolbox/TSModel/@tsm_Prem/*.m \
TS_Toolbox/TSModel/@tsm_Conc/*.m \
TS_Toolbox/TSModel/@tsm_DataSet/*.m \
TS_Toolbox/Examples/*.m \
TS_Toolbox/Examples/*.pdf \
TS_Toolbox/Examples/Data/*.mat \
TS_Toolbox/*.txt \
TS_Toolbox/*.pdf
