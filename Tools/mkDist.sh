#!/bin/sh

zip -r TS_Toolbox-$(date '+%Y-%m-%d.zip) \
TS_Toolbox/TSModel/*.m \ 
TS_Toolbox/TSModel/@tsm_Base*.m \ 
TS_Toolbox/TSModel/@tsm_Prem\*.m \ 
TS_Toolbox/TSModel/@tsm_Conc\*.m \ 
TS_Toolbox/TSModel/@tsm_DataSet\*.m \ 
TS_Toolbox/*.txt
TS_Toolbox/*.pdf
