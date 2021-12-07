@echo off
set loopcount=844
:loop
echo Iter %loopcount%
CloudCompare -SILENT -AUTO_SAVE OFF -O -GLOBAL_SHIFT -620000.00 -1000000.00 0 "D:\BCI_Spatial\UAV_Data\TiledPointClouds\BCI20Tiles_dec\BCI20d_%LOOPCOUNT%.laz"  -O -GLOBAL_SHIFT -620000.00 -1000000.00 0  "D:\BCI_Spatial\UAV_Data\TiledPointClouds\BCI18Tiles_ref2\BCI18r_%LOOPCOUNT%.laz" -ICP -ITER 800 -OVERLAP 80 -C_EXPORT_FMT LAS -SAVE_CLOUDS FILE  "D:\BCI_Spatial\UAV_Data\TiledPointClouds\BCI20Tiles_alignedto18\BCI20a_%LOOPCOUNT%.las D:\BCI_Spatial\UAV_Data\TiledPointClouds\BCI18Tiles_ref2\BCI18r_%LOOPCOUNT%.las" -CLEAR
set /a loopcount=loopcount-1
if %loopcount%==0 goto exitloop
goto loop
:exitloop
pause