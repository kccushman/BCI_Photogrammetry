@echo off
set loopcount=844
:loop
echo Iter %loopcount%
CloudCompare -SILENT -AUTO_SAVE OFF -O -GLOBAL_SHIFT -620000.00 -1000000.00 0 "D:\BCI_Spatial\UAV_Data\TiledPointClouds\BCI18Tiles_dec\BCI18d_%LOOPCOUNT%.laz"  -O -GLOBAL_SHIFT -620000.00 -1000000.00 0  "D:\BCI_Spatial\UAV_Data\TiledPointClouds\BCI15Tiles_ref2\BCI15r_%LOOPCOUNT%.las" -ICP -ITER 800 -OVERLAP 80 -C_EXPORT_FMT LAS -SAVE_CLOUDS FILE  "D:\BCI_Spatial\UAV_Data\TiledPointClouds\BCI18Tiles_alignedto15\BCI18a_%LOOPCOUNT%.las D:\BCI_Spatial\UAV_Data\TiledPointClouds\BCI15Tiles_ref2\BCI15r_%LOOPCOUNT%.las" -CLEAR
set /a loopcount=loopcount-1
if %loopcount%==0 goto exitloop
goto loop
:exitloop
pause