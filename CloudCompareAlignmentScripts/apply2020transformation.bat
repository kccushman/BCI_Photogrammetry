@echo off
set loopcount=844
:loop
echo Iter %loopcount%
CloudCompare -SILENT -AUTO_SAVE OFF -O -GLOBAL_SHIFT -620000.00 -1000000.00 0 "D:\BCI_Spatial\UAV_Data\TiledPointClouds\BCI20Tiles\BCI20_%LOOPCOUNT%.laz" -APPLY_TRANS "D:\BCI_Spatial\UAV_Data\TiledPointClouds\BCI20Tiles_dec\BCI20mat_%LOOPCOUNT%.txt"  -C_EXPORT_FMT LAS -SAVE_CLOUDS FILE  "D:\BCI_Spatial\UAV_Data\TiledPointClouds\BCI20Tiles_alignedFull\BCI20af_%LOOPCOUNT%.las" -CLEAR
set /a loopcount=loopcount-1
if %loopcount%==0 goto exitloop
goto loop
:exitloop
pause