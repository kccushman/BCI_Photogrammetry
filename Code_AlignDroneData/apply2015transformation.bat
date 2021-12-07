@echo off
set loopcount=844
:loop
echo Iter %loopcount%
CloudCompare -SILENT -AUTO_SAVE OFF -O -GLOBAL_SHIFT -620000.00 -1000000.00 0 "D:\BCI_Spatial\UAV_Data\TiledPointClouds\BCI15Tiles\BCI15_%LOOPCOUNT%.laz" -APPLY_TRANS "D:\BCI_Spatial\UAV_Data\TiledPointClouds\BCI15Tiles_dec\BCI15mat_%LOOPCOUNT%.txt"  -C_EXPORT_FMT LAS -SAVE_CLOUDS FILE  "D:\BCI_Spatial\UAV_Data\TiledPointClouds\BCI15Tiles_alignedFull\BCI15af_%LOOPCOUNT%.las" -CLEAR
set /a loopcount=loopcount-1
if %loopcount%==0 goto exitloop
goto loop
:exitloop
pause