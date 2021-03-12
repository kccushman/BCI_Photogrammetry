@echo off
set loopcount=844
:loop
echo Iter %loopcount%
CloudCompare -SILENT -AUTO_SAVE OFF -O -GLOBAL_SHIFT -620000.00 -1000000.00 0 "D:\BCI_Spatial\UAV_Data\TiledPointClouds\BCI17Tiles\BCI17_%LOOPCOUNT%.laz" -APPLY_TRANS "D:\BCI_Spatial\UAV_Data\TiledPointClouds\BCI17Tiles_dec\BCI17mat_%LOOPCOUNT%.txt"  -C_EXPORT_FMT LAS -SAVE_CLOUDS FILE  "D:\BCI_Spatial\UAV_Data\TiledPointClouds\BCI17Tiles_alignedFull\BCI17af_%LOOPCOUNT%.las" -CLEAR
set /a loopcount=loopcount-1
if %loopcount%==0 goto exitloop
goto loop
:exitloop
pause