@echo off
set loopcount=844
:loop
echo Iter %loopcount%
CloudCompare -SILENT -AUTO_SAVE OFF -O -GLOBAL_SHIFT -620000.00 -1000000.00 0 "D:\BCI_Spatial\UAV_Data\TiledPointClouds\BCI15Tiles_dec\BCI15d_%LOOPCOUNT%.laz"  -O -GLOBAL_SHIFT -620000.00 -1000000.00 0  "D:\BCI_Spatial\Lidar_Data\BCI09Tiles\BCI09_%LOOPCOUNT%.laz" -ICP -ITER 800 -OVERLAP 80 -C_EXPORT_FMT LAS -SAVE_CLOUDS FILE  "D:\BCI_Spatial\UAV_Data\TiledPointClouds\BCI15Tiles_aligned\BCI15a_%LOOPCOUNT%.las D:\BCI_Spatial\Lidar_Data\BCI09Tiles_ref\BCI09r_%LOOPCOUNT%.las" -CLEAR
set /a loopcount=loopcount-1
if %loopcount%==0 goto exitloop
goto loop
:exitloop
pause