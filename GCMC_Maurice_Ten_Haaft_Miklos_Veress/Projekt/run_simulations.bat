@echo off
if not exist results mkdir results

bin\GCMC_VS.exe results\test1.txt 0.56  100 100 100 37412
bin\GCMC_VS.exe results\test2.txt 0.84  100 100 100 87699


echo Simulations started.
pause