@echo off
setlocal

set "MATLAB_PATH=C:\Program Files\MATLAB\R2025a\bin\matlab.exe"
if not exist "%MATLAB_PATH%" (
    set "MATLAB_PATH=C:\Program Files\MATLAB\R2024a\bin\matlab.exe"
)

if not exist "%MATLAB_PATH%" (
    echo MATLAB executable not found.
    exit /b 1
)

set "CALL_DIR=%~dp0"
cd /d "%CALL_DIR%\.."

if not "%~1"=="" set "SAT_ENDPOINT_TSV_FILE=%~1"
if not "%~2"=="" set "SAT_ENDPOINT_NAV_FILE=%~2"
if not "%~3"=="" set "SAT_ENDPOINT_OUTPUT_TSV=%~3"

"%MATLAB_PATH%" -batch "cd('C:/Users/mskim/Desktop/workspace/goGPS_loadRinex/endpoint_satellite'); run_tsv_nav_endpoint_batch"

endlocal
