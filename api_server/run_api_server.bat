@echo off
setlocal

set "GNSS_CONDA_ENV=gnss-matlab-win"
if not "%~1"=="" set "GNSS_CONDA_ENV=%~1"

if not defined MATLAB_PREWARM set "MATLAB_PREWARM=1"

if defined CONDA_EXE (
    call "%CONDA_EXE%" activate "%GNSS_CONDA_ENV%"
) else if exist "%USERPROFILE%\anaconda3\condabin\conda.bat" (
    call "%USERPROFILE%\anaconda3\condabin\conda.bat" activate "%GNSS_CONDA_ENV%"
) else if exist "%USERPROFILE%\miniconda3\condabin\conda.bat" (
    call "%USERPROFILE%\miniconda3\condabin\conda.bat" activate "%GNSS_CONDA_ENV%"
) else (
    echo Conda activation script not found. Continuing with the current Python environment.
)

cd /d "%~dp0\.."

python -m uvicorn api_server.app:app --host 0.0.0.0 --port 8000 --reload

endlocal
