@echo off
setlocal

REM 외부 인자 받아오기
REM %1: input_dir
REM %2: nav_dir
REM %3: output_dir

if "%~3"=="" (
    echo run_processing.bat [input_dir] [nav_dir] [output_dir]
    exit /b 1
)

REM MATLAB 경로 (수정 필요 시)
set MATLAB_PATH="C:\Program Files\MATLAB\R2024a\bin\matlab.exe"

REM 현재 .bat 파일 위치 기준으로 goGPS_loadRinex 경로로 이동
REM 현재 디렉토리를 가져와서 저장
set "CALL_DIR=%~dp0"
cd /d "%CALL_DIR%\.."

REM 상대 경로 기준으로 goGPS_loadRinex 폴더 설정
set "SCRIPT_DIR=%CD%\goGPS_loadRinex"

REM MATLAB 함수 실행 - 인자 전달
%MATLAB_PATH% -r "process_and_save_rinex_files('%~1', '%~2', '%~3')"

endlocal
