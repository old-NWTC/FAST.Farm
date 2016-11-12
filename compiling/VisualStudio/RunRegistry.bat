@ECHO OFF

set lines=======================================================================
echo %lines%
IF "%1"=="" (
ECHO.
ECHO   The calling syntax for this script is
ECHO             RunRegistry ModuleName
ECHO.
GOTO Done
)


REM ----------------------------------------------------------------------------
REM ------------------------- LOCAL PATHS --------------------------------------
REM ----------------------------------------------------------------------------
REM -- USERS MAY EDIT THESE PATHS TO POINT TO FOLDERS ON THEIR LOCAL MACHINES. -
REM -- NOTE: do not use quotation marks around the path names!!!! --------------
REM ----------------------------------------------------------------------------
REM ----------------------------------------------------------------------------

SET Root_Loc=..\..
SET Subs_Loc=%Root_Loc%\subs\FAST\subs
SET Farm_Loc=%Root_Loc%\src
SET Registry=%Root_Loc%\subs\FAST\bin\Registry_win32.exe

SET NWTC_Lib_Loc=%Subs_Loc%\NWTC_Library\source


SET ModuleName=%1

GOTO %ModuleName%

REM ----------------------------------------------------------------------------
REM ---------------- RUN THE REGISTRY TO AUTO-GENERATE FILES -------------------
REM ----------------------------------------------------------------------------
:FarmDriver
SET CURR_LOC=%Farm_Loc%\Driver
%REGISTRY% "%CURR_LOC%\FAST_Farm_Registry.txt" -I "%NWTC_Lib_Loc%" -noextrap -O "%CURR_LOC%"
GOTO checkError


:checkError
ECHO.
IF %ERRORLEVEL% NEQ 0 (
ECHO Error running FAST Registry for %ModuleName%.
) ELSE (
ECHO Registry for %ModuleName% completed.
REM COPY /Y "%ModuleName%_Types.f90"   "%CURR_LOC%"
rem IF /I "%ModuleName%"=="MAP" COPY /Y "%ModuleName%_Types.h" "%CURR_LOC%"
)




:end
REM ----------------------------------------------------------------------------
REM ------------------------- CLEAR MEMORY -------------------------------------
REM ----------------------------------------------------------------------------
ECHO. 


SET REGISTRY=

SET NWTC_Lib_Loc=
SET Root_Loc=
SET Subs_Loc=
SET Farm_Loc=

SET ModuleName=
SET CURR_LOC=
:Done
echo %lines%
set lines=