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
REM --------      THESE PATHS SHOULD WORK "AS IS", HOWEVER ---------------------
REM -- USERS MAY EDIT THESE PATHS TO POINT TO FOLDERS ON THEIR LOCAL MACHINES. -
REM -- NOTE: do not use quotation marks around the path names!!!! --------------
REM -- IT is also NOT RECOMMENDED to use spaces in the path names!!!! ----------
REM ----------------------------------------------------------------------------
REM ----------------------------------------------------------------------------

REM Note that Root_Loc is relative to the projects that **use** this batch file. Thus, it must 
REM reflect one directory lower than the location of this batch file.
SET Root_Loc=..\..\..


SET Farm_Loc=%Root_Loc%\src

SET Driver_Loc=%Farm_Loc%\Driver
SET Wake_Loc=%Farm_Loc%\WakeDynamics
SET Wrapper_Loc=%Farm_Loc%\FASTWrapper
SET AmbWind_Loc=

:: Get all of the paths we'd normally use in FAST, but make them relative to this %FAST_Loc% instead 
:: note that Root_Loc and FAST_Loc get overwritten!
SET FAST_Loc=%Root_Loc%\subs\FAST\
SET SetFASTPaths=%FAST_Loc%\Compiling\VisualStudio\FASTlib\RunRegistry.bat
call %SetFASTPaths% PathsOnly %FAST_Loc%

SET ModuleName=%1

GOTO %ModuleName%

REM ----------------------------------------------------------------------------
REM ---------------- RUN THE REGISTRY TO AUTO-GENERATE FILES -------------------
REM ----------------------------------------------------------------------------
:FarmDriver
SET CURR_LOC=%Driver_Loc%
%REGISTRY% "%CURR_LOC%\FAST_Farm_Registry.txt" -I %Driver_Loc% -I %Wake_Loc% -I %Wrapper_Loc% %ALL_FAST_INCLUDES% -noextrap -O "%CURR_LOC%"
GOTO checkError

:FASTWrapper
SET CURR_LOC=%Wrapper_Loc%
%REGISTRY% "%CURR_LOC%\FASTWrapper_Registry.txt" -I %NWTC_Lib_Loc%  %ALL_FAST_INCLUDES% -noextrap -O "%CURR_LOC%"
GOTO checkError

:WakeDynamics
SET CURR_LOC=%Wake_Loc%
%REGISTRY% "%CURR_LOC%\WakeDynamics_Registry.txt" -I %NWTC_Lib_Loc% -noextrap -O "%CURR_LOC%"
GOTO checkError

:WakeDynamics_Driver
SET CURR_LOC=%Wake_Loc%\driver
%REGISTRY% "%CURR_LOC%\WakeDynamics_Driver_Registry.txt" -I %NWTC_Lib_Loc% -I %Wake_Loc% -noextrap -O %CURR_LOC% 



:checkError
ECHO.
IF %ERRORLEVEL% NEQ 0 (
ECHO Error running FAST Registry for %ModuleName%.
) ELSE (
ECHO Registry for %ModuleName% completed.
)


:end
REM ----------------------------------------------------------------------------
REM ------------------------- CLEAR MEMORY -------------------------------------
REM ----------------------------------------------------------------------------
call %SetFASTPaths% end

SET Farm_Loc=

SET Driver_Loc=
SET Wake_Loc=
SET AmbWind_Loc=

SET ModuleName=
SET CURR_LOC=
SET SetFASTPaths=

:Done
echo %lines%
set lines=