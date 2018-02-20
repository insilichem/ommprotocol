if not exist "%PREFIX%\Menu" mkdir "%PREFIX%\Menu"
copy "%RECIPE_DIR%\menu-windows.json" "%PREFIX%\Menu"
copy "%RECIPE_DIR%\app.ico" "%PREFIX%\Menu"

"%PYTHON%" setup.py install
if errorlevel 1 exit 1