#!/bin/bash
# MetaboGraph Launcher for Linux/Mac
# This script launches the MetaboGraph GUI application

echo "Starting MetaboGraph..."
echo ""

# Get the directory where this script is located
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"

# Run the GUI
python3 "$SCRIPT_DIR/main_script/run_gui.py"

# Check exit code
if [ $? -ne 0 ]; then
    echo ""
    echo "An error occurred. Press Enter to exit..."
    read
fi
