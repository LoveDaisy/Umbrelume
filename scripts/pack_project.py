#!/usr/bin/env python3
"""
pack_project.py - Package lunar eclipse color simulation project

Creates a zip file containing:
- All Python files from python/ folder
- Complete docs/ folder  
- CSV data files from data/ folder

Usage:
    python scripts/pack_project.py [output_name]

If no output name is provided, defaults to 'lunar_eclipse_color.zip'
"""

import os
import sys
import zipfile
import glob
from pathlib import Path

def pack_project(output_name="lunar_eclipse_color.zip"):
    """Pack the project into a zip file."""
    
    # Get project root directory (parent of scripts folder)
    script_dir = Path(__file__).parent
    project_root = script_dir.parent
    
    # Change to project root for consistent paths in zip
    os.chdir(project_root)
    
    # Remove existing zip file if it exists
    if os.path.exists(output_name):
        os.remove(output_name)
        print(f"Removed existing {output_name}")
    
    # Files to include
    files_to_pack = []
    
    # Add Python files
    python_files = glob.glob("python/*.py")
    files_to_pack.extend(python_files)
    
    # Add docs folder (all files)
    for root, dirs, files in os.walk("docs"):
        for file in files:
            files_to_pack.append(os.path.join(root, file))
    
    # Add CSV files from data folder
    csv_files = glob.glob("data/*.csv")
    files_to_pack.extend(csv_files)
    
    # Create zip file
    with zipfile.ZipFile(output_name, 'w', zipfile.ZIP_DEFLATED) as zf:
        for file_path in sorted(files_to_pack):
            if os.path.exists(file_path):
                zf.write(file_path)
                print(f"  adding: {file_path}")
    
    # Show results
    zip_size = os.path.getsize(output_name)
    print(f"\nCreated {output_name} ({zip_size/1024:.0f}K)")
    print(f"Contains {len(files_to_pack)} files:")
    print(f"  - {len(python_files)} Python files")
    print(f"  - {len(csv_files)} CSV data files")
    print(f"  - {len([f for f in files_to_pack if f.startswith('docs')])} documentation files")

if __name__ == "__main__":
    output_name = sys.argv[1] if len(sys.argv) > 1 else "lunar_eclipse_color.zip"
    pack_project(output_name)