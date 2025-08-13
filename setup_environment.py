#!/usr/bin/env python3
"""
Setup script for PyMOCAT-MC environment
Creates a clean virtual environment with all dependencies
"""

import sys
import subprocess
import os
import shutil
from pathlib import Path

def run_command(cmd, description):
    """Run a command and handle errors"""
    print(f"Installing {description}...")
    try:
        result = subprocess.run(cmd, shell=True, capture_output=True, text=True, timeout=300)
        if result.returncode == 0:
            print(f"SUCCESS: {description} completed successfully")
            return True
        else:
            print(f"ERROR {description} failed:")
            print(f"  Error: {result.stderr}")
            return False
    except subprocess.TimeoutExpired:
        print(f"ERROR {description} timed out")
        return False
    except Exception as e:
        print(f"ERROR {description} error: {e}")
        return False

def main():
    print("Setup: PyMOCAT-MC Environment Setup")
    print("="*50)
    
    # Check Python version
    print(f"Info: Python version: {sys.version}")
    if sys.version_info < (3, 8):
        print("ERROR: Python 3.8 or higher is required")
        sys.exit(1)
    
    # Get repository root
    repo_root = Path(__file__).parent
    print(f"Directory: Repository root: {repo_root}")
    
    # Remove old virtual environment if it exists
    old_venvs = ['venv', 'test_env', 'verify_env']
    for venv_name in old_venvs:
        venv_path = repo_root / venv_name
        if venv_path.exists():
            print(f"Removing old virtual environment: {venv_name}")
            shutil.rmtree(venv_path)
    
    # Create new virtual environment
    venv_path = repo_root / "pymocat_env"
    if not run_command(f'python3 -m venv "{venv_path}"', "Creating virtual environment"):
        print("ERROR: Failed to create virtual environment")
        sys.exit(1)
    
    # Determine activation script path
    if sys.platform == "win32":
        activate_script = venv_path / "Scripts" / "activate"
        python_exe = venv_path / "Scripts" / "python.exe"
    else:
        activate_script = venv_path / "bin" / "activate"
        python_exe = venv_path / "bin" / "python"
    
    # Install dependencies
    requirements_file = repo_root / "requirements.txt"
    if requirements_file.exists():
        cmd = f'"{python_exe}" -m pip install --upgrade pip'
        if not run_command(cmd, "Upgrading pip"):
            print("WARNING: pip upgrade failed, continuing...")
        
        cmd = f'"{python_exe}" -m pip install -r "{requirements_file}"'
        if not run_command(cmd, "Installing dependencies"):
            print("ERROR: Failed to install dependencies")
            sys.exit(1)
    else:
        print("WARNING: No requirements.txt found, installing basic dependencies...")
        basic_deps = ["numpy", "scipy", "pandas", "matplotlib"]
        for dep in basic_deps:
            cmd = f'"{python_exe}" -m pip install {dep}'
            if not run_command(cmd, f"Installing {dep}"):
                print(f"WARNING: Failed to install {dep}")
    
    # Test the environment
    print("\nTesting environment...")
    test_cmd = f'"{python_exe}" -c "import numpy, scipy, matplotlib; print(\\"Environment test passed!\\")"'
    if run_command(test_cmd, "Testing basic imports"):
        print("SUCCESS: Environment setup completed successfully!")
    else:
        print("ERROR: Environment test failed")
        sys.exit(1)
    
    # Provide usage instructions
    print("\nInfo: Usage Instructions:")
    print("="*50)
    print("To activate the environment:")
    if sys.platform == "win32":
        print(f"  {venv_path}\\Scripts\\activate")
    else:
        print(f"  source {venv_path}/bin/activate")
    
    print("\nTo run tests:")
    print("  python tests/run_tests.py")
    print("  # or")
    print("  ./tests/run_tests.sh")
    
    print("\nTo deactivate:")
    print("  deactivate")
    
    print(f"\nSUCCESS: Setup complete! Virtual environment created at: {venv_path}")

if __name__ == "__main__":
    main()