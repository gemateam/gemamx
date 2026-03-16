import os
import subprocess
import sys

def setup():
    print("Setting up Python virtual environment...")
    if not os.path.exists("venv"):
        subprocess.run([sys.executable, "-m", "venv", "venv"])
    
    pip_path = os.path.join("venv", "Scripts", "pip.exe") if os.name == 'nt' else os.path.join("venv", "bin", "pip")
    
    print("Installing dependencies from requirements.txt...")
    subprocess.run([pip_path, "install", "-r", "requirements.txt"])
    
    print("\nEnvironment setup complete.")
    print("To activate the environment, run:")
    print("  .\\venv\\Scripts\\activate" if os.name == 'nt' else "  source venv/bin/activate")

if __name__ == "__main__":
    setup()
