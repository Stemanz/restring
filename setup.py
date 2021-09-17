# -*- coding: utf-8 -*-
import os
from setuptools import setup, find_packages
import sys

# Python compatibility check
major, minor = sys.version_info[:2]

if major < 3:
    print("You are trying to use Python 2.")
    print("Python 3.8 or newer is required to install and use reString.")
    print("Please update your Python version, or create a suitable virtual environment.")
    print("\nTip: try 'pip3 install restring' and see if this fixes the issue.")
    sys.exit()
elif major >= 3 and minor < 8:
    print(f"You using Python {major}.{minor}.")
    print("Python 3.8 or newer is required to install and use reString.")
    print("Please update your Python version, or create a suitable virtual environment.")
    sys.exit()

with open("README.md", "r") as f:
    readme = f.read()

# thanks https://github.com/mindflayer/python-mocket/blob/master/setup.py
def get_version(packagedir):
    init_path = os.path.join(packagedir, "__init__.py")
    with open(init_path, "r") as pyinit:
        for line in pyinit:
            if line.startswith("__version__"):
                return line.split()[-1].strip().strip('"') #middle .strip() just in case
    
setup(
    name="restring",
    # this comes from: https://stackoverflow.com/questions/19534896/enforcing-python-version-in-setup-py
    python_requires=">3.8.0",
    version=get_version("restring"),
    description="Functional enrichment terms aggregator.",
    long_description=readme,
    long_description_content_type="text/markdown",
    author="Stefano Manzini",
    author_email="stefano.manzini@gmail.com",
    url="https://github.com/Stemanz/restring",
    download_url="https://github.com/Stemanz/restring/archive/master.zip",
    license="GPL-3.0",
    packages=find_packages(exclude=("data", "images", "sample_data", "sample_tables")),
    keywords=['String', 'functional enrichment', 'GO', "David", "KEGG", "pathways"],
    install_requires=[
        'matplotlib',
        'seaborn',
        'pandas',
        'requests',
      ],
    entry_points={
        'console_scripts': ['restring-gui=restring.restring:restring_gui'],
    },
    include_package_data=True, # processes MANIFEST.in
)
