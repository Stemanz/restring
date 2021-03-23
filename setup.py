# -*- coding: utf-8 -*-
import os
from setuptools import setup, find_packages

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
    keywords = ['String', 'functional enrichment', 'GO', "David", "KEGG", "pathways"],
    install_requires=[
          'matplotlib',
          'seaborn',
          'pandas',
      ],
    scripts=['bin/restring-gui'],
    include_package_data=True, # processes MANIFEST.in
)
