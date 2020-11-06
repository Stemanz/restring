# -*- coding: utf-8 -*-

from setuptools import setup, find_packages

with open("README.md", "r") as f:
    readme = f.read()

setup(
    name="restring",
    version="0.00.0",
    description="Functional enrichment terms aggregator.",
    long_description=readme,
    long_description_content_type="text/markdown",
    author="Stefano Manzini",
    author_email="stefano.manzini@gmail.com",
    url="https://github.com/Stemanz/restring",
    download_url="https://github.com/Stemanz/restring/archive/master.zip",
    license="GPL-3.0",
    packages=find_packages(exclude=("sample_data")),
    keywords = ['String', 'functional enrichment', 'GO', "David"],
    install_requires=[
          'matplotlib',
          'seaborn'
          'pandas',
      ],
)
