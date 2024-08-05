from pathlib import Path

from setuptools import  setup, find_packages

from datetime import datetime

now = datetime.now()

__version__ = f"{now.year % 100}.{now.month+1}.0"

setup(
    name="spartaabc",
    version=__version__,
    author="Elya Wygoda",
    author_email="elya.wygoda@gmail.com",
    description="Simulation based inference on indel parameters",
    long_description=open("README.md", 'r').read(),
    long_description_content_type='text/markdown',
    zip_safe=False,
    packages=find_packages(exclude=["test"]),
    install_requires=["msastats", "msasim"],
    python_requires=">=3.7",

)