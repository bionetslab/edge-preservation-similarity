from distutils.core import setup
from setuptools import find_packages

setup(
    name='edge_preservation_similarity',
    version='0.1',
    author='Jana Kiederle',
    author_email='jana.kiederle@fau.de',
    packages=find_packages(),
    install_requires=['numpy', 'networkx', 'gurobipy']
    )