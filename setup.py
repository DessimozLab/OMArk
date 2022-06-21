from setuptools import setup, find_packages


name = 'omark'
__version__ = None
with open('{:s}/__init__.py'.format(name), 'rt') as fp:
    for line in fp:
        if line.startswith('__version__'):
            exec(line.rstrip())

requirements = ['biopython', 'ete3', 'omamer>=0.2.2', 'matplotlib']

desc = 'OMArk - Proteome quality assesment based on OMAmer placements'

setup(
    name=name,
    version=__version__,
    author='Yannis Nevers',
    email='yannis.nevers@unil.ch',
    url='https://github.com/YanNevers/omark',
    description=desc,
    packages=find_packages(),
    install_requires=requirements,
    python_requires=">=3.9",
    scripts=['bin/omark'])

