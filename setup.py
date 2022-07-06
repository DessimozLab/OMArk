from setuptools import setup, find_packages


name = 'omark'
__version__ = None
with open('{:s}/__init__.py'.format(name), 'rt') as fp:
    for line in fp:
        if line.startswith('__version__'):
            exec(line.rstrip())

with open("README.md", "rt") as fh:
    readme = fh.read()

requirements = ['biopython', 'ete3', 'omamer>=0.2.2', 'matplotlib', 'jinja2']

desc = 'OMArk - Proteome quality assesment based on OMAmer placements'

setup(
    name=name,
    version=__version__,
    author='Yannis Nevers',
    email='yannis.nevers@unil.ch',
    url='https://github.com/DessimozLab/omark',
    description=desc,
    long_description=readme,
    long_description_content_type="text/markdown",
    packages=find_packages(),
    package_data={'omark':['assets/*.txt']},
    include_package_data=True,
    install_requires=requirements,
    python_requires=">=3.9",
    scripts=['bin/omark'])

