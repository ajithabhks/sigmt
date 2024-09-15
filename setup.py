from setuptools import setup, find_packages

from sigmt.__version__ import __version__

with open('requirements.txt') as f:
    requirements = f.read().splitlines()

setup(
    name='sigmt',
    version=__version__,
    packages=find_packages(),
    include_package_data=True,
    install_requires=requirements,
    entry_points={
        'console_scripts': [
            'sigmt=sigmt.main:main',
        ],
    },
    python_requires='>=3.8',
    # Metadata
    author='Ajithabh K.S.',
    author_email='ajithabhks@gmail.com',
    description='A tool for magnetotelluric time series processing',
    long_description=open('README.md').read(),
    long_description_content_type='text/markdown',
    url='https://github.com/ajithabhks/sigmt',
    license='GPL',
    classifiers=[
        'License :: OSI Approved :: GNU General Public License (GPL)',
        'Programming Language :: Python :: 3',
    ],
)