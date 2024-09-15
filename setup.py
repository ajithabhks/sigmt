from setuptools import setup

with open('requirements.txt') as f:
    requirements = f.read().splitlines()

setup(
    name='sigmt',
    version="2.0.0",
    packages=['sigmt', 'sigmt.core',
              'sigmt.gui', 'sigmt.utils',
              'sigmt.utils.edi', 'sigmt.utils.metronix',
              'sigmt.gui.project_related'],
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
        'Development Status :: 4 - Beta',
    ],
)
