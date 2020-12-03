from setuptools import find_packages, setup

setup(
    name='wrfhydropy',
    version='0.0.18',
    packages=find_packages(),
    package_data={'wrfhydropy': ['core/data/*']},
    url='https://github.com/NCAR/wrf_hydro_py',
    license='MIT',
    install_requires=[
        'pandas>=1.0.3',
        'f90nml>=1.2',
        'netCDF4>=1.5.3',
        'numpy==1.18.1',
        'deepdiff==3.3.0',
        'pathlib==1.0.1',
        'xarray==0.14.1',
        'properscoring==0.1',
        'pytest>=5.4.1',
        'pytest-html>=3.0.0',
        'pytest-datadir-ng>=1.1.1',
        'pytest-lazy-fixture>=0.6.3',
        'boltons>=20.2.1',
        'bs4>=0.0.1',
        'requests>=2.23.0',
        'dask[bag]>=2.14.0',
        'spotpy>=1.5.14'
    ],
    author='James McCreight',
    author_email='jamesmcc@ucar.edu',
    description='Crude API for the WRF-Hydro model',
)
