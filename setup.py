import sys

required_verion = (3,)
if sys.version_info < required_verion:
    raise ValueError('SurfRadPy needs at least python {}! You are trying to install it under python {}'.format('.'.join(str(i) for i in required_verion), sys.version))

# from setuptools import setup
from distutils.core import setup
setup(
    name="SurfRadPy",
    version="0.1",
    packages=['SurfRadPy'],
    author="Hagen Telg",
    author_email="hagen@hagnet.net",
    description="...",
    license="MIT",
    # keywords="matplotlib",
    url="https://github.com/hagne/SURFRAD",
    install_requires=['pandas', 'numpy', 'xarray'],
    scripts=['scripts/qcrad2ncei'],
    data_files=[('data', ['data/SURFRAD_QCrad_metadata.cdl'])]
    # extras_require={'plotting': ['matplotlib'],
    #                 'testing': ['scipy']},
    # test_suite='nose.collector',
    # tests_require=['nose'],
)
