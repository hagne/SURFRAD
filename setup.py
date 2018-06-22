import sys

required_verion = (3,)
if sys.version_info < required_verion:
    raise ValueError('SurfRadPy needs at least python {}! You are trying to install it under python {}'.format('.'.join(str(i) for i in required_verion), sys.version))

# import ez_setup
# ez_setup.use_setuptools()

from setuptools import setup
# from distutils.core import setup
setup(
    name="SurfRadPy",
    version="0.1",
    packages=['SurfRadPy'],
    author="Hagen Telg",
    author_email="hagen@hagnet.net",
    description="...",
    license="MIT",
    # keywords="matplotlib",
    url="https://github.com/hagne/SURFRAD", install_requires=['pandas', 'numpy', 'xarray']
    # install_requires=['numpy','pandas'],
    # extras_require={'plotting': ['matplotlib'],
    #                 'testing': ['scipy']},
    # test_suite='nose.collector',
    # tests_require=['nose'],
)