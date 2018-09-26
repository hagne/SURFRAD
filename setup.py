import sys

required_verion = (3,)
if sys.version_info < required_verion:
    raise ValueError('SurfRadPy needs at least python {}! You are trying to install it under python {}'.format('.'.join(str(i) for i in required_verion), sys.version))

from setuptools import setup, find_packages

setup(
    name="SurfRadPy",
    version="0.1.1", #setting this caused a huge hadeache ... basically the script wasn't found when the version was set
    packages=find_packages(),
    author="Hagen Telg",
    author_email="hagen@hagnet.net",
    description="...",
    license="MIT",
    url="https://github.com/hagne/SURFRAD",
    install_requires=['pandas', 'numpy', 'xarray'],
    scripts=['scripts/qcrad2ncei'],
    # entry_points = {'console_scripts': ['qcrad2ncei=SurfRadPy.NCEI:qcrad2ncei'],},
    package_data={'': ['*.cdl']},
    include_package_data=True,
    zip_safe=False
)
