#/urs/bin/env python

from setuptools import setup, find_packages
# import versioneer

__author__ = "Justine Debelius"
__copyright__ = "2019-- Justine Debelius"
__credits__ = ["Justine Debelius", "Noam Shental"]
__license__ = "BSD-3-Clause"
__version__ = "0.0.1-dev"
__maintainer__ = "Justine Debelius"
__email__ = "j.debelius@gmail.com"


setup(name='q2-sidle',
      version=__version__,
      description="A qiime 2 implemention of Short MUliple Regions Framework",
      license=__license__,
      author=__maintainer__,
      author_email=__email__,
      maintainer=__maintainer__,
      maintainer_email=__email__,
      # packages=find_packages(),       
      packages=['q2_sidle', 'q2_sidle.tests'],
      # modules=['py_smurf.cli'],
      # url=,
       entry_points={
        'qiime2.plugins':
        ['q2-sidle=q2_sidle.plugin_setup:plugin']
      },
      # entry_points="""
      #     [console_scripts]
      #     sidle=sidle.cli:cli
      # """,
      install_requires=['biom-format >= 2.1.6',
                        'sparse >= 0.8.0',
                        'pandas == 0.25.0',
                        'dask >= 2.0',
                        'numba == 0.48',
                        'regex',
                        ],
      )
