#/urs/bin/env python

from setuptools import setup, find_packages
# import versioneer

__author__ = "Justine Debelius"
__copyright__ = "2019-- Justine Debelius"
__credits__ = ["Justine Debelius", "Noam Shental"]
__license__ = "BSD-3-Clause"
__version__ = "2021.2-dev"
__maintainer__ = "Justine Debelius"
__email__ = "j.debelius@gmail.com"
__url__ = "https://q2-sidle.readthedocs.io/en/latest/index.html"

setup(name='q2-sidle',
      version=__version__,
      description="A qiime 2 implemention of Short MUliple Regions Framework",
      license=__license__,
      author=__maintainer__,
      author_email=__email__,
      maintainer=__maintainer__,
      maintainer_email=__email__,
      url=__url__,
      # packages=find_packages(),       
      packages=['q2_sidle', 'q2_sidle.tests'],
      package_data={
            'q2_sidle': ['citations.bib'],
       },
       entry_points={
        'qiime2.plugins':
        ['q2-sidle=q2_sidle.plugin_setup:plugin']
      },
      install_requires=['biom-format >= 2.1.6',
                        'pandas > 1.0',
                        'dask >= 2.0',
                        'regex',
                        ],
      )
