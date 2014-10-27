import os
from setuptools import setup, find_packages

setup(name='ncpipe',
      version='0.0.1',
      description='ncpipe is a Python module to analyze ncRNA NGS data',
      author='Patrick Lombard',
      author_email='ptk.lmb55@gmail.com',
      packages=['ncpipe'],
      package_data={"ncpipe":['data/*']},
      scripts=['scripts/ncalign.py', 'scripts/ncanno.py', 'scripts/ncount.py', 'scripts/ncviz.py'],
      install_requires=['pysam', 'pybedtools'],
      license='GPLv3',
      platforms='any',
      classifiers=[
          'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
          'Development Status :: 3 - Alpha',
          'Programming Language :: Python :: 2.7',
          'Environment :: Console',
      ],
      long_description="""

ncpipe is a Python module to analyze ncRNA NGS data

 Contact
=============

If you have any questions or comments about ncpipe, please feel free to contact me via
eMail: ptk.lmb55@gmail.com

""",
    )