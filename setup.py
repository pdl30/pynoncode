import os
from setuptools import setup, find_packages

setup(name='pynoncode',
      version='0.0.1',
      description='pynoncode is a Python module to analyze ncRNA NGS data',
      author='Patrick Lombard',
      author_email='ptk.lmb55@gmail.com',
      packages=['pynoncode'],
      package_data={"pynoncode":['data/*']},
      scripts=['scripts/pynon_align.py', 'scripts/pynon_anno.py', 'scripts/pynon_count.py', 'scripts/pynon_viz.py'],
      install_requires=['numpy', 
        'scipy'
        'cython', 
        'pysam', 
        'pybedtools', 
        'HTSeq'],
      license='GPLv3',
      platforms='any',
      classifiers=[
          'License :: OSI Approved :: GNU General Public License v2 (GPLv2)',
          'Development Status :: 3 - Alpha',
          'Programming Language :: Python :: 2.7',
          'Environment :: Console',
      ],
      long_description="""

pynoncode is a Python module to analyze ncRNA NGS data

 Contact
=============

If you have any questions or comments about pynoncode, please feel free to contact me via
eMail: ptk.lmb55@gmail.com

""",
    )
