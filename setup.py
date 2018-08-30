from setuptools import setup
from setuptools import find_packages

setup(name='matsim',
      version='0.0.1',
      decription='materials simulation toolkits',
      author='Chi Chen',
      author_email='chen08013@gmail.com',
      url='https://chc273.github.io',
      download_url='...',
      license='MIT',
      install_requires=['pymatgen', 'numpy', 'phonopy'],
      extras_require={
          'model_saving': ['json', 'h5py'],},
      package_data = {'matsim':['README.md']},
      package=find_packages())
