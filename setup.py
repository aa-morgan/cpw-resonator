from setuptools import setup

def readme():
    with open('README.md') as f:
        return f.read()

setup(name='rydberg',
      version='0.0.1',
      description='Package with tools useful for designing co-planar waveguide (CPW) resonators.',
      url='',
      author='Alex Morgan',
      author_email='axm108@gmail.com',
      license='BSD 3-clause',
      packages=['cpwprop'],
      include_package_data=True,
      zip_safe=False)
