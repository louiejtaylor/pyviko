from distutils.core import setup

setup(name='pyviko',
      version='1.0',
      description='Python Viral Knockouts',
      author='Louis Taylor',
      author_email='louist@upenn.edu',
      url='',
      py_modules = ['pyviko.core', 'pyviko.mutation', 'pyviko.restriction'],
      data_files=[('tests', ['tests/core_test.py'])]
     )