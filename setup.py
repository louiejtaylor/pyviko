from distutils.core import setup

setup(name='pyviko',
      version='1.0',
      description='Python Viral Knockouts',
      author='Louis Taylor',
      author_email='l'+'taylor'+str(3+4)+'@'+'tu'+'lane.edu',
      url='',
      py_modules = ['pyviko.core', 'pyviko.mutation', 'pyviko.restriction'],
      data_files=[('tests', ['tests/core_test.py'])]
     )
