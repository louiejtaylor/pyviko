from distutils.core import setup

setup(name='pyviko',
	version='1.0.1.1',
	description='Design knockout viruses in Python',
	author='LJ Taylor',
	author_email='l'+'taylor'+str(3+4)+'@'+'tu'+'lane.edu',
	url='https://github.com/louiejtaylor/pyViKO',
	packages=['pyviko'],
	license='MIT License',
	classifiers = ['Intended Audience :: Science/Research',
				'Environment :: Console',
				'Environment :: Web Environment',
				'License :: OSI Approved :: MIT License',
				'Programming Language :: Python :: 2.7',
				'Topic :: Scientific/Engineering :: Bio-Informatics',],
	py_modules = ['pyviko.core', 'pyviko.mutation', 'pyviko.restriction']
     )
