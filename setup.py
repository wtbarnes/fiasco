import os
from distutils.core import setup

# create home directory
if not os.path.isdir(os.path.join(os.environ['HOME'], '.fiasco')):
    os.mkdir(os.path.join(os.environ['HOME'], '.fiasco'))

setup(
    name='fiasco',
    version='0.1dev',
    author='Will Barnes',
    url='https://github.com/wtbarnes/fiasco',
    #package_data={'fiasco':[]},
    packages=['fiasco','fiasco.io','fiasco.io.sources']
)