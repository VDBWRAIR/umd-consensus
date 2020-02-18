from glob import glob
from setuptools import setup, find_packages
import consensus 
# Run setuptools setup

setup(
    name = "consensus",
    packages = find_packages(),
    scripts = glob('bin/*'),
    entry_points = {
        'console_scripts':
            'consensus_generator = consensus.consensus:main'

        } )

