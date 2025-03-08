"""Download example scripts to a local folder, `cocoex-example` by default or a given argument::

    python -m cocoex.download_example

"""
import os
import sys

local_folder = 'cocoex-example'

remote_folder = 'https://raw.githubusercontent.com/numbbo/coco-experiment/main/build/python/example/'

file_names = ['mkl_bugfix.py', 
        'example_experiment_simple.py',
        'example_experiment_complete.py',
        'example_experiment_non_anytime.py',
        'example_experiment_parameter_sweep.py',
        ]

if __name__ == '__main__':
    try:
        import urllib.request
    except ImportError:
        print('Install urllib (`pip install urllib`) to download example')
        raise
    folder = sys.argv[1] if len(sys.argv) > 1 else local_folder        
    if os.path.exists(folder):
        raise ValueError('local target folder {0} already exists. Remove the folder or\n'
                         'pass a (new) folder name as argument to `cocoex.download_example`'
                         .format(folder))
    os.mkdir(folder)
    for file_name in file_names:
        with urllib.request.urlopen(remote_folder + file_name) as f:
            d = f.read().decode('utf-8')
        with open(folder + '/' + file_name, 'w') as f:
            f.write(d)
    print('  example scripts have been downloaded to {0}'.format(folder))
