"""Download example scripts from GitHub to a local folder, by default `cocoex-example`
in the current folder.

Usage::

    python -m cocoex.download_example [folder_name]

"""
import os
import sys

default_local_folder = 'cocoex-example'

remote_folder = 'https://raw.githubusercontent.com/numbbo/coco-experiment/main/build/python/example/'

file_names = ['mkl_bugfix.py', 
        'example_experiment_simple.py',
        'example_experiment_complete.py',
        'example_experiment_non_anytime.py',
        'example_experiment_parameter_sweep.py',
        ]

def main(folder=None):
    """download example scripts from GitHub into `folder`, by default `./cocoex-example`"""
    try:
        import urllib.request
    except ImportError:
        print('Install urllib (`pip install urllib`) to download example')
        raise
    if folder is None:
        folder = default_local_folder
    if os.path.exists(folder):
        raise ValueError('local target folder {0} already exists. Remove the folder\n'
                         'or pass another folder name as argument to `cocoex.download_example`'
                         .format(folder))
    os.mkdir(folder)
    for file_name in file_names:
        with urllib.request.urlopen(remote_folder + file_name) as f:
            d = f.read().decode('utf-8')
        with open(folder + '/' + file_name, 'w') as f:
            f.write(d)
    print('  example scripts have been downloaded to {0}'.format(folder))

if __name__ == '__main__':
    main(sys.argv[1] if len(sys.argv) > 1 else None)