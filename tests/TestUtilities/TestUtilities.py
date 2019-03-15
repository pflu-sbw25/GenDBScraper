""" :module TestUtilities: Hosting various utilities useful for testing. """

import os, shutil

def _remove_test_files(files):
    """ """
    """ Remove all files and directories listed.

    :param files: Files and directories to remove.
    :type files: list

    """

    # Loop over files
    for f in files:
        # Check if file.
        if os.path.isfile(f):
            os.remove(f)
        # Check if dir.
        elif os.path.isdir(f):
            shutil.rmtree(f)


