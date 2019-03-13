""" :module TestUtilities: Hosting various utilities useful for testing. """

def _remove_test_files(files):
    """ """
    """ Remove all files and directories listed.

    :param files: Files and directories to remove.
    :type files: list

    """

    # Loop over files
    for f in files:
        # Check if file.
        if os.is_file(f):
            os.remove(f)
        # Check if dir.
        elif os.is_dir(f):
            shutil.rmtree(f)


