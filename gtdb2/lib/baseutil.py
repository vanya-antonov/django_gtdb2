import fileinput


def get_all_ids(fn):
    """Returns a list  of integers from 'fn'.
    
    Arguments:
        - fn - a string representing a input file name. If None or equal to '-' - read
          from STDIN; fn is treated as a single ID if it is a digit string.
    """
    if fn is None:
        fn = '-'
    
    if fn.isdigit():
        return [int(fn)]   # just a single ID
    
    all_ids = []
    for line in fileinput.input(files=fn):
        line = line.rstrip()
        if line.isspace():
            continue
        if not line.isdigit():
            raise Exception("Wrong ID '%s'!" % line)
        all_ids.append(int(line))
    
    return all_ids
