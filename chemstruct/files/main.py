# Copyright 2019 Pedro G. Demingos

"""Basic functions and class for files."""

import os


def find_between(line: str, char1: str, char2: str, nth_time=1):
    """
    Goes through the given line, looking for what's between char1 and char2.
    Returns substring found. If char1...char2 is present more than once in
    line, the nth_time parameter can be used to specify which occurrence should
    be returned. The char1 and char2 are NOT returned along with the substring.

    Parameter
    ---------
    line : str
        Line where substring shall be looked for.
    char1 : str
        Char (single char!) at the beginning of the desired substring.
    char2 : str
        Char (single char!) at the end of the desired substring.
    nth_time : int, optional
        Which occurrence of the char1...char2 pattern must be returned.
        if no number is given, nth_time is set to 1.

    Returns
    -------
    between : str
        Substring found between char1 and char2.
        Does NOT include char1 and char2.

    """

    if line.count(char1) < nth_time or line.count(char2) < nth_time:
        raise ValueError("bad arguments given to find_between()")

    between = ""
    this_time = 1
    line = iter(line)
    char = next(line)

    while True:
        while char != char1:
            char = next(line)
        char = next(line)
        while char != char2:
            between += char
            char = next(line)
        if this_time == nth_time:
            return between
        else:
            this_time += 1
            between = ""


def clear_end(line: str, chars: list):
    """
    Clears line's end from unwanted chars.

    Parameters
    ----------
    line : str
        Line to be cleared.
    chars : list of chars
        Unwanted chars.

    Returns
    -------
    line : str
        Given line, cleared from unwanted chars.

    """
    while any(line.endswith(char) for char in chars):
        line = line[:-1]
    return line


def clear_start(line: str, chars: list):
    """
    Clears line's beginning from unwanted chars.

    Parameters
    ----------
    line : str
        Line to be cleared.
    chars : list of chars
        Unwanted chars.

    Returns
    -------
    line : str
        Given line, cleared from unwanted chars.

    """
    while any(line.startswith(char) for char in chars):
        line = line[1:]
    return line


class File:
    """Base class for files."""

    def __init__(self, absolute_path=None):
        """

        Parameters
        ----------
        absolute_path : str
            Absolute path to the file.

        """
        self._absolute_path = None
        self.dir_path = None
        if absolute_path is not None:
            self.absolute_path = absolute_path
            self.dir_path = "/".join(absolute_path.split("/")[:-1])
        self._content = None
        self._complete_content = None
        if self._absolute_path is not None:
            self.read_file()
        self.tags = set()

    @property
    def content(self):
        return self._content

    def classify(self):
        # find out what type of file it is
        # instantiate new obj with right class
        # call new obj's classify() method
        pass

    @property
    def absolute_path(self):
        return self._absolute_path

    @absolute_path.setter
    def absolute_path(self, absolute_path: str):
        if not isinstance(absolute_path, str):
            raise TypeError("absolute_path must be str")
        if not os.path.isfile(absolute_path):
            raise NameError("absolute_path must be a path")
        self._absolute_path = absolute_path

    def read_file(self):
        """Reads file in self's absolute_path."""
        if self._absolute_path is None:
            raise NameError("absolute_path not defined")
        else:
            with open(self._absolute_path, "r") as F:
                self._content = F.readlines()

        self._complete_content = self._content
        self._remove_comments()

    def _remove_comments(self, comment_symbol="#"):
        # CARE: '#' may be relevant in some files.
        self._content = []
        for line in self._complete_content:
            if comment_symbol in line:
                clear_line = ""
                for char in line:
                    if char == comment_symbol:
                        break
                    else:
                        clear_line += char
                self._content.append(clear_line)
            else:
                self._content.append(line)

    def _check_read(self):
        if self._content is None:
            self.read_file()

    def print_file(self):
        """Prints lines in self's content.
        If self._content is None, calls self.read_file()."""
        self._check_read()
        for line in self._content:
            print(line, end='')

    def about(self):
        """Prints object's tags."""
        print(self.tags)

    def _find(self, string):
        """
        Goes through all lines of self.content, looking for given string.
        Returns list with the indexes of the lines where the given string
        was found.

        Parameters
        ----------
        string : str
            String to look for.

        Returns
        -------
        indexes : list or ints
            All indexes (in order) where the given string was found in
            self.content.

        """
        self._check_read()
        indexes = []
        for (index, line) in enumerate(self._content):
            if string in line:
                indexes.append(index)
        return indexes

    def simple_write(self, filename):
        with open(filename, "w") as F:
            for line in self.content:
                F.write(line)
