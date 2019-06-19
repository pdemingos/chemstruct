# Copyright 2019 Pedro G. Demingos

"""Simple class for Comma-Separated Values files."""

from files.main import *


class Csv(File):
    """Class for csv (comma-separated values) files.
    Internally used for some tables. Not very efficient."""

    def __init__(self, absolute_path=None):
        super().__init__(absolute_path)
        self.tags.add("csv")

        self.data = []  # should be a list of lists (the columns)
        self.headers = []  # should be a list of strings (name of columns)
        self.max_len = 0

        if absolute_path is not None:
            if absolute_path.endswith(".csv"):
                self.read_csv()
            else:
                self.read_generic()

    def read_csv(self):
        pass

    def read_generic(self):
        for line in self.content:
            for (index, number) in enumerate(line.split()):
                try:
                    self.data[index].append(number)
                except IndexError:
                    self.data.append([])
                    self.data[index].append(number)

    def add_array(self, array: list):
        """Adds a column to the table.
        Must be given in the desired order."""
        if not isinstance(array, list):
            raise TypeError("array must be list, got {}".format(type(array)))
        self.data.append(array)
        if len(array) > self.max_len:
            self.max_len = len(array)

    def add_header(self, header: str):
        """Adds a header to the next (last) column.
        Must be given in the desired order."""
        try:
            header = str(header)
        except TypeError:
            raise TypeError("header should be a string, got {}".format(
                type(header)))
        else:
            self.headers.append(header)

    def fill_columns(self):
        """Fills columns with 'NaN' so they all have the same
        length i.e. number of lines."""
        for i in range(self.max_len):
            for d in self.data:
                try:
                    _ = d[i]
                except IndexError:
                    d.append("NaN")
                else:
                    continue

    def write_csv(self, path: str):
        """Writes a csv file with all info present in the object."""
        if not self.data:
            raise ValueError("no data in Csv object")
        if any(len(d) < self.max_len for d in self.data):
            print("WARNING: filling columns with 'NaN'...")
            self.fill_columns()
        with open(path, "w") as F:
            if self.headers:
                F.write(",".join(self.headers) + "\n")
            for i in range(self.max_len):
                line = ",".join(str(d[i]) for d in self.data) + "\n"
                F.write(line)

    def sort_by(self, index: int):
        """Re-orders rows so the values in the column whose index is given
        are in ascending order."""
        value_to_index = dict()
        for (i, value) in enumerate(self.data[index]):
            value_to_index[value] = i
        sorted_values = sorted(self.data[index])
        new_data = []
        for _ in self.data:
            new_data.append([])
        for value in sorted_values:
            for (column_i, column) in enumerate(new_data):
                if column_i == index:
                    column.append(value)
                else:
                    column.append(self.data[column_i][value_to_index[value]])
        self.data = new_data
