import numpy as np
from AnyQt.QtCore import Qt

import Orange.data
import Orange.misc
from Orange.widgets.utils.sql import check_sql_input
from Orange.widgets.widget import OWWidget, Msg


class OWMultipleSequenceAlignment(OWWidget):
    name = "Multiple Sequence Alignment"
    description = "Compute a matrix of pairwise sequence alignment distance."
    icon = "icons/MultipleSeqAlignment.svg"

    inputs = [("Data", Orange.data.Table, "set_data")]
    outputs = [("Distances", Orange.misc.DistMatrix)]

    raw_output = Orange.misc.DistMatrix(np.array([]))

    class Error(OWWidget.Error):
        no_continuous_features = Msg("No continuous features")
        dense_metric_sparse_data = Msg("Selected metric does not support sparse data")
        empty_data = Msg("Empty data (shape = {})")
        too_few_observations = Msg("Too few observations for the number of dimensions")

    class Warning(OWWidget.Warning):
        ignoring_discrete = Msg("Ignoring discrete features")
        imputing_data = Msg("Imputing missing values")

    def __init__(self):
        super().__init__()

        self.data = None
        # No gui

    @check_sql_input
    def set_data(self, data):
        """
        Set the input data set from which to compute the distances
        """
        self.data = data
        self.commit()

    def commit(self):
        self.send("Distances", self.compute_alignment(self.data))

    def computeDistance(self, s, p):
        return -1.0

    def compute_alignment(self, data):
        #HACK
        n = len([row[0] for row in data])
        print([row[0].value for row in data])
        outdata = np.zeros([n, n])
        for i, row in enumerate(data):
            for j, rowCompare in enumerate(data):
                if i != j:
                    outdata[i, j] = self.computeDistance(str(row[0].value), str(rowCompare[0].value))



        """ Mock """
        #data = [range(0, 4), range(4, 8), range(8, 12), range(12, 16)]
        self.raw_output = Orange.misc.DistMatrix(np.array(outdata))
        # TODO align given sequences and return edit distance matrix
        return self.raw_output

"""
 Reads some data from a table file and sets it as
 the input of this widget.
 Then pipes the output of this widget into the
 distance matrix widget to show the results

 NOTE: if you wish to see anything the widget
    displays (currently nothing), then
    unncoment ow.show() line

 TODO: figure out, why it throws error after
       finishing,
       figure out how to properly pipe output
       without using raw_output
"""
if __name__ == "__main__":
    import sys
    from AnyQt.QtWidgets import QApplication
    from Orange.data import Table
    from Orange.data.domain import Domain, DiscreteVariable, ContinuousVariable
    from Orange.widgets.unsupervised.owdistancematrix import OWDistanceMatrix

    a = QApplication(sys.argv)
    ow = OWMultipleSequenceAlignment()

    #setup test data
    domain = Domain([DiscreteVariable(name="dna", values=["TTAAACTGAA", "ACTGTATAACTG", "ACTGACTG"])])
    data = np.array([[0], [1], [2], [2]]) #this data MUST be a 2d array -> otherwise id doesn't work
    d = Table.from_numpy(domain=domain, X=data)

    #set the data
    ow.set_data(d)

    #show this widget -> currently empty
    #ow.show()

    #setup and show distance matrix
    disp = OWDistanceMatrix()
    disp.set_distances(ow.raw_output)
    disp.show()
    a.exec_()


