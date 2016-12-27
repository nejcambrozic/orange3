import numpy as np
from AnyQt.QtCore import Qt

import Orange.data
import Orange.misc
from Orange.widgets.utils.sql import check_sql_input
from Orange.widgets.widget import OWWidget, Msg


class OWDistances(OWWidget):
    name = "Multiple Sequence Alignment"
    description = "Compute a matrix of pairwise sequence alignment distance."
    icon = "icons/MultipleSeqAlignment.svg"

    inputs = [("Data", Orange.data.Table, "set_data")]
    outputs = [("Distances", Orange.misc.DistMatrix)]

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

    def compute_alignment(self, data):
        """ Mock """
        d = Orange.misc.DistMatrix(np.array([[1, 2, 3], [4, 5, 6]], np.int32))
        # TODO align given sequences and return edit distance matrix
        return d
