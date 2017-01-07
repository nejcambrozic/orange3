import numpy as np
from AnyQt.QtCore import Qt


from Orange.data import Table, Domain, StringVariable
import Orange.misc
from Orange.widgets.utils.sql import check_sql_input
from Orange.widgets.widget import OWWidget, Msg
from Orange.widgets import gui
from Orange.widgets.settings import Setting


class OWMultipleSequenceAlignment(OWWidget):
    name = "Multiple Sequence Alignment"
    description = "Compute a matrix of pairwise sequence alignment distance."
    icon = "icons/MultipleSeqAlignment.svg"

    align_score = 0
    align_setting = Setting(False)
    misalign_score = 1
    misalign_setting = Setting(False)
    indel_score = 1
    indel_setting = Setting(False)

    # Spinner arguments: label, value, minval, maxval, checked
    score_settings = (('Custom alignment score', 'align_score', -100, 0, 'align_setting'),
                      ('Custom misalignment score', 'misalign_score', 0, 100, 'misalign_setting'),
                      ('Custom indel score', 'indel_score', 0, 100, 'indel_setting'))

    inputs = [("Data", Orange.data.Table, "set_data")]
    outputs = [("Distances", Orange.misc.DistMatrix),
               ("Strings", Orange.data.Table)]

    raw_output = Orange.misc.DistMatrix(np.array([]))

    class Error(OWWidget.Error):
        no_discrete_features = Msg("No discrete features")
        empty_data = Msg("Empty data (shape = {})")
        no_meta_data = Msg("One column must be meta")

    class Warning(OWWidget.Warning):
        ignoring_discrete = Msg("Ignoring discrete features")
        imputing_data = Msg("Imputing missing values")

    def __init__(self):
        super().__init__()

        self.data = None

        box = gui.vBox(self.controlArea, "Parameters")

        for lab, val, minval, maxval, chck in self.score_settings:
            gui.spin(box, self, val, minval, maxval, label=lab, checked=chck)

        gui.button(box, self, 'Update parameters', callback=self._invalidate, default=True)

    def _invalidate(self):
        if not self.misalign_setting:
            self.misalign_score = 1
        if not self.align_setting:
            self.align_score = 0
        if not self.indel_setting:
            self.indel_score = 1
        self.commit()

    @check_sql_input
    def set_data(self, data):
        """
        Set the input data set from which to compute the distances
        """
        self.data = data
        self.commit()

    def commit(self):
        self.send("Distances", self.compute_alignment(self.data))
        self.send("Strings", self.data)

    def edit_distance(self, s, p):
        """
        Returns the minimum edit distance of alignment of s and p
        """
        # init matrices
        dp = np.zeros((len(s) + 1, len(p) + 1), dtype=int)

        # init first column and row of dynamic programming table
        for col in range(1, len(s) + 1):
            dp[col][0] = col
        for row in range(1, len(p) + 1):
            dp[0][row] = row

        # compute dynamic programming table
        for i in range(1, len(s) + 1):
            for j in range(1, len(p) + 1):
                dp[i, j] = min(dp[i - 1, j] + self.indel_score,
                               dp[i, j - 1] + self.indel_score,
                               # dp[i - 1, j - 1] + (s[i - 1] != p[j - 1]))
                               dp[i - 1, j - 1] + (s[i - 1] != p[j - 1]) * self.misalign_score + (s[i - 1] == p[j - 1])
                               * self.align_score)

        # min edit distance is most bottom right element of dp
        min_distance = dp[len(s), len(p)]

        return min_distance

    def compute_alignment(self, data):
        # Check data
        self.clear_messages()
        if data is None:
            return
        if len(data.domain.metas) < 1:
            self.Error.no_meta_data()
            return
        if not any(a.is_discrete for a in data.domain.attributes):
            self.Error.no_discrete_features()
            return

        # HACK
        n = data.approx_len()
        outdata = np.zeros([n, n])
        for i, row in enumerate(data):
            for j, rowCompare in enumerate(data[i+1::], i+1):
                dist = self.edit_distance(str(row[0].value), str(rowCompare[0].value))
                outdata[i, j] = dist
                outdata[j, i] = dist

        labels = Table.from_list(
            Domain([], metas=[StringVariable("label")]),
            [[item] for item in data.get_column_view(data.domain.metas[0])[0]])

        """ Mock """
        self.raw_output = Orange.misc.DistMatrix(data=np.array(outdata), row_items=labels)
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
    from Orange.data.domain import Domain, DiscreteVariable, StringVariable
    from Orange.widgets.unsupervised.owdistancematrix import OWDistanceMatrix

    a = QApplication(sys.argv)
    ow = OWMultipleSequenceAlignment()

    # setup test data
    domain = Domain([DiscreteVariable(name="dnaSeq", values=["TTAAACTGAA", "ACTGTATAACTG", "ACTGACTG"])], [],
                        [StringVariable(name="dnaName")])
    data = np.array([[0], [1], [2], [2]])  # this data MUST be a 2d array -> otherwise id doesn't work
    metas = np.array([["dna1"], ["dna2"], ["dna3"], ["dna3"]])
    d = Table.from_numpy(domain=domain, X=data, metas=metas)

    # set the data
    ow.set_data(d)

    # show this widget -> currently empty
    # ow.show()

    # setup and show distance matrix
    disp = OWDistanceMatrix()
    disp.set_distances(ow.raw_output)
    disp.show()
    a.exec_()
