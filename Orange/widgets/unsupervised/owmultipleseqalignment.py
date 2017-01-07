import numpy as np
from AnyQt.QtCore import Qt

from Orange.data import Table, Domain, StringVariable
import Orange.misc
from Orange.widgets.utils.sql import check_sql_input
from Orange.widgets.widget import OWWidget, Msg
from Orange.widgets import gui, settings
from PyQt4.QtGui import QGridLayout


class OWMultipleSequenceAlignment(OWWidget):
    name = "Multiple Sequence Alignment"
    description = "Compute a matrix of pairwise sequence alignment distance."
    icon = "icons/MultipleSeqAlignment.svg"

    inputs = [("Data", Orange.data.Table, "set_data")]
    outputs = [("Distances", Orange.misc.DistMatrix),
               ("Strings", Orange.data.Table)]

    raw_output = Orange.misc.DistMatrix(np.array([]))

    miss_penalty = 1
    skip_penalty = 1

    #DELETE THIS
    update_skip = settings.Setting(0)
    update_miss = settings.Setting(0)

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
        # No gui

        smbox = gui.vBox(None, margin=0)
        ssbox = gui.vBox(None, margin=0)
        self.spin_miss = gui.spin(
            smbox, self, "miss_penalty", minv=1, maxv=100,
            controlWidth=80, alignment=Qt.AlignRight, callback=self.update_miss)

        self.spin_skip = gui.spin(
            ssbox, self, "skip_penalty", minv=1, maxv=100,
            controlWidth=80, alignment=Qt.AlignRight, callback=self.update_skip)

        buttonbox = gui.vBox(None, margin=0)
        self.apply_button = gui.button(
            buttonbox, self, "Run", callback=self.commit)


        self.layout().addWidget(gui.widgetLabel(smbox, "Miss penalty: "))
        self.layout().addWidget(self.spin_miss)
        self.layout().addWidget(gui.widgetLabel(ssbox, "Skip penalty: "))
        self.layout().addWidget(self.spin_skip)
        self.layout().addWidget(self.apply_button)

        self.layout().setSizeConstraint(self.layout().SetFixedSize)


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

        if self is not None:
            skip = self.skip_penalty
            miss = self.miss_penalty
        else:
            skip = 1
            miss = 1

        # compute dynamic programming table
        for i in range(1, len(s) + 1):
            for j in range(1, len(p) + 1):
                dp[i, j] = min(dp[i - 1, j] + skip,
                               dp[i, j - 1] + skip,
                               dp[i - 1, j - 1] + (s[i - 1] != p[j - 1])*miss)

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
    domain = Domain([DiscreteVariable(name="dnaSeq", values=["TTAACTGAA", "ACTGTATAACTG", "ACTGACTG"])], [],
                        [StringVariable(name="dnaName")])
    data = np.array([[0], [1], [2], [2]])  # this data MUST be a 2d array -> otherwise id doesn't work
    metas = np.array([["dna1"], ["dna2"], ["dna3"], ["dna3"]])
    d = Table.from_numpy(domain=domain, X=data, metas=metas)

    # set the data
    ow.set_data(d)

    # show this widget -> currently empty
    ow.show()


    # setup and show distance matrix
    disp = OWDistanceMatrix()
    disp.set_distances(ow.raw_output)
    disp.show()
    a.exec_()
