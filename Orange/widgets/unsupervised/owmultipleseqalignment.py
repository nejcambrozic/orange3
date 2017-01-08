import numpy as np
from AnyQt.QtCore import Qt

import Orange.misc
from Orange.data import Table, Domain, StringVariable, DiscreteVariable
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
               ("Alignments", Orange.data.Table)]

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
        self.alignment = None

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
        sc, al = self.compute_alignment(self.data)
        self.send("Distances", sc)
        self.send("Alignments", al)

    def edit_distance(self, s, t):
        def sig(c1, c2):
            """Sigma function - returns cost"""
            if c1 == "" or c2 == "":
                ret = self.indel_score
            elif c1 == c2:
                ret = self.align_score
            else:
                ret = self.misalign_score
            return ret

        # create matrix
        # M = [[0 for _ in range(len(t) + 1)] for _ in range(len(s) + 1)]
        M = np.zeros((len(s) + 1, len(t) + 1), dtype=int)
        directions = {}

        # init sides
        for i in range(0, len(s) + 1):
            M[i, 0] = sum([sig(s[k], "") for k in range(i)])
            directions[i,0] = 2
        for j in range(0, len(t) + 1):
            M[0, j] = sum([sig("", t[k]) for k in range(j)])
            directions[0, j] = 1

        directions[0,0] = 4

        # fill matrix
        for i in range(1, len(s) + 1):
            for j in range(1, len(t) + 1):
                val,dire = min((M[i - 1, j] + sig(s[i - 1], ""), 2),
                               (M[i, j - 1] + sig("", t[j - 1]), 1),
                               (M[i - 1, j - 1] + sig(s[i - 1], t[j - 1]), 0))
                M[i, j] = val
                directions[i, j] = dire

        i = len(s)
        j = len(t)

        bufs = ""
        buft = ""

        while i > 0 and j > 0:
            d = directions[i, j]
            if d == 4: break
            elif d == 0:
                i -= 1; j -= 1
                bufs += s[i]
                buft += t[j]
            elif d == 1:
                j-= 1
                bufs += "_"
                buft += t[j]
            elif d == 2:
                i -= 1
                bufs += s[i]
                buft += "_"

        bufs = bufs[::-1]
        buft = buft[::-1]
        alignstr = bufs + "\n" + buft
        min_distance = M[len(s), len(t)]

        return min_distance, alignstr

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

        n = data.approx_len()
        outdata = np.zeros([n, n])
        alignment = np.zeros([n*n, 1])
        counter = 0

        values = []
        for i, row in enumerate(data):
            for j, rowCompare in enumerate(data[i + 1::], i + 1):
                dist, al = self.edit_distance(str(row[0].value), str(rowCompare[0].value))
                outdata[i, j] = dist
                outdata[j, i] = dist
                alignment[i*n+j, 0] = counter
                alignment[j*n+i, 0] = counter
                values.append(al)  # Append the alignement to domain values
                counter += 1


        labels = Table.from_list(
            Domain([], metas=[StringVariable("label")]),
            [[item] for item in data.get_column_view(data.domain.metas[0])[0]])

        domain = Domain([DiscreteVariable(name="Alignments", values=values)])

        self.alignment = Table.from_numpy(domain=domain, X=alignment)
        self.raw_output = Orange.misc.DistMatrix(data=np.array(outdata), row_items=labels)
        return self.raw_output, self.alignment


"""
 Reads some data from a table file and sets it as
 the input of this widget.
 Then pipes the output of this widget into the
 distance matrix widget to show the results
"""
if __name__ == "__main__":
    import sys
    from AnyQt.QtWidgets import QApplication
    from Orange.data import Table
    from Orange.data.domain import Domain, DiscreteVariable, StringVariable
    from Orange.widgets.unsupervised.owalignmentdistancematrix import OWAlignmentDistanceMatrix

    a = QApplication(sys.argv)
    ow = OWMultipleSequenceAlignment()

    # setup test data
    domain = Domain([DiscreteVariable(name="dnaSeq", values=["ABCD", "ABC", "AAAD"])], [],
                    [StringVariable(name="dnaName")])
    data = np.array([[0], [1], [2], [2]])  # this data MUST be a 2d array -> otherwise id doesn't work
    metas = np.array([["dna1"], ["dna2"], ["dna3"], ["dna4"]])
    d = Table.from_numpy(domain=domain, X=data, metas=metas)

    ow.align_score = -1
    ow.misalign_score = 1
    ow.indel_score = 1
    # set the data
    ow.set_data(d)

    # show this widget -> NOTE: changing values here will not update AlignDistMat below
    # ow.show()

    # setup and show distance matrix
    disp = OWAlignmentDistanceMatrix()
    disp.set_distances(ow.raw_output)
    disp.set_alignments(ow.alignment)
    disp.show()
    a.exec_()

