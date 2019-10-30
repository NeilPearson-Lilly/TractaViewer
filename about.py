# -*- coding: utf-8 -*-

# Created by: PyQt5 UI code generator 5.5
#
# WARNING! All changes made in this file will be lost!

from PyQt5 import QtCore, QtGui, QtWidgets

class Ui_About(object):
    def setupUi(self, About):
        About.setObjectName("About")
        About.resize(400, 300)
        self.buttonBox = QtWidgets.QDialogButtonBox(About)
        self.buttonBox.setGeometry(QtCore.QRect(30, 240, 341, 32))
        self.buttonBox.setOrientation(QtCore.Qt.Horizontal)
        self.buttonBox.setStandardButtons(QtWidgets.QDialogButtonBox.Ok)
        self.buttonBox.setObjectName("buttonBox")
        self.label = QtWidgets.QLabel(About)
        self.label.setGeometry(QtCore.QRect(120, 10, 161, 21))
        self.label.setTextFormat(QtCore.Qt.AutoText)
        self.label.setObjectName("label")
        self.label_logo = QtWidgets.QLabel(About)
        self.label_logo.setGeometry(QtCore.QRect(16, 14, 91, 91))
        self.label_logo.setText("")
        self.label_logo.setObjectName("label_logo")
        self.label_2 = QtWidgets.QLabel(About)
        self.label_2.setGeometry(QtCore.QRect(120, 40, 191, 151))
        self.label_2.setWordWrap(True)
        self.label_2.setObjectName("label_2")

        self.retranslateUi(About)
        self.buttonBox.accepted.connect(About.accept)
        self.buttonBox.rejected.connect(About.reject)
        QtCore.QMetaObject.connectSlotsByName(About)

    def retranslateUi(self, About):
        _translate = QtCore.QCoreApplication.translate
        About.setWindowTitle(_translate("About", "About TractaViewer"))
        self.label.setText(_translate("About", "<html><head/><body><p><span style=\" font-size:12pt; font-weight:600;\">TractaViewer</span></p></body></html>"))
        self.label_2.setText(_translate("About", "<html><head/><body><p>A tool to mine for, organise and visualise data relating to the suitability of genes and their protein products as drug targets.</p><p>Created by:</p><p>-Neil Pearson</p><p>-Karim Malki</p><p>-Nathan Lawless</p><p>-David Collier</p></body></html>"))

