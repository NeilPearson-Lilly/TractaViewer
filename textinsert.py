# -*- coding: utf-8 -*-

# Created by: PyQt5 UI code generator 5.5
#
# WARNING! All changes made in this file will be lost!

from PyQt5 import QtCore, QtGui, QtWidgets

class Ui_Dialog(object):
    def setupUi(self, Dialog):
        Dialog.setObjectName("Dialog")
        Dialog.resize(379, 239)
        Dialog.setMinimumSize(QtCore.QSize(379, 239))
        self.gridLayout = QtWidgets.QGridLayout(Dialog)
        self.gridLayout.setObjectName("gridLayout")
        self.textEdit = QtWidgets.QTextEdit(Dialog)
        self.textEdit.setObjectName("textEdit")
        self.gridLayout.addWidget(self.textEdit, 0, 0, 2, 1)
        self.groupBox = QtWidgets.QGroupBox(Dialog)
        self.groupBox.setMinimumSize(QtCore.QSize(110, 91))
        self.groupBox.setMaximumSize(QtCore.QSize(110, 91))
        self.groupBox.setObjectName("groupBox")
        self.verticalLayout_2 = QtWidgets.QVBoxLayout(self.groupBox)
        self.verticalLayout_2.setObjectName("verticalLayout_2")
        self.radioButton_symbol = QtWidgets.QRadioButton(self.groupBox)
        self.radioButton_symbol.setChecked(True)
        self.radioButton_symbol.setObjectName("radioButton_symbol")
        self.buttonGroup = QtWidgets.QButtonGroup(Dialog)
        self.buttonGroup.setObjectName("buttonGroup")
        self.buttonGroup.addButton(self.radioButton_symbol)
        self.verticalLayout_2.addWidget(self.radioButton_symbol)
        self.radioButton_ensemblID = QtWidgets.QRadioButton(self.groupBox)
        self.radioButton_ensemblID.setObjectName("radioButton_ensemblID")
        self.buttonGroup.addButton(self.radioButton_ensemblID)
        self.verticalLayout_2.addWidget(self.radioButton_ensemblID)
        self.radioButton_uniprotID = QtWidgets.QRadioButton(self.groupBox)
        self.radioButton_uniprotID.setObjectName("radioButton_uniprotID")
        self.buttonGroup.addButton(self.radioButton_uniprotID)
        self.verticalLayout_2.addWidget(self.radioButton_uniprotID)
        self.gridLayout.addWidget(self.groupBox, 0, 1, 1, 1)
        self.buttonBox = QtWidgets.QDialogButtonBox(Dialog)
        self.buttonBox.setOrientation(QtCore.Qt.Vertical)
        self.buttonBox.setStandardButtons(QtWidgets.QDialogButtonBox.Cancel|QtWidgets.QDialogButtonBox.Ok)
        self.buttonBox.setObjectName("buttonBox")
        self.gridLayout.addWidget(self.buttonBox, 1, 1, 1, 1)
        self.lineEdit = QtWidgets.QLineEdit(Dialog)
        self.lineEdit.setReadOnly(True)
        self.lineEdit.setObjectName("lineEdit")
        self.gridLayout.addWidget(self.lineEdit, 2, 0, 1, 1)

        self.retranslateUi(Dialog)
        self.buttonBox.accepted.connect(Dialog.accept)
        self.buttonBox.rejected.connect(Dialog.reject)
        QtCore.QMetaObject.connectSlotsByName(Dialog)

    def retranslateUi(self, Dialog):
        _translate = QtCore.QCoreApplication.translate
        Dialog.setWindowTitle(_translate("Dialog", "Gene input"))
        self.groupBox.setTitle(_translate("Dialog", "Input type"))
        self.radioButton_symbol.setText(_translate("Dialog", "HGNC Symbol"))
        self.radioButton_ensemblID.setText(_translate("Dialog", "Ensembl ID"))
        self.radioButton_uniprotID.setText(_translate("Dialog", "Uniprot ID"))

