#-*- coding:utf-8 -*-
import sys
from os.path import expanduser
from PyQt4 import QtGui, QtCore

class mainWin(QtGui.QMainWindow):
    """
    main window
    """
    def __init__(self, width, height):
        super(mainWin, self).__init__()
        self.init(width, height)

    def init(self, width, height):
        self.resize(width, height)
        self.setWindowTitle("Video Painterly Rendering System")

        # move window to the center of screen
        qr = self.frameGeometry()
        cp = QtGui.QDesktopWidget().availableGeometry().center()
        qr.moveCenter(cp)
        self.move(qr.topLeft())

        # init menu bar actions
        openAction = QtGui.QAction("Open", self)
        openAction.triggered.connect(self.open)

        saveAction = QtGui.QAction("Save", self)
        saveAction.triggered.connect(self.save)
        
        exitAction = QtGui.QAction("Exit", self)
        exitAction.triggered.connect(QtGui.qApp.quit)

        aboutAction = QtGui.QAction("About", self)
        aboutAction.triggered.connect(QtGui.qApp.quit)
        
        # init menu bar
        menubar = self.menuBar()
        fileMenu = menubar.addMenu("&File")
        fileMenu.addAction(openAction)
        fileMenu.addAction(saveAction)
        fileMenu.addAction(exitAction)

        helpMenu = menubar.addMenu("&Help")
        helpMenu.addAction(aboutAction)

        self.IRPanel = IRPanel(self, width, height)
        # self.VRPanel = VRPanel(self, width, height)

        self.showIRPanel()
        self.show()

    def open(self):
        fname = QtGui.QFileDialog.getOpenFileName(self, "Open Video", expanduser("~"),
                                                  "Videos(*.mp4 *.avi)")
        print "open video: " + fname

    def save(self):
        fname = QtGui.QFileDialog.getSaveFileName(self, "Save Video", expanduser("~"),
                                                  "Videos(*.mp4 *.avi)")
        print "save to: " + fname

    def showIRPanel(self):
        self.setCentralWidget(self.IRPanel)

    def showVRPanel(self):
        self.setCentralWidget(self.VRPanel)

class IRPanel(QtGui.QWidget):
    """
    Image Rendering panel
    """
    def __init__(self, parent, width, height):
        super(IRPanel, self).__init__(parent)
        self.init(width, height)

    def init(self, width, height):
        IRlayout = QtGui.QHBoxLayout()

        # init control panel
        ctrlLayout = QtGui.QVBoxLayout()
        ctrlLayout.setContentsMargins(0, 0, 0, 0)
        ctrlLayout.setSpacing(0)
        
        # min stroke slider and label
        self.minStrokeLbl = QtGui.QLabel(self)
        self.setMinStrokeLbl(2)
        minStrokeSLD = QtGui.QSlider(QtCore.Qt.Horizontal, self)
        minStrokeSLD.setSingleStep(1)
        minStrokeSLD.setRange(2, 32)
        minStrokeSLD.valueChanged[int].connect(self.setMinStrokeLbl)
        
        ctrlLayout.addWidget(self.minStrokeLbl)
        ctrlLayout.addWidget(minStrokeSLD)

        # max stroke slider and label
        self.maxStrokeLbl = QtGui.QLabel(self)
        self.setMaxStrokeLbl(2)
        maxStrokeSLD = QtGui.QSlider(QtCore.Qt.Horizontal, self)
        maxStrokeSLD.setSingleStep(1)
        maxStrokeSLD.setRange(2, 32)
        maxStrokeSLD.valueChanged[int].connect(self.setMaxStrokeLbl)
        
        ctrlLayout.addWidget(self.maxStrokeLbl)
        ctrlLayout.addWidget(maxStrokeSLD)

        # stroke step slider and label
        self.strokeLbl = QtGui.QLabel(self)
        self.setStrokeLbl(2)
        strokeSLD = QtGui.QSlider(QtCore.Qt.Horizontal, self)
        strokeSLD.setSingleStep(1)
        strokeSLD.setRange(2, 32)
        strokeSLD.valueChanged[int].connect(self.setStrokeLbl)
        
        ctrlLayout.addWidget(self.strokeLbl)
        ctrlLayout.addWidget(strokeSLD)

        # min stroke len slider and label
        self.minStrokeLenLbl = QtGui.QLabel(self)
        self.setMinStrokeLenLbl(2)
        minStrokeLenSLD = QtGui.QSlider(QtCore.Qt.Horizontal, self)
        minStrokeLenSLD.setSingleStep(1)
        minStrokeLenSLD.setRange(2, 32)
        minStrokeLenSLD.valueChanged[int].connect(self.setMinStrokeLenLbl)
        
        ctrlLayout.addWidget(self.minStrokeLenLbl)
        ctrlLayout.addWidget(minStrokeLenSLD)

        # max stroke len slider and label
        self.maxStrokeLenLbl = QtGui.QLabel(self)
        self.setMaxStrokeLenLbl(2)
        maxStrokeLenSLD = QtGui.QSlider(QtCore.Qt.Horizontal, self)
        maxStrokeLenSLD.setSingleStep(1)
        maxStrokeLenSLD.setRange(2, 32)
        maxStrokeLenSLD.valueChanged[int].connect(self.setMaxStrokeLenLbl)
        
        ctrlLayout.addWidget(self.maxStrokeLenLbl)
        ctrlLayout.addWidget(maxStrokeLenSLD)

        IRlayout.addLayout(ctrlLayout)
        
        self.setLayout(IRlayout)

    def setMinStrokeLbl(self, value):
        self.minStrokeLbl.setText("minimum stroke size: %d" %value)

    def setMaxStrokeLbl(self, value):
        self.maxStrokeLbl.setText("maximum stroke size: %d" %value)

    def setStrokeLbl(self, value):
        self.strokeLbl.setText("stroke size step: %d" %value)

    def setMinStrokeLenLbl(self, value):
        self.minStrokeLenLbl.setText("minimum stroke length: %d" %value)

    def setMaxStrokeLenLbl(self, value):
        self.maxStrokeLenLbl.setText("maximum stroke length: %d" %value)

class VRPanel(QtGui.QWidget):
    """
    Video Rendering panel
    """
    def _init__(self, parent, width, height):
        super(VRPanel, self).__init__(parent)
        self.init(width, height)

    def init(self, width, height):
        pass
    
if __name__ == "__main__":
    app = QtGui.QApplication(sys.argv)
    width = 800
    height = 600
    win = mainWin(width, height)
    sys.exit(app.exec_())
