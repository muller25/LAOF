#-*- coding:utf-8 -*-
import sys
from os.path import expanduser
from PyQt4 import QtGui, QtCore
from PyQt4.phonon import Phonon

class mainWin(QtGui.QMainWindow):
    """
    main window
    """
    def __init__(self, width, height):
        super(mainWin, self).__init__()
        self.init(width, height)

    def init(self, width, height):
        self.setFixedSize(width, height)
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

        self.IRPanel = IRPanel(width, height)
        self.VRPanel = VRPanel(width, height)

        tab = QtGui.QTabWidget(self)
        tab.addTab(self.IRPanel, "Image Render")
        tab.addTab(self.VRPanel, "Video Render")

        self.setCentralWidget(tab)
        self.show()

    def open(self):
        fname = QtGui.QFileDialog.getOpenFileName(self, "Open Video", expanduser("~"))#,
#                                                  "Videos(*.mp4 *.avi)")

        if not fname.isEmpty():
            print "open video: " + fname
            self.VRPanel.showVideo(fname)
            # self.IRPanel.showImage(fname)
        else:
            print "do not choose any video"

    def save(self):
        fname = QtGui.QFileDialog.getSaveFileName(self, "Save Video", expanduser("~"),
                                                  "Videos(*.mp4 *.avi)")
        if not fname.isEmpty():
            print "save to: " + fname
        else:
            print "do not choose any video"

class IRPanel(QtGui.QWidget):
    """
    Image Rendering panel
    """
    def __init__(self, width, height):
        super(IRPanel, self).__init__()
        self.init(width, height)

    def init(self, width, height):
        IRlayout = QtGui.QHBoxLayout()

        # init control panel
        ctrlLayout = QtGui.QVBoxLayout()
        
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

        # add scale factor to push buttons to the bottom
        ctrlLayout.addStretch(1)

        # show button
        runBtn = QtGui.QPushButton("Run", self)
        runBtn.clicked.connect(self.runBtnClick)
        ctrlLayout.addWidget(runBtn)

        IRlayout.addLayout(ctrlLayout)
        
        # add image view
        scrollArea = QtGui.QScrollArea(self)
        scrollArea.setAlignment(QtCore.Qt.Alignment(QtCore.Qt.AlignCenter))
        scrollArea.setFixedWidth(width - 220)
        self.imgView = QtGui.QLabel(self)
        scrollArea.setWidget(self.imgView)
        IRlayout.addWidget(scrollArea)

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

    def runBtnClick(self):
        print "Image rendering..."
        
    def showImage(self, fname):
        img = QtGui.QPixmap(fname)
        self.imgView.resize(img.size())
        self.imgView.setPixmap(img)
        
class VRPanel(QtGui.QWidget):
    """
    Video Rendering panel
    """
    def __init__(self, width, height):
        super(VRPanel, self).__init__()
        self.init(width, height)

    def init(self, width, height):
        VRlayout = QtGui.QHBoxLayout()

        # init control panel
        ctrlLayout = QtGui.QVBoxLayout()

        lblWidth = 130
        # spatial smoothness
        tmpHBox = QtGui.QHBoxLayout()
        spatialLbl = QtGui.QLabel("sptial smooth", self)
        spatialLbl.setFixedWidth(lblWidth)
        self.spatial = QtGui.QLineEdit(self)
        tmpHBox.addWidget(spatialLbl)
        tmpHBox.addWidget(self.spatial)
        ctrlLayout.addLayout(tmpHBox)
        
        # temporal smoothness
        tmpHBox = QtGui.QHBoxLayout()
        temporalLbl = QtGui.QLabel("temporal smooth", self)
        temporalLbl.setFixedWidth(lblWidth)
        self.temporal = QtGui.QLineEdit(self)
        tmpHBox.addWidget(temporalLbl)
        tmpHBox.addWidget(self.temporal)
        ctrlLayout.addLayout(tmpHBox)

        # segment similarity
        tmpHBox = QtGui.QHBoxLayout()
        segmentLbl = QtGui.QLabel("segment smooth", self)
        segmentLbl.setFixedWidth(lblWidth)
        self.segment = QtGui.QLineEdit(self)
        tmpHBox.addWidget(segmentLbl)
        tmpHBox.addWidget(self.segment)
        ctrlLayout.addLayout(tmpHBox)

        # run button
        ctrlLayout.addStretch(1)
        runBtn = QtGui.QPushButton("Run", self)
        runBtn.clicked.connect(self.runBtnClick)
        ctrlLayout.addWidget(runBtn)

        VRlayout.addLayout(ctrlLayout)
        
        # add video player
        vlayout = QtGui.QVBoxLayout()
        self.player = Phonon.VideoPlayer(self)
        self.seekSlider = Phonon.SeekSlider(self)

        vlayout.addWidget(self.player, 1)
        vlayout.addWidget(self.seekSlider)
        VRlayout.addLayout(vlayout, 1)

        self.setLayout(VRlayout)

    def runBtnClick(self):
        print "video rendering..."

    def showVideo(self, fname):
        self.player.load(Phonon.MediaSource(fname))
        self.seekSlider.setMediaObject(self.player.mediaObject())
        # self.player.play()
    
if __name__ == "__main__":
    app = QtGui.QApplication(sys.argv)
    app.setApplicationName("Video Painterly Rendering System")
    width = 1000
    height = 800
    win = mainWin(width, height)
    sys.exit(app.exec_())
