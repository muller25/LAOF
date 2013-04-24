#-*- coding:utf-8 -*-
import sys, os, subprocess, shutil
from PyQt4 import QtGui, QtCore
from PyQt4.phonon import Phonon

cacheDir = os.path.join(os.path.expanduser("~"), "cache")
imDir = os.path.join(cacheDir, "im/")
flowDir = os.path.join(cacheDir, "flow/")
renderDir = os.path.join(cacheDir, "render/")
outVideo = os.path.join(cacheDir, "render.avi")
imPattern = os.path.join(imDir, "%03d.jpg")
flowPattern = os.path.join(flowDir, "%s%03d.yml")
renderPattern = os.path.join(renderDir, "out%03d.jpg")

if not os.path.exists(imDir):
    os.makedirs(imDir)
if not os.path.exists(flowDir):
    os.makedirs(flowDir)
if not os.path.exists(renderDir):
    os.makedirs(renderDir)

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
        tab.addTab(self.VRPanel, "Video Render")
        tab.addTab(self.IRPanel, "Image Render")
        tab.currentChanged.connect(self.tabChanged)
        self.setCentralWidget(tab)
        self.show()

    def open(self):
        fname = QtGui.QFileDialog.getOpenFileName(self, "Open Video", os.path.expanduser("~"),
                                                  "Videos(*.mp4 *.avi)")
        if not fname.isEmpty():
            print "open video: " + fname
            self.VRPanel.showVideo(fname)
        else:
            print "do not choose any video"

    def save(self):
        fname = QtGui.QFileDialog.getSaveFileName(self, "Save Video", os.path.expanduser("~"),
                                                  "Videos(*.avi)")
        if not fname.isEmpty():
            print "save to: " + fname
            if os.path.exists(outVideo):
                shutil.copyfile(outVideo, fname)
            else:
                print "you should generate video first"
        else:
            print "do not choose any video"

    def tabChanged(self, curIdx):
        # ir panel
        if curIdx == 1:
            self.IRPanel.showImage(self.VRPanel.screenshot())

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

        # algorithm switcher
        tmpLbl = QtGui.QLabel("choose algorithm to use", self)
        ctrlLayout.addWidget(tmpLbl)
        
        algoSwitcher = QtGui.QListWidget(self)
        algoSwitcher.addItem(self.tr("CDIR"))
        
        ctrlLayout.addWidget(algoSwitcher)
        
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
        self.setStrokeLbl(1)
        strokeSLD = QtGui.QSlider(QtCore.Qt.Horizontal, self)
        strokeSLD.setSingleStep(1)
        strokeSLD.setRange(1, 30)
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
        pr = os.path.join("bin", "pr")
        inFile = os.path.join(cacheDir, "in.jpg")
        outFile = os.path.join(cacheDir, "out.jpg")

        img = self.imgView.pixmap()
        if img is None: return

        img.save(inFile, "JPG")
        proc = subprocess.Popen([pr, inFile, outFile])
        proc.wait()

        res = QtGui.QPixmap()
        res.load(outFile)
        self.showImage(res)
        
    def showImage(self, img):
        self.imgView.resize(img.size())
        self.imgView.setPixmap(img)
        
class VRPanel(QtGui.QWidget):
    """
    Video Rendering panel
    """
    def __init__(self, width, height):
        super(VRPanel, self).__init__()
        self.init(width, height)
        self.videoName = None
        
    def init(self, width, height):
        VRlayout = QtGui.QHBoxLayout()

        # init control panel
        ctrlLayout = QtGui.QVBoxLayout()
        lblWidth = 130

        # algorithm switcher
        tmpLbl = QtGui.QLabel("choose algorithm to use", self)
        ctrlLayout.addWidget(tmpLbl)
        
        algoSwitcher = QtGui.QListWidget(self)
        algoSwitcher.addItem(self.tr("LAOF"))
        
        ctrlLayout.addWidget(algoSwitcher)

        # spatial smoothness
        tmpHBox = QtGui.QHBoxLayout()
        spatialLbl = QtGui.QLabel("sptial smooth", self)
        spatialLbl.setFixedWidth(lblWidth)
        self.spatial = QtGui.QLineEdit(self)
        self.spatial.setText(str(0.05))
        tmpHBox.addWidget(spatialLbl)
        tmpHBox.addWidget(self.spatial)
        ctrlLayout.addLayout(tmpHBox)
        
        # temporal smoothness
        tmpHBox = QtGui.QHBoxLayout()
        temporalLbl = QtGui.QLabel("temporal smooth", self)
        temporalLbl.setFixedWidth(lblWidth)
        self.temporal = QtGui.QLineEdit(self)
        self.temporal.setText(str(0.015))
        tmpHBox.addWidget(temporalLbl)
        tmpHBox.addWidget(self.temporal)
        ctrlLayout.addLayout(tmpHBox)

        # segment similarity
        tmpHBox = QtGui.QHBoxLayout()
        segmentLbl = QtGui.QLabel("segment similarity", self)
        segmentLbl.setFixedWidth(lblWidth)
        self.segment = QtGui.QLineEdit(self)
        self.segment.setText(str(0.8))
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
        self.player.setFixedWidth(width-220)
        self.seekSlider = Phonon.SeekSlider(self)
        self.media = None
        
        playout = QtGui.QHBoxLayout()
        playout.addStretch(1)
        self.playIcon = QtGui.QIcon("icons/play.png")
        self.pauseIcon = QtGui.QIcon("icons/pause.png")
        self.playBtn = QtGui.QPushButton(self.playIcon, self.tr(""), self)
        self.playBtn.clicked.connect(self.playBtnClick)
        playout.addWidget(self.playBtn)
        playout.addStretch(1)
        
        vlayout.addWidget(self.player, 1)
        vlayout.addWidget(self.seekSlider)
        vlayout.addLayout(playout)
        
        VRlayout.addLayout(vlayout, 1)
        self.setLayout(VRlayout)

    def runBtnClick(self):
        vr = os.path.join("bin", "vr")
        flow = os.path.join("bin", "biof")
        video2im = os.path.join("bin", "video2images")
        im2video = os.path.join("bin", "image2video")

        if self.videoName is None:
            print "please choose a video"
            return

        # video to image
        v2iProc = subprocess.Popen([video2im, self.videoName, imDir])
        v2iProc.wait()

        # running laof
        nfile = len([name for name in os.listdir(imDir) if os.path.isfile(os.path.join(imDir, name))])
        ofProc = subprocess.Popen([flow, imPattern, flowDir, str(0), str(nfile)])
        ofProc.wait()

        # video rendering
        vrProc = subprocess.Popen([vr, imPattern, flowPattern, renderDir, str(0), str(nfile)])
        vrProc.wait()

        i2vProc = subprocess.Popen([im2video, renderPattern, str(0), str(nfile), outVideo, str(15)])
        i2vProc.wait()

        self.showVideo(outVideo)
        
        print "done"
        
    def playBtnClick(self):
        if self.player.isPlaying():
            self.player.pause()
            self.playBtn.setIcon(self.playIcon)
            return
        
        if self.media is None:
            self.player.stop()
            self.playBtn.setIcon(self.playIcon)
            return

        self.player.play()
        self.playBtn.setIcon(self.pauseIcon)
            
    def showVideo(self, fname):
        self.videoName = fname
        self.player.load(Phonon.MediaSource(fname))
        self.media = self.player.mediaObject()
        
        self.playBtn.setIcon(self.playIcon)
        self.seekSlider.setMediaObject(self.media)
            
    def screenshot(self):
        return QtGui.QPixmap.grabWindow(self.player.winId())
    
if __name__ == "__main__":
    app = QtGui.QApplication(sys.argv)
    app.setApplicationName("Video Painterly Rendering System")
    width = 1000
    height = 800
    win = mainWin(width, height)
    sys.exit(app.exec_())
