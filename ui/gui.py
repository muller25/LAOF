#-*- coding:utf-8 -*-
import sys
from os.path import expanduser
from PyQt4 import QtGui

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
        saveAction.triggered.connect(QtGui.qApp.quit)
        
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
        
        self.show()

    def open(self):
        fname = QtGui.QFileDialog.getOpenFileName(self, "Open Video", expanduser("~"))

        print "open video: " + fname

    def save(self):
        fname = QtGui.QFileDialog.getSaveFileName(self, "Save Video", expanduser("~"))

        print "save to: " + fname

class IRPanel(QtGui.QWidget):
    """
    Image Rendering panel
    """
    def __init__(self, parent, width, height):
        super(IRPanel, self).__init__(parent)
        self.init(width, height)

    def init(self, width, height):
        pass

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
    width = 500
    height = 500
    win = mainWin(width, height)
    sys.exit(app.exec_())
