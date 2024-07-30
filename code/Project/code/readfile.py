from multiprocessing.dummy import current_process
from pydoc import pathdirs
import sys
import os
import re
from PyQt5 import QtCore
from PyQt5.Qt import *
import shutil
from ffmpy3 import FFmpeg
global filePath


def mkdir(file_name):
    global filePath
    file_path = filePath+file_name
    folder_exist=os.path.exists(filePath+file_name)

    if not folder_exist:                   #判断是否存在文件夹如果不存在则创建为文件夹
        os.makedirs(filePath+file_name)            #makedirs 创建文件时如果路径不存在会创建这个路径
        return filePath+file_name
    else:
        print("OK")
        return filePath+file_name

def process(FrameDngList, DirPos, FinalDir):
    global filePath
    print(FrameDngList)
    for i in FrameDngList:
        shutil.move(filePath+"/"+i, DirPos)

    # 写text
    f = open(DirPos + "/reference_frame.txt","w")   #设置文件对象
    f.write("0")
    f.close() #关闭文件

    # 开始操作
    cmd = "python runHdrplus.py -i " + DirPos + " -o " + FinalDir + " -m full -v 2"
    print(cmd)
    # os.system("conda list")
    os.system(cmd)

    for i in FrameDngList:
        shutil.move(DirPos+"/"+i, filePath)
    print(DirPos)

    os.remove(DirPos + "//reference_frame.txt")  # path是文件的路径，如果这个路径是一个文件夹，则会抛出OSError的错误，这时需用用rmdir()来删除
    os.rmdir(DirPos) # 需要为空

class MainWidget(QWidget):
    def __init__(self, parent=None):
        super(MainWidget, self).__init__(parent)
 
        #获取系统所有文件
        self.model01 = QFileSystemModel()
        #进行筛选只显示文件夹，不显示文件和特色文件
        self.model01.setFilter(QtCore.QDir.Dirs|QtCore.QDir.NoDotAndDotDot)
        self.model01.setRootPath('')
 
        #定义创建左边窗口
        self.treeView1 = QTreeView(self)
        self.treeView1.setModel(self.model01)
        for col in range(1, 4):
            self.treeView1.setColumnHidden(col, True)
        self.treeView1.doubleClicked.connect(self.initUI)
 
        #定义创建右边窗口
        self.model02 = QStandardItemModel()
        self.treeView2 = QTreeView(self)
        self.treeView2.setModel(self.model02)
 
        #将创建的窗口进行添加
        self.main_layout = QVBoxLayout()
        self.layout = QHBoxLayout()

        self.bt1 = QPushButton("Start Making Video",self)
        self.bt1.clicked.connect(self.start)
        # self.layout.addWidget(self.bt1)

        self.bt2 = QLabel("No .Dng in Processing",self)
        # self.layout.addWidget(self.bt2)

        self.layout.addWidget(self.treeView1)
        self.layout.addWidget(self.treeView2)
        self.main_layout.addLayout(self.layout)
        self.main_layout.addWidget(self.bt2)
        self.main_layout.addWidget(self.bt1)

        self.setLayout(self.main_layout)
 
 
    def initUI(self, Qmodelidx):
        global filePath
        #每次点击清空右边窗口数据
        self.model02.clear()
        #定义一个数组存储路径下的所有文件
        PathData = []
        #获取双击后的指定路径
        filePath = self.model01.filePath(Qmodelidx)
        # List窗口文件赋值
        PathDataName = self.model02.invisibleRootItem()
        #拿到文件夹下的所有文件
        PathDataSet = os.listdir(filePath)
        #进行将拿到的数据进行排序
        PathDataSet.sort()
        #遍历判断拿到的文件是文件夹还是文件，Flase为文件，True为文件夹
        for Data in range(len(PathDataSet)):
            if os.path.isdir(filePath + '\\' + PathDataSet[Data]) == False:
                PathData.append(PathDataSet[Data])
            elif os.path.isdir(filePath + '\\' + PathDataSet[Data]) == True:
                print('2')
        #将拿到的所有文件放到数组中进行右边窗口赋值。
        for got in range(len(PathData)):
            gosData = QStandardItem(PathData[got])
            PathDataName.setChild(got, gosData)
        print(filePath)

    def update_info(self, now_index, total_index):
        self.bt2.setText("Processing: " + str(now_index) + "/" + str(total_index))

    def spilt2dir(self, filePath, DngList):
        for i in range(len(DngList)-10):
            print(i)
        tempFinal = mkdir("/tempFinal")
        for i in range(len(DngList)-10):
            current_list = []
            tempDirPos = mkdir("/tempDir%06d" %i)
            for j in range(10):
                current_list.append(DngList[i+j])
            process(current_list, tempDirPos, tempFinal)
            self.update_info(i, len(DngList)-10)
            QApplication.processEvents()
        self.bt2.setText("Process End")
        self.bt1.setEnabled(True)
        self.make_video(tempFinal)
        print("End Process")

    def make_video(self, DirPos):
        print(DirPos)
        ffin = FFmpeg(inputs={DirPos+"/final_tempDir%06d.jpg": '-f image2 -r 30'}, outputs={DirPos + '/Final.mp4': None})
        print(ffin.cmd)
        os.system(ffin.cmd)


    def start(self):
        global filePath
        PathDataSet = os.listdir(filePath)
        print(PathDataSet)
        DngList = []
        for file in PathDataSet:
            if(".dng" in file):
                DngList.append(file)
        print(DngList)
        DngList.sort()
        print(DngList)
        if (len(DngList)==0):
            print("There is no \".dng\" file!")
        else:
            self.bt1.setEnabled(False)
            self.spilt2dir(filePath, DngList)

 
if __name__ == "__main__":
    app = QApplication(sys.argv)
    window = MainWidget()
    window.resize(600, 400)
    window.show()
    sys.exit(app.exec_())
