#! /usr/bin/env python
from nptdms import TdmsFile
import nptdms
from matplotlib import pyplot
from Cylinder_Geometry import Cylinder_Geometry


#This opens and reads tdms files
tdms_file=TdmsFile('4_17_2018_t3_rpm9000.tdms')
channel_dev1_ai7 = tdms_file.object('This is the initial test on the Dyno, We are testing the code and see what the data produces. ambiat pressure 101.52 Kpa','Dev1/ai7')
#channel_time = tdms_file.object('This is the initial test on the Dyno, We are testing the code and see what the data produces. ambiat pressure 101.52 Kpa','Time*')
#time = channel_time.data
data = channel_dev1_ai7.data

C=Cylinder_Geometry()
volume=C.cylinder_volume_func()
o=2340
pyplot.plot(volume,data[o:7200+o])
#pyplot.plot(data[o:7200+o])
pyplot.show()

#print (data)


#for name, value in channel_object.properties.items():
  #  print("{0}: {1}".format(name, value))








