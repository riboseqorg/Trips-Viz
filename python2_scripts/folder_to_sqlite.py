import os
import subprocess

basepath = "/home/DATA/www/tripsviz/tripsviz/trips_shelves/riboseq/homo_sapiens/Calviello16"


total_files = 0
for item in os.walk(basepath):
    folder = item[0]
    files = item[2]
    for filename in files:
        if "shelf" in filename:
            total_files += 1
              
            print "calling conversion with {}/{}".format(folder,filename)
            print "total_files", total_files
            filepath =  "{}/{}".format(folder,filename)
            subprocess.call("python shelf_to_sqlite.py {}".format(filepath),shell=True)
