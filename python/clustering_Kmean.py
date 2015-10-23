from numpy import array
from scipy.cluster.vq import vq, kmeans, whiten
import os
import sys
import os.path
import numpy as np
from scipy.cluster.vq import vq,kmeans,whiten,kmeans2
from copy import deepcopy
# this code is used to cluster my Tibet result using the K-mean clustering method

#====================
class dataFile:
  #-constructor, define instance variable here
  #-the 1st two numbers in columnlist should be the column number for lon&lat
  def __init__(self,filename,columnlist):
    self.filename=filename;
    self.columnlist=columnlist;

    colLon=columnlist[0];
    colLat=columnlist[1];
    if(not os.path.exists(filename)):
      print "file %s does not exist\n"%(filename);
      sys.exit();
    #-find the lonmin/max and latmin/max, and the step, so that I can initiate the array  
    lonmin=1e10;latmin=1e10;lonmax=-1e10;latmax=-1e10;step=1e10;
    count=0;
    for line in open(filename):
      count+=1;    
      l=line.rstrip().split();
      lon=float(l[colLon]);
      lat=float(l[colLat]);
      lonmin=lon if lon<lonmin else lonmin;
      latmin=lat if lat<latmin else latmin;
      lonmax=lon if lon>lonmax else lonmax;
      latmax=lat if lat>latmax else latmax;
      if(count>=2 ): # this is not the best way to find step, if the data is not stored in order, then this would have problem
        dif=np.fabs(lon-tlon)
	if(dif>1e-3):
	  step=dif if dif<step else step;
      tlon=lon;	        
    Nlon=int((lonmax-lonmin+step*0.1)/step)+1;  
    Nlat=int((latmax-latmin+step*0.1)/step)+1;  

    self.lonmin=lonmin;
    self.lonmax=lonmax;
    self.latmin=latmin;
    self.latmax=latmax;
    self.step=step;
    self.Nlon=Nlon;
    self.Nlat=Nlat;
    self.Ncol=len(columnlist);
    
    #-initialize the array [lon][lat][column]
    self.dataArray=[[[-999. for k in range(self.Ncol)] for j in range(Nlat)] for i in range(Nlon)]
    self.dataArray=np.array(self.dataArray)

  #-read in the data into array, each row is one sample, # of column represent # of features
  def readData(self):
    for line in open(self.filename):
      l=line.rstrip().split();
      lon=float(l[self.columnlist[0]]);
      lat=float(l[self.columnlist[1]]);
      ilon=int((lon-self.lonmin)/self.step+1e-2);
      ilat=int((lat-self.latmin)/self.step+1e-2);
      
      for i in range(self.Ncol):
        idata=self.columnlist[i];
        if(l[idata]=="nan"):
          continue
	value=float(l[idata]);      
        if(np.isnan(value)):
          continue
        self.dataArray[ilon][ilat][i]=value;

  #--get a data at a give location
  def getData(self,lon,lat):
    ilon=int((lon-self.lonmin)/self.step+1e-2);
    ilat=int((lat-self.latmin)/self.step+1e-2);
    return self.dataArray[ilon][ilat];
   
  #--print a data (1d array)
  def printData(self,lon,lat):
    print "printData ",lon,lat
    array1d=self.getData(lon,lat);
    str=" "
    for i in range(self.Ncol):
      str=str+"%.3f "%(array1d[i])
    print str 
  
  #-- remove points that has incomplete misfit record
  #-- produce a M*N matrix, M samples, N features
  def cleanDataForKmean(self):
    self.cleanMatrix=[]
    for i in range(self.Nlon):
      for j in range(self.Nlat):
          count=0;
	  for k in range(self.Ncol):
	      if(np.fabs(self.dataArray[i][j][k]+999)<1e-2):
	        count+=1
          if(count==0):
	      self.cleanMatrix.append(self.dataArray[i][j])     
    self.cleanMatrix=np.array(self.cleanMatrix)     

  #-- remove points that has incomplete misfit record
  #-- remove points inside the given file
  #-- produce a M*N matrix, M samples, N features
  def cleanDataForKmean(self,fileBasinLoc):
    self.cleanMatrix=[]
    if(not os.path.exists(fileBasinLoc)):
      print "basin point location file %s does not exists! "%(fileBasinLoc);
      sys.exit();
    idBasin=[];  
    for line in open(fileBasinLoc):
      l=line.rstrip().split();
      lon=float(l[0]);
      lat=float(l[1]);
      ilon=int((lon-self.lonmin)/self.step+1e-2);
      ilat=int((lat-self.latmin)/self.step+1e-2);
      idBasin.append(" %d_%d "%(ilon,ilat));
    for i in range(self.Nlon):
      for j in range(self.Nlat):
        str=" %d_%d "%(i,j);
        if(str in idBasin):
          continue
        count=0;
        for k in range(self.Ncol):
          if(np.fabs(self.dataArray[i][j][k]+999)<1e-2):
	    count+=1
        if(count==0):
          self.cleanMatrix.append(self.dataArray[i][j])
    self.cleanMatrix=np.array(self.cleanMatrix)


  #--write out labeled data
  def writeDataLabel(self,outfilenm):
      fout=open(outfilenm,"w");
      Ninfo=len(self.cleanMatrix[0]);
      if(Ninfo==self.Ncol):
        print"label not atted yet! cannot print label"
	sys.exit();
      for i in range(len(self.cleanMatrix)):
        info=self.cleanMatrix[i];
      	fout.write("%8.2f %8.2f "%(info[0],info[1]))
	for j in range(Ninfo-self.Ncol):
	  fout.write(" %d"%(info[self.Ncol+j]))
	fout.write("\n")
      fout.close();
      
  #--
  def doKmean(self,matrix,Ncentroid,fileinfo):
    matrixWhitened=whiten(matrix);
    centroid1,distortion=kmeans(matrixWhitened,Ncentroid,iter=1000)
    label,_=vq(matrixWhitened,centroid1)
    #centroids,label=kmeans2(matrixWhitened,centroid1,minit='matrix');
    #centroids,label=kmeans2(matrixWhitened,Ncentroid,minit='points');
    label=np.reshape(label,(len(label),1))
    fileinfo.cleanMatrix=np.hstack((file1.cleanMatrix,label));
    return label


  #--toString	
  def __repr__(self):
    str="=============\nfilename=%s\nlon=(%.2f,%.1f), lat=(%.2f,%.2f), step=%.2f\nNlon=%d, Nlat=%d, Ncol=%d\n=============\n"%(self.filename,self.lonmin,self.lonmax,self.latmin,self.latmax,self.step,self.Nlon,self.Nlat,self.Ncol)
    return str

#====================
#def combineArray():
  #maybe this is not a method for dataFile

#--prepare data. read in all the info, location, misfit 1-lay, misfit 2-lay, misfit 3-lay, misfit difference, amplitude of ani
# # of row is # of samples, # of column is # of features
# go through each file, read in all the information, and then combine those array together

#-----Main-------------

#--read in the misfit difference information
subname="scale3.freeta.vpvs.3";
date="aug11";
type="best";
filename1="/projects/jixi7887/work/YunNanTibet/Analyze/get_map_info/info_misfit_1lay2lay_diff_%s.%s.%s.txt"%(subname,date,type); #difference for misfit all, R, Ramp, Rphi, L
columnlist1=[0,1,3,4,5,6]; # column number for lon, lat, MD(misfit_diff)_all, MD_Riso, MD_Ramp, MD_Rphi, MD_Liso
Ncentroid=3; # number of centroids
############
file1=dataFile(filename1,columnlist1);
file1.readData();
#file1.cleanDataForKmean();
file1.cleanDataForKmean("pointInBasin.txt");
print file1
#
matrix=deepcopy(file1.cleanMatrix);
misfitMatrix=matrix[:,2:] # exclude the 1st two columns
locationMatrix=matrix[:,0:2] # the 1st two columns, 0,1
misfitlocMatrix=matrix[:,0:(file1.Ncol-1)]

#--do K-mean clustering 
labelMisfit=file1.doKmean(misfitMatrix,Ncentroid,file1); # the label is already appended to file1.cleanMatrix

labelMisLoc=file1.doKmean(misfitlocMatrix,Ncentroid,file1);

hybridMatrix=np.hstack((locationMatrix,labelMisfit));
labelHybrid=file1.doKmean(hybridMatrix,Ncentroid,file1);

file1.writeDataLabel("DataKMeanLabel.txt");

print file1.cleanMatrix
print hybridMatrix


# cluster data based on their features (misfit diff for LoveIso, RayIso, RayAmp, RayPhi) and location (lon,lat)

# cluster data feature-clustered data based on their location

#--apply K-mean

#--output result and make a plot (?)

