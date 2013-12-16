# -*- coding: utf-8 -*-
"""
Created on Sun Nov 24 15:08:58 2013

@author: Ldalab
"""
import os,csv,scipy,cv2
from scipy import signal
import matplotlib.pyplot as plt
import matplotlib.mlab as ml
from numpy import *
import numpy as np
import pylab as pyl
imgDir=os.listdir('imgs')
img_no=len(imgDir)
import matplotlib.colors as mcolors
import matplotlib.cm as cm


def markPLIFparticle(IFull):
    filteredImage = cv2.medianBlur(IFull,5)
    binaryImageComplement = cv2.threshold(filteredImage,22,255,cv2.THRESH_BINARY_INV)[1]
#    filledImage = ndimage.morphology.binary_fill_holes(binaryImageComplement[1]).astype(int)
#    creating the structuring element
    se = cv2.getStructuringElement(cv2.MORPH_ELLIPSE,(5,5));iterate = 3
    imgEr = cv2.erode(binaryImageComplement,se,iterations = iterate)
    imgDi = cv2.dilate(imgEr,se,iterations = iterate)
    return imgDi
    
def prepareVec(x,y,u,v,chk):
#    Find the spiky vectores based on global threshold: Global Validation
    mag = scipy.sqrt(u**2+v**2)
    spike = mag>400
    u[spike]=0;v[spike]=0;
    
#    Reshape the column based on the criteria of repition of 'x'
    
    ind = ml.find(x==min(x))
    wid = ind[2]-ind[1]
    x2d = mat(scipy.reshape(x,(wid,-1)))
    y2d = mat(scipy.reshape(y,(wid,-1)))
    u2d = mat(scipy.reshape(u,(wid,-1)))
    v2d = mat(scipy.reshape(v,(wid,-1)))
    c2d = mat(scipy.reshape(chk,(wid,-1)))

    return x2d,y2d,u2d,v2d,c2d

def sub2ind(shapeOfArray,dim):
        n,m = shapeOfArray
        return (n*(dim[0])+dim[1])

for i in [98]:
#    Reading data from CSV file and converting to vectors
    print 'img:'+str(i)
    os.chdir('650mm_vecs')
    if i <10:
        csvfile=open('650mm_vecs.33iyhgwy.00000'+str(i)+'.csv','rb')
    else: 
        csvfile=open('650mm_vecs.33iyhgwy.0000'+str(i)+'.csv','rb')
    wb=csv.reader(csvfile)
    a=[];row_no=0
    os.chdir('..')
    os.chdir('imgs')
#    img=scipy.ndimage.imread(imgDir[i])
    if i <10:
        img = cv2.imread('650mm_imgs.33iuhov0.00000'+str(i)+'.tif',cv2.CV_LOAD_IMAGE_GRAYSCALE)
    else:
        img = cv2.imread('650mm_imgs.33iuhov0.0000'+str(i)+'.tif',cv2.CV_LOAD_IMAGE_GRAYSCALE)
    iFinal = markPLIFparticle(img)
#    Reading in the vectors
    for row in wb:
        row_no+=1
        if row_no>6:
            a.append(map(float,row))
    a=scipy.array(a)
    x,y,u,v,chk=a[:,0],a[:,1],a[:,2],a[:,3],a[:,4]
    os.chdir('..')    
    
    x_min=min(x);x_max=max(x)
    y_min=min(y);y_max=max(y)
    u_mean=scipy.mean(u)
    v_mean=scipy.mean(v)
    x2d,y2d,u2d,v2d,c2d=prepareVec(x,y,u,v,chk)
    u2d=signal.wiener(u2d);v2d=signal.wiener(v2d)
    x_width=x_max-x_min;y_width=y_max-y_min
    
    XI = (x2d-x_min)/x_width*1599
    YI = (y2d-y_min)/y_width*1599
    
    
    rows = floor(YI)
    cols = floor(XI)
    
    partInterp = np.zeros(shape(XI))
    
    IFinal=iFinal.T.reshape(2560000,1)
#    IFinal=iFinal.reshape(2560000,1)
    for j in range(shape(XI)[1]):
        linInd = sub2ind(shape(iFinal),(scipy.array(rows[:,j]),scipy.array(cols[:,j])))
        partInterp[:,j] = IFinal[map(int,floor(linInd))].T
    
        
    partMask = np.fliplr(partInterp.T)
    
    partMask = np.flipud(partMask)
    partMask = np.fliplr(partMask)    
    
    u2d[partMask.astype(bool)] = 0; v2d[partMask.astype(bool)] = 0
    
    XI,YI = np.meshgrid(np.linspace(x_min,x_max,1600),np.linspace(y_min,y_max,1600))
    skip = 3
#    plt.pcolor(XI[1:-1:skip][1:-1:skip],YI[1:-1:skip][1:-1:skip],np.flipud(img[1:-1:skip][1:-1:skip]))  
    plt.hold(True);
    Q=pyl.quiver(x2d,y2d,u2d,v2d,color='g')   
    plt.show()
    
    plt.axis([0,50,0,50])
    plt.title('processed image'+str(i))  
    os.chdir('ScipyFinalResult')
    
    plt.hold();
    plt.savefig('img:'+str(i),)
    os.chdir('..')
    
    