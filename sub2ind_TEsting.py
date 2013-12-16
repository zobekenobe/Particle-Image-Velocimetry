# -*- coding: utf-8 -*-
"""
Created on Tue Dec 10 01:56:48 2013

@author: zobekenobe
"""

import numpy as np

def sub2ind(shapeOfArray,dim):
        n,m = shapeOfArray
        return (n*(dim[0])+dim[1]+1)

a = np.random.random_integers(0,25,(5,5)) # iFinal = 1600 *1600
p = np.ones((3,3)) # partInterp= 99 * 99 
# Size = same as p and range should be [0,24]
x = np.random.random_integers(0,4,(3,3))
y = np.random.random_integers(0,4,(3,3))



linInd = sub2ind(a.shape,(x[:,1],y[:,1]))
#A = a.T.reshape(25,1)
#
#print p
#p[:,0] = A[linInd].T  
#print p

import PIL as pl
import matplotlib.pyplot as plt
import matplotlib.cm as cm


imgg = plt.imread('/home/zobekenobe/Downloads/imgs/650mm_imgs.33iuhov0.000001.tif')
#plt.imshow(imgg, cmap = cm.Greys_r)
plt.imshow(imgg)
plt.show()

#img = pl.Image.open('/home/zobekenobe/Downloads/imgs/650mm_imgs.33iuhov0.000001.tif','r')
#im2 = pl.Image.open('/home/zobekenobe/Downloads/img:0.png','r')
#im2 = im2.convert('L')
#img = img.resize((800,600))
#imgFinal= pl.Image.blend(img,im2,0.15)
plt.savefig('zobe')