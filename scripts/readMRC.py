import struct
from numpy import fromfile, meshgrid, power, sqrt, exp, vectorize, divide, multiply
from numpy.fft import fft2, fftshift, rfft2, ifft2, ifftshift
import sys
import numpy as np
import scipy
from joblib import Parallel, delayed
import multiprocessing


#from math import exp, sqrt

class ReadMRC:
    def __init__(self,filename):
         self.numbytes1=56           # 56 long ints
         self.numbytes2=80*10          # 10 header lines of 80 chars each
         self.filename=filename
         self.frames = list()

    def bfactor(self, framenum, bfactor):
        #function[bimg] = apply_bfact(img,bfactor)
        #%add padding in FFT to minimize interpolation error
        #c = ceil(size(img))./2;
        #x = 1:1:size(img,1);
        #y = 1:1:size(img,2);
        #[X,Y] = meshgrid(x,y);
        #f = (sqrt((X-c(1)).^2 + (Y-c(2)).^2))./(size(img,1));
        #z = exp(-0.5*bfactor*(f.^2));
        #bimg = real(ifft2(ifftshift(  fftshift(fft2(img)).*z  )));
        #page number 259

        framedata = self.frames[framenum]
        finalsize = min(len(framedata), len(framedata[0]))

        def squarifyFun(a):
            cut1 = a[:finalsize]
            for i in cut1:
                i = i[:finalsize]
            return cut1

        fft = squarifyFun(framedata)

        #print('finalsize is ' + str(finalsize))
        c = finalsize / 2

        xv, yv = meshgrid([x for x in range(finalsize)], [x for x in range(finalsize)], sparse=False, indexing='ij')

        f = divide(sqrt(pow(xv-c, 2) + pow(yv - c, 2)), finalsize)
        #print('f is ' + str(f))
        z = exp(-0.5*bfactor*pow(f,2));
        bimg = ifft2(ifftshift(  multiply(fftshift(fft2(fft)), z)  )).real
        print('bimg is ' + str(bimg))
        return bimg

        

    def read(self):
        input_image=open(self.filename,'rb')
        self.header1=input_image.read(self.numbytes1*4)
        self.header2=input_image.read(self.numbytes2)
        byte_pattern='=' + 'l'*self.numbytes1   #'=' required to get machine independent standard size
        self.dim=struct.unpack(byte_pattern,self.header1)[:3]   #(dimx,dimy,dimz)
        self.imagetype=struct.unpack(byte_pattern,self.header1)[3]  #0: 8-bit signed, 1:16-bit signed, 2: 32-bit float, 6: unsigned 16-bit (non-std)
        imtype='unknown'

        if (self.imagetype == 0):
            imtype='b'
        elif (self.imagetype ==1):
            imtype='h'
        elif (self.imagetype ==2):
            imtype='f4'
        elif (self.imagetype ==6):
            imtype='H'

        print('imtype is ' + imtype)

        input_image_dimension=(self.dim[1],self.dim[0])  #2D images assumed
        print('input_image_dimension is ' + str(input_image_dimension))

        #print(len(fromfile(file=input_image,dtype=imtype,count=self.dim[0]*self.dim[1])))
        xysize = self.dim[0]*self.dim[1]
        
        #self.image_data=fromfile(file=input_image,dtype=imtype,count=xysize*20)[2*xysize:3*xysize].reshape(input_image_dimension)

        self.frames.append(fromfile(file=input_image,dtype=imtype,count=xysize*20)[0:xysize].reshape(input_image_dimension))
        input_image.close()

##########
# Start  #
##########
test = ReadMRC('test-stack.mrc')
test.read()
#print('self.dim is ' + str(test.dim))
#print(test.frames[0])
print(test.bfactor(0, 150))