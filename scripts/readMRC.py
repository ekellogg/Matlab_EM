import struct
from numpy import fromfile, meshgrid, power, sqrt, exp, vectorize, divide, multiply, subtract, add, log, mean
from numpy.fft import fft2, fftshift, rfft2, ifft2, ifftshift, irfft2, fftn, ifftn
import sys
import numpy as np
import scipy
from joblib import Parallel, delayed
import multiprocessing
from PIL import Image
import pylab
import time

class ReadMRC:
    def __init__(self,filename):
         self.numbytes1=56           # 56 long ints
         self.numbytes2=80*10          # 10 header lines of 80 chars each
         self.filename=filename
         self.frames = list()
         self.test = False

    def bfactorReal(self, framenum, bfactor):
        return self.bfactor(self.frames[framenum], bfactor)
    #For Car Test
    def bfactorTest(self, image, bfactor):
        return self.bfactor(image, bfactor)

    #function[corr_coeff] = calculate_corr_vector(elems_to_compare,frstack)
    def calculate_corr_vector(elements, framestack):

        #n = length(elems_to_compare);
        n = len(elements)

        #corr_coeff = zeros((n*(n-1))/2,2);
        corr_coef = np.zeros(n * (n - 1) / 2, 2)

        #vectorize for speed
        #ndx=1;
        ndx = 1

        #for( i = 1:(n-2) )
        for i in range(0, n - 3):

            #array1 = (frstack(:,:,i) - mean(mean(frstack(:,:,i))));
            array1 = framestack[:][:][i] - np.mean(np.mean(framestack[:][:][i]))

            #[xx,yy] = arrayfun(@(ind) findmax(xcorr2_fft(sub_mean(frstack(:,:,i)),sub_mean(frstack(:,:,ind)) )),...
            [xx, yy]
                                              (i+2):n);
          %  tic; [xx,yy] = arrayfun(@(ind) refine_autocorr_peak(xcorr2_fft(sub_mean(frstack(:,:,i)),sub_mean(frstack(:,:,ind)))),...
          %                             (i+1):n); toc;
            corr_coeff(ndx:(ndx+length(xx)-1),:) = [-(xx-size(frstack,1))' -(yy-size(frstack,2))'];
            ndx=ndx+length(xx);
        end
        corr_coeff = corr_coeff(1:(ndx-1),:);
    end



    def create_coeff_matrix(self, numframes):
        pass
        #function[A] = create_coeff_matrix(i)

        #inputLength = len(input)

        #A = gpuArray(zeros(length(i),length(i)-1));
        A = np.zeros((numframes*2, numframes -1))

        print('start is ' + str(A))
        #ii = i(1);
        ii = 0

        #ndx = 1;
        ndx = 0

        #while(ii <= i(length(i)))
        for ii in range(numframes/2):

            #jj = ii+1;
            #jj = ii 

            #while(jj < i(length(i))  )
            for jj in range(ii+1, numframes):

                #A(ndx,ii:jj) = 1;
                print('ndx is ' + str(ndx) + ' ii jj is  ' + str(ii) + ' ' + str(jj))

                A[ndx][ii:jj] = 1
                print('a is now ' + str(A) + '\n')

                #%do not align frames that are adjacent
                #jj = jj + 1;
                #jj += 1

                #ndx = ndx + 1;
                ndx += 1
            #ii = ii + 1;



        return A


    def driftCorrection(self, image):
        pass
        #nfr = size(fr,3);
        nfr = len(image)

        #create matrix A
        #A = create_coeff_matrix(1:nfr);
        A = self.create_coeff_matrix() ######HMMMMMM

        return A

        #b = calculate_corr_vector_v2(1:nfr,fr);
        #b = 
        #img_shift = inv(transpose(A)*A)*transpose(A)*b;

    def bfactor(self, framedata, bfactor):
        framedata = np.array(framedata)
        finalsize = min(len(framedata), len(framedata[0]))

        def squarifyFun(a):
            cut1 = a[:finalsize]
            for i in cut1:
                i = i[:finalsize]
            return cut1

        framedata = squarifyFun(framedata)

        #print('finalsize is ' + str(finalsize))
        c = finalsize / 2

        xv, yv = meshgrid([x for x in range(finalsize)], [x for x in range(finalsize)],indexing='xy')

        f = divide(sqrt(pow(subtract(xv, c), 2) + pow(subtract(yv, c), 2)), finalsize)
        z = exp(multiply(pow(f,2), -0.5*bfactor));
        startline = time.time()
        #######
        bimg = ifft2(ifftshift(  multiply(fftshift(fft2(framedata)), z)          )).real
        #######

        print('Slow Line line took ' + str(time.time() - startline))
        #print('bimg is ' + str(bimg))
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


#Car Example
#test = ReadMRC('notimportant')
#im = Image.open("f.png").convert('L')
#a = np.asarray(im)
# Before #
#pylab.imshow(a, cmap = pylab.get_cmap('gray'))
#pylab.show()
#print(a)
#b = test.bfactorTest(a, 550)
# After #
#pylab.imshow(b, cmap = pylab.get_cmap('gray'))
#pylab.show()


## Example on frame one of test-stack.mrc
#st = time.time()
#test = ReadMRC('test-stack.mrc')
#test.read()
#a = test.bfactorReal(0, 150)
#ft = time.time() - st
#print(ft)
#pylab.imshow(a, cmap = pylab.get_cmap('gray'))
#pylab.show()