import struct
from numpy import fromfile, meshgrid, power, sqrt, exp
from numpy.fft import fft2, fftshift, rfft2, ifft2, ifftshift
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
            fft = fft2(framedata)
            #print(fft)
            print(type(fft))
            finalsize = min(len(fft), len(fft[0]))
            print('finalsize is ' + str(finalsize))
            c = finalsize / 2

            #xgrid, ygrid = meshgrid(x, y)
            ffconstant = fftshift(fft2(fft))
            print('ffconstant is ' + str(ffconstant))
            for i in range(finalsize):
                for j in range(finalsize):
                    #treat xv[j,i], yv[j,i]
                    tempf = sqrt(power(fft[i][j] - c, 2) + power(fft[i][j] - c, 2)) / finalsize
                    tempz = exp(-0.5 * bfactor * power(tempf, 2));
                    fft[i][j] = ffconstant[i][j] * tempz


            return ifft2(ifftshift(fft)).real


        

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
test.bfactor(0, 150)