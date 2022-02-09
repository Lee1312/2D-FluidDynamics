import numpy as np
from functions import *

class Fluid():

    def __init__(self,n,dt,diff,visc):      #INITIALIZE THE FLUID
        self.dt=dt      #TIME STEP
        self.diff=diff  #DIFFUSION 
        self.visc=visc  #VISCOSITY
        self.n=int(n)   #SIZE

        #MATRIXES NEEDED FOR REPRESENTATION AND CALCULATION
        self.densityZERO = np.zeros(shape=(self.n,self.n))      
        self.density = np.zeros(shape=(self.n,self.n))
        
        self.Vx = np.zeros(shape=(self.n,self.n))
        self.Vy = np.zeros(shape=(self.n,self.n))
        
        self.Vx0 = np.zeros(shape=(self.n,self.n))
        self.Vy0 = np.zeros(shape=(self.n,self.n))


    def Step(self,iter):
        difuse(1, self.Vx0, self.Vx, self.visc, self.dt, iter, self.n)      #DIFUSES VELOCITY VECTORS
        difuse(2, self.Vy0, self.Vy, self.visc, self.dt, iter, self.n)      

        project(self.Vx0, self.Vy0, self.Vx, self.Vy, iter, self.n)         #PROJECTS THE VELOCITY

        advect(1, self.Vx, self.Vx0, self.Vx0, self.Vy0, self.dt, self.n)   #ADVECTS THE VELOCITY
        advect(2, self.Vy, self.Vy0, self.Vx0, self.Vy0, self.dt, self.n)
        
        project(self.Vx, self.Vy, self.Vx0, self.Vy0, iter, self.n)         #PROJECTS THE VELOCITY AGAIN
        
        difuse(0, self.densityZERO, self.density, self.diff, self.dt, iter, self.n)     #DIFUSES THE DYE
        advect(0, self.density, self.densityZERO, self.Vx, self.Vy, self.dt, self.n)    #ADVECTS THE DYE

    def addDensVel(self,mfcx,mfcy):
        #Matrix of mouse position is like:
                        

            #                    neg
            #     1.st Quadrant   | 3.nd Quadrant
            #         x=negative  |  x=positive
            #         y=negative  |  y=negative
            #                     |
            #  neg --------------------------------- pos 
            #     2.rd Quadrant   | 4.th Quadrant
            #         x=negative  |  x=positive
            #         y=positive  |  y=positive
            #                     |
            #                    pos

        if(mfcx<=0):
            if(mfcy<=0):
                #PRVI KVADRANT
                teta=np.arctan2(mfcx,mfcy)
                x=np.sin(teta)
                y=np.cos(teta)
            else:
                #DRUGI KVADRANT
                teta=np.arctan2(mfcx,mfcy)
                x=np.sin(teta)
                y=np.cos(teta)
        if(mfcx>0):
            if(mfcy<=0):
                #TRECI KVADRANT
                teta=np.arctan2(mfcx,mfcy)
                x=np.sin(teta)
                y=np.cos(teta)
            else:
                #CETVRTI KVADRANT
                teta=np.arctan2(mfcx,mfcy)
                x=np.sin(teta)
                y=np.cos(teta)

        self.addVelocity(32,32,x*10,y*10)       #ADDS VELOCITY VECTOR TO MIDDLE 
        self.addVelocity(32,33,x*10,y*10)       #BASED ON MOUSE POSITION
        self.addVelocity(33,32,x*10,y*10)
        self.addVelocity(33,33,x*10,y*10)
        self.addDye(32,32,500.0)                #ADDS DYE TO MIDDLE
        self.addDye(32,33,500.0)
        self.addDye(33,32,500.0)
        self.addDye(33,33,500.0)

    def addDye(self,x,y,kolicina):
        self.density[x,y]+=kolicina

    def addVelocity(self,x,y,kolX,kolY):
        self.Vx[x,y]+=kolX
        self.Vy[x,y]+=kolY