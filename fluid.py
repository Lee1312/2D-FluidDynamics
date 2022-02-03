import numpy as np

def difuse(b,x,x0,diff,dt,iter,n):
    a = dt * diff * (n-2) * (n-2)   
    lin_solve(b,x,x0,a,1+4*a,iter,n)
    return



def lin_solve(b,x,x0,a,c,iter,n):
    cRecip=1.0/c
    for k in range(iter):
        for i in range(1,n-1):
            for j in range(1,n-1):
                x[i,j] = (x0[i,j] + a*(x[i-1,j]+x[i+1,j]+x[i,j-1]+x[i,j+1]))*cRecip
        set_bnd(b,x,n)
    return

def set_bnd(b,x,n):
    for i in range(1,n-1):
        if(b==1):
            x[0,i]=-x[1,i]
        else:
            x[0,i]=x[1,i]
        if(b==1):
            x[n-1,i]=-x[n-2,i]
        else:
            x[n-1,i]=x[n-2,i]
             
    #for j in range(1,n-1):  
        if(b==2):
            x[i,0]=-x[i,1]
        else:
            x[i,0]=x[i,1]
        if(b==2):
            x[i,n-1]=-x[i,n-2]
        else:
            x[i,n-1]=x[i,n-2]

    x[0,0] = 0.5*(x[1,0]+x[0,1])
    x[0,n-1] = 0.5*(x[1,n-1]+x[0,n-2])
    x[n-1,0] = 0.5*(x[n-2,0]+x[n-1,1])
    x[n-1,n-1] = 0.5*(x[n-2,n-1]+x[n-1,n-2])
    return

def project(x,y,x0,y0,iter,n):

    for i in range(1,n-1):
        for j in range(1,n-1):
            y0[i,j] = -0.5*(x[i+1,j]-x[i-1,j]+y[i,j+1]-y[i,j-1])/n
            x0[i,j] = 0
        
    set_bnd(0,y0,n)
    set_bnd(0,x0,n)
    lin_solve(0,x0,y0,1,4,iter,n)

    for i in range(1,n-1):
        for j in range(1,n-1):
            x[i,j] -= 0.5*(x0[i+1,j]-x0[i-1,j])*n
            y[i,j] -= 0.5*(x0[i,j+1]-x0[i,j-1])*n
    
    set_bnd (1,x,n)
    set_bnd (2,y,n)
    return



def advect (b,d,d0,Vx,Vy,dt,n):
    dt0 = dt*(n-2)
    
    for i in range(1,n-2):
        for j in range(1,n-2):
            x = float(i-dt0*Vx[i,j])   #Oduzme vektor brzine od 
            y = float(j-dt0*Vy[i,j])   #trenutnog polozaja

            if (x<0.5): 
                x=0.5
            if (x>n+0.5): 
                x=n+0.5
            i0=np.floor(x) 
            i1=i0+1.0

            if (y<0.5): 
                y=0.5 
            if (y>n+0.5): 
                y=n+0.5 
            j0=np.floor(y)
            j1=j0+1.0
            
            s1 = x-i0
            s0 = 1.0-s1
            t1 = y-j0
            t0 = 1.0-t1
            
            i0i=int(i0)
            i1i=int(i1)
            j0i=int(j0)
            j1i=int(j1) 
            
            # if(i0i>=n-1):
            #     i0i=n-1
            # if(i1i>=n-1):
            #     i0i=n-1
            # if(j0i>=n-1):
            #     i0i=n-1
            # if(j1i>=n-1):
            #     i0i=n-1
            # if(i>=n-1):
            #     i=n-1
            # if(j>=n-1):
            #     j=n-1

            d[i,j] = s0*(t0*d0[i0i,j0i]+t1*d0[i0i,j1i])+ s1*(t0*d0[i1i,j0i]+t1*d0[i1i,j1i])

    set_bnd (b,d,n)
    return

class Fluid():
    def __init__(self,n,dt,diff,visc):
        self.dt=dt  #TIME STEP
        self.diff=diff  #DIFFUSION 
        self.visc=visc  #VISCOSITY
        self.n=int(n)

        self.densityZERO = np.zeros(shape=(self.n,self.n))
        self.density = np.zeros(shape=(self.n,self.n))
        
        self.Vx = np.zeros(shape=(self.n,self.n))
        self.Vy = np.zeros(shape=(self.n,self.n))
        
        self.Vx0 = np.zeros(shape=(self.n,self.n))
        self.Vy0 = np.zeros(shape=(self.n,self.n))

    def Step(self,iter):
        N       = self.n
        visc    = self.visc
        diff    = self.diff
        dt      = self.dt
        Vx      = self.Vx
        Vy      = self.Vy
        Vx0     = self.Vx0
        Vy0     = self.Vy0
        s       = self.densityZERO
        density = self.density

        difuse(1, Vx0, Vx, visc, dt, iter, N)
        difuse(2, Vy0, Vy, visc, dt, iter, N)
        
        project(Vx0, Vy0, Vx, Vy, iter, N)
        
        advect(1, Vx, Vx0, Vx0, Vy0, dt, N)
        advect(2, Vy, Vy0, Vx0, Vy0, dt, N)
        
        project(Vx, Vy, Vx0, Vy0, iter, N)
        
        difuse(0, s, density, diff, dt, iter, N)
        advect(0, density, s, Vx, Vy, dt, N)

    def addDensVel(self):
        self.addVelocity(32,32,20,0)
        self.addVelocity(32,33,20,0)
        self.addVelocity(33,32,20,0)
        self.addVelocity(33,33,20,0)
        self.addDye(32,32,500.0)
        self.addDye(32,33,500.0)
        self.addDye(33,32,500.0)
        self.addDye(33,33,500.0)

    def addDye(self,x,y,kolicina):
        self.density[x,y]+=kolicina

    def difuseFluid(self,iter):
        difuse(0,self.density,self.densityZERO,self.diff,self.dt,iter,self.n)
        
    def advectFluid(self):
        advect(0,self.density,self.densityZERO,self.Vx,self.Vy,self.dt,self.n)

    def DyeStep(self):
        self.addDye(32,32,500.0)
        self.addDye(32,33,500.0)
        self.addDye(33,32,500.0)
        self.addDye(33,33,500.0)

        #SWAP
        temp=np.copy(self.density)
        self.density=np.copy(self.densityZERO)
        self.densityZERO=np.copy(temp)

        self.difuseFluid(iter)

        #SWAP
        temp=np.copy(self.density)
        self.density=np.copy(self.densityZERO)
        self.densityZERO=np.copy(temp)

        self.advectFluid()

    def addVelocity(self,x,y,kolX,kolY):
        self.Vx[x,y]+=kolX
        self.Vy[x,y]+=kolY

    def difuseVelocityX(self,iter):
        difuse(1,self.Vx,self.Vx0,self.visc,self.dt,iter,self.n)

    def difuseVelocityY(self,iter):
        difuse(2,self.Vy,self.Vy0,self.visc,self.dt,iter,self.n)

    def project(self,iter):
        project(self.Vx,self.Vy,self.Vx0,self.Vy0,iter,self.n)

    def advectVelocityX(self):
        advect(1,self.Vx,self.Vx0,self.Vx0,self.Vy0,self.dt,self.n)

    def advectVelocityY(self):
        advect(2,self.Vy,self.Vy0,self.Vx0,self.Vy0,self.dt,self.n)

    def VelStep(self,iter):
        self.addVelocity(32,32,2,2)
        self.addVelocity(32,33,2,2)
        self.addVelocity(33,32,2,2)
        self.addVelocity(33,33,2,2)
 
        print(self.Vx[32,32])
        print(self.Vx[32,33])
        print(self.Vy[32,32])
        print(self.Vy[32,33])

        A=self.Vx[[29,30,31,32,33,34,35,36]][:,[29,30,31,32,33,34,35,36]]
        B=self.Vx0[[29,30,31,32,33,34,35,36]][:,[29,30,31,32,33,34,35,36]]
        C=self.Vy[[29,30,31,32,33,34,35,36]][:,[29,30,31,32,33,34,35,36]]
        D=self.Vy0[[29,30,31,32,33,34,35,36]][:,[29,30,31,32,33,34,35,36]]
        print("BRZINE PRE SMENE")
        
        print(A)
        print(B)
        print(C)
        print(D)


        #SWAP Vx
        temp=np.copy(self.Vx0)
        self.Vx0=np.copy(self.Vx)
        self.Vx=np.copy(temp)

        A=self.Vx[[29,30,31,32,33,34,35,36]][:,[29,30,31,32,33,34,35,36]]
        B=self.Vx0[[29,30,31,32,33,34,35,36]][:,[29,30,31,32,33,34,35,36]]
        
        print("BRZINE POSLE SMENE Vx")
        print(A)
        print(B)

        self.difuseVelocityX(iter)

        A=self.Vx[[29,30,31,32,33,34,35,36]][:,[29,30,31,32,33,34,35,36]]
        B=self.Vx0[[29,30,31,32,33,34,35,36]][:,[29,30,31,32,33,34,35,36]]

        print("BRZINA POSLE DIFUSA Vx")
        print(A)
        print(B)

        #SWAP Vy
        temp=np.copy(self.Vy0)
        self.Vy0=np.copy(self.Vy)
        self.Vy=np.copy(temp)
        
        C=self.Vy[[29,30,31,32,33,34,35,36]][:,[29,30,31,32,33,34,35,36]]
        D=self.Vy0[[29,30,31,32,33,34,35,36]][:,[29,30,31,32,33,34,35,36]]
        

        print("BRZINE POSLE SMENE Vy")
        print(C)
        print(D)

        self.difuseVelocityY(iter)

        C=self.Vy[[29,30,31,32,33,34,35,36]][:,[29,30,31,32,33,34,35,36]]
        D=self.Vy0[[29,30,31,32,33,34,35,36]][:,[29,30,31,32,33,34,35,36]]
        

        print("BRZINA POSLE DIFUSA Vy")
        print(C)
        print(D)

        self.project(iter)

        A=self.Vx[[29,30,31,32,33,34,35,36]][:,[29,30,31,32,33,34,35,36]]
        B=self.Vx0[[29,30,31,32,33,34,35,36]][:,[29,30,31,32,33,34,35,36]]
        C=self.Vy[[29,30,31,32,33,34,35,36]][:,[29,30,31,32,33,34,35,36]]
        D=self.Vy0[[29,30,31,32,33,34,35,36]][:,[29,30,31,32,33,34,35,36]]
        print("BRZINE POSLE PROJECT")
        
        print(A)
        print(B)
        print(C)
        print(D)

        #SWAP Vx
        temp=np.copy(self.Vx0)
        self.Vx0=np.copy(self.Vx)
        self.Vx=np.copy(temp)

        #SWAP Vy
        temp=np.copy(self.Vy0)
        self.Vy0=np.copy(self.Vy)
        self.Vy=np.copy(temp)

        A=self.Vx[[29,30,31,32,33,34,35,36]][:,[29,30,31,32,33,34,35,36]]
        B=self.Vx0[[29,30,31,32,33,34,35,36]][:,[29,30,31,32,33,34,35,36]]
        C=self.Vy[[29,30,31,32,33,34,35,36]][:,[29,30,31,32,33,34,35,36]]
        D=self.Vy0[[29,30,31,32,33,34,35,36]][:,[29,30,31,32,33,34,35,36]]
        print("BRZINE POSLE SWAPA")
        
        print(A)
        print(B)
        print(C)
        print(D)

        self.advectVelocityX()
        self.advectVelocityY()

        A=self.Vx[[29,30,31,32,33,34,35,36]][:,[29,30,31,32,33,34,35,36]]
        B=self.Vx0[[29,30,31,32,33,34,35,36]][:,[29,30,31,32,33,34,35,36]]
        C=self.Vy[[29,30,31,32,33,34,35,36]][:,[29,30,31,32,33,34,35,36]]
        D=self.Vy0[[29,30,31,32,33,34,35,36]][:,[29,30,31,32,33,34,35,36]]
        print("BRZINE POSLE ADVECTA")
        
        print(A)
        print(B)
        print(C)
        print(D)


        # self.project(iter)


    def printMatrixAndValues(self):
        # print("DT= ")
        # print(self.dt)
        # print("DIFF= ")
        # print(self.diff)
        # print("VISC= ")
        # print(self.visc)

        # print("DensityZERO= ")
        # print(np.round(self.densityZERO,3))
        # print("Density= ")
        # print(np.round(self.density,3))
    
        # print("Vx= ")
        # print(self.Vx)
        # print("Vy= ")
        # print(self.Vy)
        # print("Vx0= ")
        # print(self.Vx0)
        # print("Vy0= ")
        # print(self.Vy0)

        A=self.Vx[[29,30,31,32,33,34,35,36]][:,[29,30,31,32,33,34,35,36]]
        B=self.Vx0[[29,30,31,32,33,34,35,36]][:,[29,30,31,32,33,34,35,36]]
        C=self.Vy[[29,30,31,32,33,34,35,36]][:,[29,30,31,32,33,34,35,36]]
        D=self.Vy0[[29,30,31,32,33,34,35,36]][:,[29,30,31,32,33,34,35,36]]
        print("BRZINE Vx,Vx0,Vy,Vy0")
        
        print(A)
        print(B)
        print(C)
        print(D)


