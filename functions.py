import numpy as np

#FUNCTION FOR THE DIFFUSION OF MATRIXES
def difuse(b,x,x0,diff,dt,iter,n):      
    a = dt * diff * (n-2) * (n-2)   
    lin_solve(b,x,x0,a,1+4*a,iter,n)
    return


#SOLVES EVERY SINGLE POINT IN 2D MATRIX
#USING THE GAUSS-SEIDEL RELAXATION METHOD
def lin_solve(b,x,x0,a,c,iter,n):       
    cRecip=1.0/c                        
    for k in range(iter):
        for i in range(1,n-1):
            for j in range(1,n-1):
                x[i,j] = (x0[i,j] + a*(x[i-1,j]+x[i+1,j]+x[i,j-1]+x[i,j+1]))*cRecip
        set_bnd(b,x,n)
    return


#SETS THE BOUNDARY CONDITION FOR THE SOLVER
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
             
    #for j in range(1,n-1):  #NOT NECESSARY
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



#FUNCTION USED TO PROJECT THE FLOW AND GET
#THE DIVERGENCE FREE(CURL ONLY) PART
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


#FUNCTION FOR ADVECTING THE FLUID OR DYE
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

            d[i,j] = s0*(t0*d0[i0i,j0i]+t1*d0[i0i,j1i])+ s1*(t0*d0[i1i,j0i]+t1*d0[i1i,j1i])

    set_bnd (b,d,n)
    return