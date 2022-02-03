import numpy as np
import pygame
from fluid import Fluid



#Drawing the Fluid
def redraw(screen):             
    for i in range(1,n-1):
        k=i-1
        for j in range(1,n-1):
            m=j-1
            if fluid.densityZERO[i,j]>=255:
                color=255
            else:
                color=fluid.densityZERO[i,j]
            pygame.draw.rect(screen,(color,color,color),(k*8,m*8,8,8))
    pygame.display.update()



#Initialize the Pygame
res=512
pygame.init()
screen=pygame.display.set_mode((res,res))
pygame.display.set_caption("2D Dinamika Fluida")
icon=pygame.image.load("water-drop.png")
pygame.display.set_icon(icon)

#Initialize the Fluid

n=64+2 # +2 for border case! 
iter=4     
dt=0.01
diff=0.001
visc=0.0001
fluid=Fluid(n,dt,diff,visc)

# A=np.array([[0,1,2,3],
#     [4,5,6,7],
#     [8,9,10,11],
#     [12,13,14,15]])
# fluid.addVelocity(32,32,1,1)
# A=fluid.Vx0[[29,30,31,32,33,34,35,36]][:,[29,30,31,32,33,34,35,36]]
# B=fluid.Vx[[29,30,31,32,33,34,35,36]][:,[29,30,31,32,33,34,35,36]]

# print(A)
# print(B)

#fluid.VelStep(iter)
# fluid.VelStep(iter)
# fluid.VelStep(iter)
# fluid.VelStep(iter)
# fluid.VelStep(iter)
#fluid.printMatrixAndValues()

#Game Loop    
running=True
while running:
    screen.fill((0,0,0))
    mx,my=pygame.mouse.get_pos()        #Mouse Position
    for event in pygame.event.get():
        if event.type==pygame.QUIT:
            running =False
    
    mfcx=(mx-256)/8              #Mouse from Center,divide by 8
    mfcy=(my-256)/8              #because Matrix is scaled by 8             
    # print(mfcx)
    # print(mfcy)
    fluid.addDensVel(mfcx,mfcy)
    fluid.Step(iter)
    redraw(screen)
    



