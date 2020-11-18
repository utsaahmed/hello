from numpy import *
import numpy as np 
import time
from matplotlib import pyplot as plt
# For graphical presentation
plt.rcParams.update({'font.size':22}) 
plt.rc('legend',fontsize=16)
plt.rcParams['figure.figsize']=(20,10)

N=int(input("Enter the number of grids in x "))
N1=N
rho=np.ones((N1,N));U=np.ones((N1,N,4));p=np.ones((N1,N));E=np.ones((N1,N));T=np.ones((N1,N));\
a=np.ones((N1,N));H=np.ones((N1,N));
F=np.ones((N1,N,4));


lmax=np.ones((N1,N));\
dt=np.ones(N);
s=np.ones((N1,N));
Mx=np.ones((N1,N));
D=np.ones((N));eae=0.6;

u=np.zeros((N1,N))
v=np.zeros((N1,N))
My=np.zeros_like(Mx);
G=np.ones_like(F)
g=1.4;r=287;CFL=0.6;cv=float(r/(g-1))
ny=ex=1
sy=wx=-1
nx=ey=sx=wy=0
t=0
L=1; Ly=1 ;dx=(L/(N-1));dy=dx
#D[-1]=L
k=np.ones((N1,N,3,3))
fa=np.zeros((N1-1,N-1,4))
fp=np.zeros((N1-1,N-1,4))
ga=np.zeros((N1-1,N-1,4))
gp=np.zeros((N1-1,N-1,4))
t0=time.time()
y= np.ones(N1)

# input
for _ in range(10):
    plt.close()
    plt.close()
    plt.close()
    plt.close()

lim=lambda a,b:min(abs(a),abs(b)) if a*b>0 else 0

def fluxs(U):
    
    U_half_lx=np.zeros_like(U);U_half_ly=np.zeros_like(U);
    U_half_rx=np.zeros_like(U);U_half_ry=np.zeros_like(U)
    F_e=np.zeros_like(U);G_n=np.zeros_like(U)
    U_halfx=np.zeros_like(U);U_halfy=np.zeros_like(U)
    F_w=np.zeros_like(U);G_s=np.zeros_like(U)
    for j in range(0,N1):
        for i in range(1,N-2):
            for k in range(4):
                U_half_lx[j,i,k]=U[j,i,k]+0.5*lim(U[j,i+1,k]-U[j,i,k],U[j,i,k]-U[j,i-1,k])
                U_half_rx[j,i,k]=U[j,i+1,k]-0.5*lim(U[j,i+2,k]-U[j,i+1,k],U[j,i+1,k]-U[j,i,k])
    for j in range(1,N1-2):
        for i in range(0,N):
            for k in range(4):
                U_half_ly[j,i,k]=U[j,i,k]+0.5*lim(U[j+1,i,k]-U[j,i,k],U[j,i,k]-U[j-1,i,k])
                U_half_ry[j,i,k]=U[j+1,i,k]-0.5*lim(U[j+2,i,k]-U[j+1,i,k],U[j+1,i,k]-U[j,i,k])
    for j in range(0,N1):
        for i in range(1,N-2):
            U_halfx[j,i,:]=0.5*(U[j,i,:]+U[j,i+1,:])
    for j in range(2,N1-2):
        for i in range(2,N-2):
            p_ex=(g-1)*(U_halfx[j,i,3]-0.5*((((U_halfx[j,i,1]**2)+(U_halfx[j,i,2])**2))/U_halfx[j,i,0]))
            p_wx=(g-1)*(U_halfx[j,i-1,3]-0.5*((((U_halfx[j,i-1,1]**2)+(U_halfx[j,i-1,2])**2))/U_halfx[j,i-1,0]))
            rho_ex=U_halfx[j,i,0]; 
            rho_wx=U_halfx[j,i-1,0]; 
            
            u_ex=U_halfx[j,i,1]/U_halfx[j,i,0]
            u_wx=U_halfx[j,i-1,1]/U_halfx[j,i-1,0]
            
            s_ex=np.sqrt(g*p_ex/(rho_ex))
            s_wx=np.sqrt(g*p_wx/(rho_wx))
            
            v_ex=U_halfx[j,i,2]/U_halfx[j,i,0]
            v_wx=U_halfx[j,i-1,2]/U_halfx[j,i-1,0]
            
            V_ex=0.5*(v_ex**2+u_ex**2)
            V_wx=0.5*((v_wx**2+u_wx**2))
            # print('vex',V_ex,'vwx',V_wx,'i',i)
            
            H_ex=V_ex+(s_ex**2)/(g-1)
            H_wx=V_wx+(s_wx**2)/(g-1)
            
            Kex=np.array([[1,1,0,1],[u_ex-s_ex*ex,u_ex,0,u_ex+s_ex*ex],[v_ex,v_ex,ex,v_ex],[H_ex-u_ex*s_ex*ex,V_ex,v_ex*ex,H_ex+u_ex*s_ex*ex]])
            
            Kexi=((g-1)/(2*s_ex**2))*np.array([[H_ex+s_ex*((u_ex*ex-s_ex)/(g-1)),-(u_ex+s_ex*ex/(g-1)),-v_ex,1],\
            [-2*H_ex+(4/(g-1))*(s_ex**2),2*u_ex,2*v_ex,-2],[-2*v_ex*ex*(s_ex**2)/(g-1),0,2*(s_ex**2)*ex/(g-1),0],\
            [H_ex-s_ex*(u_ex*ex+s_ex)/(g-1),-u_ex+s_ex*ex/(g-1),-v_ex,1]])
            
            Kwx=np.array([[1,1,0,1],[u_wx-s_wx*wx,u_wx,0,u_wx+s_wx*wx],[v_wx,v_wx,1*wx,v_wx],[H_wx-u_wx*s_wx*wx,V_wx,v_wx*wx,H_wx+u_wx*s_wx*wx]])
            
            Kwxi=((g-1)/(2*s_wx**2))*np.array([[H_wx+s_wx*((u_wx*wx-s_wx)/(g-1)),-(u_wx+s_wx*wx/(g-1)),-v_wx,1],\
            [-2*H_wx+(4/(g-1))*(s_wx**2),2*u_wx,2*v_wx,-2],[-2*v_wx*wx*(s_wx**2)/(g-1),0,2*wx*(s_wx**2)/(g-1),0],\
            [H_wx-s_wx*(u_wx*wx+s_wx)/(g-1),-u_wx+s_wx*wx/(g-1),-v_wx,1]])
            #print('kwx',Kwx,'kex',Kex,'i',i)
            lamda_pl_e=np.array([[0.5*(u_ex*ex-s_ex+abs(u_ex*ex-s_ex)),0,0,0],[0,0.5*(u_ex*ex+abs(u_ex*ex)),0,0],[0,0,0.5*(u_ex*ex+abs(u_ex*ex)),0],\
            [0,0,0,0.5*(u_ex*ex+s_ex+abs(u_ex*ex+s_ex))]])
            #print('kwxi',Kwxi,'kexi',Kexi,'i',i)
            lamda_mi_e=np.array([[0.5*(u_ex*ex-s_ex-abs(u_ex*ex-s_ex)),0,0,0],[0,0.5*(u_ex*ex-abs(u_ex*ex)),0,0],[0,0,0.5*(u_ex*ex-abs(u_ex*ex)),0],\
            [0,0,0,0.5*(u_ex*ex+s_ex-abs(u_ex*ex+s_ex))]])
            
            lamda_pl_w=np.array([[0.5*(u_wx*wx-s_wx+abs(u_wx*wx-s_wx)),0,0,0],[0,0.5*(u_wx*wx+abs(u_wx*wx)),0,0],[0,0,0.5*(u_wx*wx+abs(u_wx*wx)),0],\
            [0,0,0,0.5*(u_wx*wx+s_wx+abs(u_wx*wx+s_wx))]])
            
            lamda_mi_w=np.array([[0.5*(u_wx*wx-s_wx-abs(u_wx*wx-s_wx)),0,0,0],[0,0.5*(u_wx*wx-abs(u_wx*wx)),0,0],[0,0,0.5*(u_wx*wx-abs(u_wx*wx)),0],\
            [0,0,0,0.5*(u_wx*wx+s_wx-abs(u_wx*wx+s_wx))]])
            #print('lplw',lamda_mi_w,'\nlple',lamda_mi_e,'i',i)
            F_e[j,i,:]=np.dot(np.dot(np.dot(Kex,lamda_pl_e),Kexi),U_half_lx[j,i])+np.dot(np.dot(np.dot(Kex,lamda_mi_e),Kexi),U_half_rx[j,i])
            F_w[j,i,:]=np.dot(np.dot(np.dot(Kwx,lamda_pl_w),Kwxi),U_half_rx[j,i-1])+np.dot(np.dot(np.dot(Kwx,lamda_mi_w),Kwxi),U_half_lx[j,i-1])
            #print('fw',F_w[j,i],'\nfe',F_e[j,i-1],'i',i)
    for j in range(0,N1-1):
        for i in range(0,N):
            U_halfy[j,i]=0.5*(U[j,i,:]+U[j+1,i])
            
    for j in range(2,N1-2):
        for i in range(2,N-2):
            p_ny=(g-1)*(U_halfy[j,i,3]-0.5*((((U_halfy[j,i,1]**2)+(U_halfy[j,i,2])**2))/U_halfy[j,i,0])) 
            p_sy=(g-1)*(U_halfy[j-1,i,3]-0.5*((((U_halfy[j-1,i,1]**2)+(U_halfy[j-1,i,2])**2))/U_halfy[j-1,i,0])) 
            
            rho_ny=U_halfy[j,i,0]; 
            rho_sy=U_halfy[j-1,i,0]; 
            
            u_ny=U_halfy[j,i,1]/U_halfy[j,i,0]
            u_sy=U_halfy[j-1,i,1]/U_halfy[j-1,i,0]
            
            s_ny=np.sqrt(g*p_ny/(rho_ny))
            s_sy=np.sqrt(g*p_sy/(rho_sy))
            
            v_ny=U_halfy[j,i,2]/U_halfy[j,i,0]
            v_sy=U_halfy[j-1,i,2]/U_halfy[j-1,i,0]
            
            V_ny=0.5*(v_ny**2+u_ny**2)
            V_sy=0.5*(v_sy**2+u_sy**2)
            
            H_ny=V_ny+s_ny**2/(g-1)
            H_sy=V_sy+s_sy**2/(g-1)
            
            Kny=np.array([[1,0,1,1],[u_ny,1*ny,u_ny,u_ny],[v_ny-s_ny*ny,0,v_ny,v_ny+s_ny*ny],[H_ny-v_ny*s_ny*ny,u_ny*ny,V_ny,H_ny+v_ny*s_ny*ny]])
            
            Knyi=((g-1)/(2*s_ny**2))*np.array([[H_ny+s_ny*((v_ny*ny-s_ny)/(g-1)),-u_ny,-(v_ny+s_ny*ny/(g-1)),1],\
            [-2*u_ny*ny*(s_ny**2)/(g-1),2*ny*(s_ny**2)/(g-1),0,0],[-2*H_ny+(4/(g-1))*(s_ny**2),2*u_ny,2*v_ny,-2],\
            [H_ny-s_ny*(v_ny*ny+s_ny)/(g-1),-u_ny,-v_ny+s_ny*ny/(g-1),1]])
            
            Ksy=np.array([[1,0,1,1],[u_sy,1*sy,u_ny,u_sy],[v_sy-s_sy*sy,0,v_sy,v_sy+s_sy*sy],[H_sy-v_sy*s_sy*sy,u_sy*sy,V_sy,H_sy+v_sy*s_sy*sy]])
            
            Ksyi=((g-1)/(2*s_sy**2))*np.array([[H_sy+s_sy*((v_sy*sy-s_sy)/(g-1)),-u_sy,-(v_sy+s_sy*sy/(g-1)),1],\
            [-2*u_sy*sy*(s_sy**2)/(g-1),2*sy*(s_sy**2)/(g-1),0,0],[-2*H_sy+(4/(g-1))*(s_sy**2),2*u_sy,2*v_sy,-2],\
            [H_sy-s_sy*(v_sy*sy+s_sy)/(g-1),-u_sy,-v_sy+s_sy*sy/(g-1),1]])
            
            lamda_pl_n=np.array([[0.5*(v_ny*ny-s_ny+abs(v_ny*ny-s_ny)),0,0,0],[0,0.5*(v_ny*ny+abs(v_ny*ny)),0,0],[0,0,0.5*(v_ny*ny+abs(v_ny*ny)),0],\
            [0,0,0,0.5*(v_ny*ny+s_ny+abs(v_ny*ny+s_ny))]])
            
            lamda_mi_n=np.array([[0.5*(v_ny*ny-s_ny-abs(v_ny*ny-s_ny)),0,0,0],[0,0.5*(v_ny*ny-abs(v_ny)),0,0],[0,0,0.5*(v_ny*ny-abs(v_ny*ny)),0],\
            [0,0,0,0.5*(v_ny*ny+s_ny-abs(v_ny*ny+s_ny))]])
            
            lamda_pl_s=np.array([[0.5*(v_sy*sy-s_sy+abs(v_sy*sy-s_sy)),0,0,0],[0,0.5*(v_sy*sy+abs(v_sy*sy)),0,0],[0,0,0.5*(v_sy*sy+abs(v_sy*sy)),0],\
            [0,0,0,0.5*(v_sy*sy+s_sy+abs(v_sy*sy+s_sy))]])
            
            lamda_mi_s=np.array([[0.5*(v_sy*sy-s_sy-abs(v_sy*sy-s_sy)),0,0,0],[0,0.5*(v_sy*sy-abs(v_sy*sy)),0,0],[0,0,0.5*(v_sy*sy-abs(v_sy*sy)),0],\
            [0,0,0,0.5*(v_sy*sy+s_sy-abs(v_sy*sy+s_sy))]])
            
            G_n[j,i]=np.dot(np.dot(np.dot(Kny,lamda_pl_n),Knyi),U_half_ly[j,i])+np.dot(np.dot(np.dot(Kny,lamda_mi_n),Knyi),U_half_ry[j,i])
            G_s[j,i]=np.dot(np.dot(np.dot(Ksy,lamda_pl_s),Ksyi),U_half_ry[j-1,i])+np.dot(np.dot(np.dot(Ksy,lamda_mi_s),Ksyi),U_half_ly[j-1,i])

            #print('gn2',G_n[3:-3,3:-3,2])
    return (F_e, F_w, G_n, G_s)
            

def SW_2ndorder(U,u,v,s):
    Utemp=np.zeros_like(U)
    
    lm=np.max(np.max(np.sqrt(u[:,:]**2+v[:,:]**2)+s[:,:],axis=1))
    dtm= ((CFL*dx)/lm)  
    
    pred=np.zeros_like(U)
    (F_e,F_w,G_n,G_s)=fluxs(U)
    
    
    for j in range(2,N1-2):
        for i in range(2,N-2):
            pred[j,i,:]=U[j,i,:]-(dtm/(2*dx))*(F_e[j,i,:]+F_w[j,i,:])-(dtm/(2*dy))*(G_n[j,i,:]+G_s[j,i,:])
            
    
    pred[-2,:,:]=pred[-3,:,:];
    pred[-1,:,:]=pred[-2,:,:]
    
    # U[2,:,:]=U[3,:,:];
    pred[1,:,:]=pred[2,:,:];
    pred[0,:,:]=pred[1,:,:]
    
    # U[:,-3,:]=U[:,-4,:];
    pred[:,-2,:]=pred[:,-3,:];
    pred[:,-1,:]=pred[:,-2,:]
    
    # U[:,2,:]=U[:,3,:];
    pred[:,1,:]=pred[:,2,:];
    pred[:,0,:]=pred[:,1,:]
    
    (F_pe,F_pw,G_pn,G_ps)=fluxs(pred)
    
    for j in range(2,N1-2):
        for i in range(2,N-2):
            Utemp[j,i,:]=U[j,i,:]-(dtm/(dx))*(F_pe[j,i,:]+F_pw[j,i,:])-(dtm/(dy))*(G_pn[j,i,:]+G_ps[j,i,:])
        
    
    # print('gn',G_n[3:-3,3:-3,2])
    #print(fluxs(U)[2])
    # print('gs',G_s[3:-3,3:-3,2])
    # print(G_n[3:-3,3:-3,2]+G_s[3:-3,3:-3,2])
    # 
    # print('fe',F_e[3:-3,3:-3,3])
    # print('fw',F_w[3:-3,3:-3,3])
    # print(F_e[3:-3,3:-3,3]+F_w[3:-3,3:-3,3])
    
    # print('lamdapie',lamda_pl_e[])
    # print('lamdapiw',lamda_pl_w)
    
    #print('u half',U_halfx[3:-3,2])
            
    return Utemp,dtm
            

for i in range(0,N):
    if i==0:
        D[i]=0
        
    else:
        D[i]=(i)*dx
        
for i in range(N1):
    if i==0:
        y[i]=0
    else:
        y[i]=(i)*dy
for i in range(N):    
    if (D[i]<=L/2):
        u[:,i]=0
        p[:,i]=101325*5;
        T[:,i]=300
    else:
        u[:,i]=0
        p[:,i]=101325;
        T[:,i]=300
# for i in range(N):    
#     if (D[i]<=L/2):
#         u[:,i]=0
#         p[:,i]=101325;
#         T[:,i]=300
#     else:
#         u[:,i]=0
#         p[:,i]=101325*5;
#         T[:,i]=300
    
    
#print(T);
#print(p)

for i in range(0,N):
    
    rho[:,i]=p[:,i]/(T[:,i]*r);
    s[:,i]=np.sqrt(g*r*T[:,i]);
    Mx[:,i]=u[:,i]/s[:,i];
    My[:,i]=v[:,i]/s[:,i];
    E[:,i]=rho[:,i]*((cv*T[:,i])+0.5*((u[:,i]**2)+(v[:,i]**2)));
    H[:,i]=(E[:,i]+p[:,i])/rho[:,i]
    
       
    
    U[:,i,0]=rho[:,i];
    U[:,i,1]=rho[:,i]*u[:,i];
    U[:,i,2]=rho[:,i]*v[:,i]
    U[:,i,3]=E[:,i];
#print(np.min(np.min(E,axis=1)))


#Execution part
xime = 75*(10**-5)

while(t< xime or t==xime):
    
    # for i in range(N): 
    #         U[:,i,0]=rho[:,i];
    #         U[:,i,1]=rho[:,i]*u[:,i];
    #         U[:,i,2]=rho[:,i]*v[:,i]
    #         U[:,i,3]=E[:,i];

    # for j in range(N1):   
    #     for i in range(N): 
            # if i==0 or i==N-1:
    # U[:,-1,:]=U[:,-2,:];
    # U[:,0,:]=U[:,1,:];
    #         # if j==0 or j==N1-1:
    # U[-1,:,:]=U[-2,:,:];
    # U[0,:,:]=U[1,:,:];
                
    
    
    
    
    #print("the time ",dtm)
    
    
    # for i in range(N):
    p[:,:]=(g-1)*(U[:,:,3]-0.5*((((U[:,:,1]**2)+(U[:,:,2])**2))/U[:,:,0]))     
    rho[:,:]=U[:,:,0];
    u[:,:]=U[:,:,1]/U[:,:,0];
    v[:,:]=U[:,:,2]/U[:,:,0]
    E[:,:]=U[:,:,3]
    T[:,:]=p[:,:]/(rho[:,:]*r);
    s[:,:]=np.sqrt(g*r*T[:,:]);
    
    
    F[:,:,0]=U[:,:,1];
    F[:,:,1]=((U[:,:,1]**2)/U[:,:,0])+p[:,:];
    F[:,:,2]=(U[:,:,2]*U[:,:,1])/(U[:,:,0])
    F[:,:,3]=(U[:,:,1]/U[:,:,0])*(U[:,:,3]+p[:,:])
    # for j in range(N1):
    G[:,:,0]=U[:,:,2]
    G[:,:,1]=(U[:,:,1]*U[:,:,2])/(U[:,:,0])
    G[:,:,2]=((U[:,:,2]**2)/U[:,:,0])+p[:,:];
    G[:,:,3]=(U[:,:,3]+p[:,:])*(U[:,:,2]/U[:,:,0])
    # Boundary condition
    
                
    #The scheme
    Utemp,dtm= SW_2ndorder(U,u,v,s)
    
    U=Utemp
            # if j==0 or j==N1-1:
    # U[-3,:,:]=U[-4,:,:];
    U[-2,:,:]=U[-3,:,:];
    U[-1,:,:]=U[-2,:,:]
    
    # U[2,:,:]=U[3,:,:];
    U[1,:,:]=U[2,:,:];
    U[0,:,:]=U[1,:,:]
    
    # U[:,-3,:]=U[:,-4,:];
    U[:,-2,:]=U[:,-3,:];
    U[:,-1,:]=U[:,-2,:]
    
    # U[:,2,:]=U[:,3,:];
    U[:,1,:]=U[:,2,:];
    U[:,0,:]=U[:,1,:]
    # for i in range(0,N):        
    p[:,:]=(g-1)*(U[:,:,3]-0.5*(((((U[:,:,1])**2)+((U[:,:,2])**2))/U[:,:,0])))     
    rho[:,:]=U[:,:,0];
    u[:,:]=U[:,:,1]/U[:,:,0];
    v[:,:]=U[:,:,2]/U[:,:,0]
    E[:,:]=U[:,:,3]
    
    T[:,:]=p[:,:]/(rho[:,:]*r);
    for j in range(N1):
        for i in range(N):
            if T[j,i]<0:
                print(j,i)
                raise Exception('Arbit!')
    s[:,:]=np.sqrt(g*r*T[:,:]);
    Mx[:,:]=u[:,:]/s[:,:];
    My[:,:]=v[:,:]/s[:,:];
        
    
    
    t=t+dtm
    print("total time",CFL*(dx/dtm))
    
xlist=np.linspace(0,1,N)
ylist=np.linspace(0,1,N)
x_,y_=np.meshgrid(xlist,ylist)

Tn=np.zeros(N);un=np.zeros(N); pn=np.zeros(N) ;Mn=np.zeros(N);rhon=np.zeros(N) 
for i in range(0,N):
    # Tn[i]=T[int(N1/2),i]/np.max(T[int(N1/2)])
    # un[i]=u[int(N1/2),i]/-np.max(u[int(N1/2)])
    # pn[i]=p[int(N1/2),i]/np.max(p[int(N1/2)])
    # Mn[i]=Mx[int(N1/2),i]/-np.max(Mx[int(N1/2)])
    # rhon[i]=rho[int(N1/2),i]/np.max(rho[int(N1/2)])  
    Tn[i]=T[int(N1/2),i]/np.min(T[int(N1/2)])
    un[i]=u[int(N1/2),i]/-np.min(u[int(N1/2)])
    pn[i]=p[int(N1/2),i]/np.min(p[int(N1/2)])
    Mn[i]=Mx[int(N1/2),i]/-np.min(Mx[int(N1/2)])
    rhon[i]=rho[int(N1/2),i]/np.min(rho[int(N1/2)])  

t1=time.time()    
# print("Temperatre",T)
# print("Distance",D)
# print("Mach number",Mx)
print("Execution time",(t1-t0))
plt.figure(1)

plt.plot(D,Tn,label='grid points='+str(N1))
plt.xlabel('X');plt.ylabel('$T_{n}$')
plt.title("Normalised Temperature plot")
plt.legend()
plt.figure(2)
plt.plot(D,-Mn,label='grid points='+str(N1))
plt.xlabel('X');plt.ylabel('$M_{n}$')
plt.title("Normalised Mach number plot")
#plt.legend([str(N1)])
plt.figure(3)
plt.plot(D,-un,label='grid points='+str(N1))
plt.xlabel('X');plt.ylabel('$u_{n}$')
plt.title("Normalised Velocity plot")
#plt.legend([str(N1)])
plt.figure(4)
plt.plot(D,pn,label='grid points='+str(N1))
plt.xlabel('X');plt.ylabel('$P_{n}$')
plt.title("Normalised Pressure plot")
#plt.legend([str(N1)])
plt.figure(5)
plt.plot(D,rhon,label='grid points='+str(N1))
plt.xlabel('X');plt.ylabel('$Rho_{n}$')
plt.title("Normalised Density plot")
plt.legend([str(N1)])
plt.show()
# plt.figure(6)
# plt.contour(x_,y_,np.sqrt(u**2+v**2))
# plt.colorbar()
# plt.title('Velocity contour plot')
# plt.xlabel('X');plt.ylabel('Y')
# 
# plt.figure(7)
# plt.contour(x_,y_,rho)
# plt.colorbar()
# plt.title('Density contour plot')
# plt.xlabel('X');plt.ylabel('Y')
# 
# plt.figure(8)
# plt.contour(x_[1:-1,1:-1],y_[1:-1,1:-1],p[1:-1,1:-1])
# plt.colorbar()
# plt.title('Pressure contour plot')
# plt.xlabel('X');plt.ylabel('Y')

plt.show()
