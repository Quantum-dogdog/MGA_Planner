import sequence
from izzo import izzo2015
from jplephem.spk import SPK
import numpy as np
import math
import constant

np.set_printoptions(precision=2)  #小数点后两位

kernel = SPK.open('de405.bsp')

jd = 2459053.5
mu_sun = 1.32712440018e11  # [km ** 3 / s ** 2]

xulie, n =sequence.sequence(3,4,5,6)  #在这里输入引力弹弓飞行序列
leg_time = [202,500,700]      #各leg之间飞行时间

i = list(range(0,n))    #0,1,2...n
j = [0]
j.extend(leg_time)
results = []
jd_next = jd
for i,j in zip(i,j):
   
  jd_next = jd_next + j    
  position = kernel[0,xulie[i]].compute(jd_next)    
  r_i= np.array(position) 
  results.append(r_i)
  print('节点的儒略日',jd_next)
print()


i = list(range(0,n))
j = list(range(1,n))
sudu =[]
for i,j in zip(i,j):

 v1,v2 = izzo2015(mu_sun, results[i], results[j], leg_time[i]*86400, prograde=True)
 sudu.append(v2)
print('兰伯特求的v2',sudu)    #单位km/s

position, velocity = kernel[0,xulie[0]].compute_and_differentiate(jd + leg_time[0]) 

print('第2个节点的日心坐标',position)    #单位km

vp = velocity / 86400.0

print('第2个节点的绕日速度vp',vp)  #单位km/s

print()


leg = [xulie[i:i+2] for i in range(len(xulie) - 1)]
print('leg[0][0]',leg[0][0])


def bplane(mu, rv, vinfvec, bmag, theta, BdotT, BdotR, status_ok):
    """
    Converts the given Fortran code to Python.
    
    Parameters:
    mu (float): central body gravity parameter (km^3/s^2)
    rv (np.ndarray): state vector (km, km/s)
    vinfvec (np.ndarray): incoming V-infinity vector (km/s)
    bmag (float): magnitude of B vector (km)
    theta (float): aim point orientation (rad)
    BdotT (float): B dot T (km)
    BdotR (float): B dot R (km)
    status_ok (bool): False if there were errors (non-hyperbolic or degenerate state)
    
    Returns:
    The resulting Python code.
    """
    
    # Local variables
    r = rv[:3]
    v = rv[3:]
    h = np.cross(r, v)
    hhat = h / np.linalg.norm(h)
    
    if np.all(hhat == 0):
        print('error: degenerate state.')
        status_ok = False
        return
    
    rdv = np.dot(r, v)
    rmag = np.linalg.norm(r)
    vmag = np.linalg.norm(v)
    vmag2 = vmag * vmag
    vinf2 = vmag2 - 2 * mu / rmag
    vinf = np.sqrt(vinf2)
    evec = np.cross(v, h) / mu - r / rmag
    e = np.linalg.norm(evec)
    
    if e <= 1:
        print('error: state is not hyperbolic.')
        status_ok = False
        return
    
    ehat = evec / e
    a = 1 / (2 / rmag - vmag2 / mu)
    hehat = np.cross(hhat, ehat)
    sd2 = 1 / e
    cd2 = np.sqrt(1 - sd2 * sd2)
    Shat = cd2 * hehat + sd2 * ehat
    That = np.cross(Shat, np.array([0, 0, 1]))
    
    if np.all(That == 0):
        print('error: vinf vector is parallel to z-axis.')
        status_ok = False
        return
    
    Rhat = np.cross(Shat, That)
    Bhat = np.cross(Shat, h)
    ct = np.dot(Bhat, That)
    st = np.dot(Bhat, Rhat)
    
    # Outputs
    bmag = abs(a) * np.sqrt(e * e - 1)
    theta = np.arctan2(st, ct)
    vinfvec = vinf * Shat
    BdotT = bmag * ct
    BdotR = bmag * st
    
    return vinfvec, bmag, theta, BdotT, BdotR, status_ok

  
r_pe = constant.r_mars + 200    #飞跃近心点

vinf_in = np.linalg.norm(sudu[0]-vp)  #待修改为v2-vp

mu = constant.mu_mars

jiao1 = 2*math.asin(1/(1+(r_pe*vinf_in*vinf_in/mu)))  #求vinf incoming和vinf outcoming的夹角，r_pe是双曲线近心点的高度，mu_mars是flyby星球的mu


print('vinf incoming和outcoming的夹角',jiao1*57.3)
