import math
import numpy as np
from jplephem.spk import SPK
import constant


np.set_printoptions(precision=6)  #小数点后两位

kernel = SPK.open('de405.bsp')


jd = 2459053.5

mu = constant.mu_sun

r_v, v_v = kernel[0,3].compute_and_differentiate(jd)   #km km/d

v_v = v_v / 86400   #km/s

print(r_v,'\n', v_v)


t1 = 202









class Twobody:
    
   def statetoelement(self,r_vec, v_vec, mu):     #把xyz、vxvyvz转换成轨道六根数

     r = np.linalg.norm(r_vec)
     v = np.linalg.norm(v_vec)
     v_r = np.dot(r_vec / r, v_vec)
     v_p = np.sqrt(v**2 - v_r**2)

     h_vec = np.cross(r_vec, v_vec)
     h = np.linalg.norm(h_vec)             #km2/s
     i = np.arccos(h_vec[2] / h)    #弧度，另外有转换成度数正负的问题，比如文档里aop是300度，这里乘以57.3直接转换的是-50度，不过不影响用
     K = np.array((0, 0, 1))
     N_vec = np.cross(K, h_vec)
     N = np.linalg.norm(N_vec)
     if N_vec[1]<0: 
        raan = 2 * np.pi - np.arccos(N_vec[0] / N)
     else:
        raan = np.arccos(N_vec[0] / N)
     e_vec = np.cross(v_vec, h_vec) / mu - r_vec / r
     e = np.linalg.norm(e_vec)
     if e_vec[2]<0:
        
        aop = 2 * np.pi - np.arccos(np.dot(N_vec, e_vec) / (N * e))
     else:
        aop = np.arccos(np.dot(N_vec, e_vec) / (N * e))
     if v_r<0:
        
        ta = 2 * np.pi - np.arccos(np.dot(r_vec / r, e_vec / e))
     else:
        ta = np.arccos(np.dot(r_vec / r, e_vec / e))
     
     return h, i, e, raan, aop, ta
   
   def elementtostate(self, h, i, e, raan, aop, ta):
     r_w = h**2 / mu / (1 + e * np.cos(ta)) * np.array((np.cos(ta), np.sin(ta), 0))
     v_w = mu / h * np.array((-np.sin(ta), e + np.cos(ta), 0))
     R1 = np.array([[np.cos(-aop), -np.sin(-aop), 0], 
              [np.sin(-aop), np.cos(-aop), 0],
              [0, 0, 1]])

     R2 = np.array([[1, 0, 0],
              [0, np.cos(-i), -np.sin(-i)],
              [0, np.sin(-i), np.cos(-i)]])

     R3 = np.array([[np.cos(-raan), -np.sin(-raan), 0],
              [np.sin(-raan), np.cos(-raan), 0],
              [0, 0, 1]])

     r_rot = r_w @ R1 @ R2 @ R3
     v_rot = v_w @ R1 @ R2 @ R3
     return r_rot, v_rot

   def getsma(self, r_vec, v_vec, mu):    #获得半长轴
       
     r = np.linalg.norm(r_vec)
     v = np.linalg.norm(v_vec)
  
     sma = 1/(2/r - v*v/mu)
     return sma
   '''
   def eafromta(self, ta, e):       #获得偏近点角

     if e<1:
        ea = ta - e * math.sin(ta)
     else:
        ea = e * math.sin(ta) - ta
     return ea
     

   def taafterdt(self, sma, e, ta, ea, t1, mu):        #求t1时间后的真近点角  单位是天
     
     dum1 = math.sqrt((1-e)/(1+e))
     ma = 2*math.atan(dum1 * math.tan(ta/2)) - (e*math.sqrt(1-e*e)*math.sin(ta))/(1+e*math.cos(ta))
     n = math.sqrt(mu/(sma*sma*sma))
     t0 = ma/n
     t = t0 + t1 * 86400
     ma = n*t
     ta1 = 2*math.atan(1/dum1*math.tan(ea/2))
     return ta1
   '''
   
 

      
tb = Twobody()


h, i, e, raan, aop, ta = tb.statetoelement(r_v, v_v, mu)

i_du = i*57.3
raan_du = raan*57.3
aop_du = aop*57.3
ta_du = ta*57.3


print(h,'\n', i_du,'\n', e,'\n', raan_du,'\n', aop_du,'\n', ta_du)

r_b, v_b = tb.elementtostate(h, i, e, raan, aop, ta)

print(r_b,'\n', v_b)


r_2, v_2 = kernel[0,3].compute_and_differentiate(jd+t1)   #km km/d

v_2 = v_2 / 86400   #km/s

print(r_2,'\n', v_2)


















