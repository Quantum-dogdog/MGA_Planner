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
leg_time = [200,500,700]      #各leg之间飞行时间

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
print('xulie',xulie)

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
print('leg[0][1]',leg[0][1])
print('leg[1][0]',leg[1][0])
print('leg[1][1]',leg[1][1])
print('leg[2][0]',leg[2][0])
print('leg[2][1]',leg[2][1])




  
r_pe = constant.r_mars + 200    #飞跃近心点

vinf_in = np.linalg.norm(sudu[0]-vp)  #待修改为v2-vp

mu = constant.mu_mars

jiao1 = 2*math.asin(1/(1+(r_pe*vinf_in*vinf_in/mu)))  #求vinf incoming和vinf outcoming的夹角，r_pe是双曲线近心点的高度，mu_mars是flyby星球的mu


print('vinf incoming和outcoming的夹角',jiao1*57.3)





'''
上面的都是test，不用看

'''











