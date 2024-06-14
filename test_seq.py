import sequence

from izzo import izzo2015

from jplephem.spk import SPK

import numpy as np
np.set_printoptions(precision=2)  #小数点后两位


kernel = SPK.open('de405.bsp')

jd = 2457261.5


results = []

xulie, n =sequence.sequence(3,4,5,6)  #在这里输入引力弹弓飞行序列

i = list(range(0,n))

leg_time = [200,500,700]      #各leg之间飞行时间

j = [0]

j.extend(leg_time)

for i,j in zip(i,j):

  jd = jd + j
    
  position = kernel[0,xulie[i]].compute(jd)  
  
  r_i= np.array(position)


  #print(i)
  #print(j)
    
  #print(jd) 
  #print('r%i'% xulie[i],r_i)

  
  results.append(r_i)
#print(results)
print()
r1 = results[0]

print(r1)


mu_sun = 1.32712440018e11  # [km ** 3 / s ** 2]


i = list(range(0,n))
j = list(range(1,n))

for i,j in zip(i,j):

 v1,v2 = izzo2015(mu_sun, results[i], results[j], leg_time[i], prograde=True)

 print(v1)
 #print(i,j,i)
 


