from jplephem.spk import SPK
import numpy as np
np.set_printoptions(precision=2)  #小数点后两位


kernel = SPK.open('de405.bsp')
#print(kernel)


jd = 2457261.5  #2457061.5（2015年2月8日）是时间

position, velocity = kernel[0,4].compute_and_differentiate(jd)  #4是火星barycenter 

print('x y z', position)    #单位km

print('dx dy dz', velocity)  #每天行进的距离

velocity_per_second = velocity / 86400.0

print('vx vy vx',velocity_per_second)  #单位km/s

position = kernel[4,499].compute(2457061.5)   #火星和火星barycenter的相对距离

print('x y z', position)    #单位km


position = kernel[0,4].compute(2457061.5)
position -= kernel[0,3].compute(2457061.5)
position -= kernel[3,399].compute(2457061.5)       #火星相对地球的位置

print('x y z', position)    #单位km

print()

position, velocity = kernel[0,4].compute_and_differentiate(jd)  #4是火星barycenter 

print('mars x y z', position)    #单位km


position = kernel[0,3].compute(2457061.5)
position -= kernel[3,399].compute(2457061.5)      

print('earth x y z', position)    #单位km


position = kernel[0,3].compute(2457061.5)
      

print('earth x y z', position)    #单位km
