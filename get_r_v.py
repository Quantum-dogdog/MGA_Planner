from jplephem.spk import SPK
import numpy as np

#np.set_printoptions(precision=2)  #小数点后两位

kernel = SPK.open('de405.bsp')

def randv(shuzi,riqi):
    r,v = kernel[0,shuzi].compute_and_differentiate(riqi)   #输出日心坐标系下星球的r和v，单位是km和km/d
    v = v / 86400  #单位km/s
    return r,v
