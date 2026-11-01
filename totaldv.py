from izzo import izzo2015
import math
import numpy as np
import constant
from jplephem.spk import SPK


##算的是天问1号的total_dv结果是4.6km/s,与nasa的trajectory browser计算结果4.57km/s十分接近，完成交叉验证



np.set_printoptions(precision=2)  #小数点后两位

def totaldv(a, b, mu_leave, mu_arrive, r_leave, r_arrive, fashe_date, flytime, gaodu1, gaodu2):
    kernel = SPK.open('de405.bsp')

    mu_sun = constant.mu_sun



    arrive_date = fashe_date + flytime

    position1 = kernel[0,a].compute(fashe_date)
    position2 = kernel[0,b].compute(arrive_date)

    r1 = np.array(position1)  # [km]  #太阳到地球矢量
    r2 = np.array(position2)  # [km]  #太阳到火星矢量



#r1 = np.array([7.59e+07,-1.19e+08,-5.18e+07])  # [km]  #太阳到地球矢量
#r2 = np.array([1.42e+07,2.13e+08,9.71e+07])  # [km]  #太阳到火星矢量

    tof = flytime * 86400  # [s]

    v1, v2 =izzo2015(mu_sun, r1, r2, tof, prograde=True, low_path=True)



##print(v1) # km/s 日心坐标系速度
##print(v2) # km/s
    print()
 


    position, velocity = kernel[0,a].compute_and_differentiate(fashe_date)


    velocity_per_second = velocity / 86400.0

#print('leave vx vy vx',velocity_per_second)  #单位km/s

    v_infinite_leave = np.linalg.norm(v1-velocity_per_second)


    r_park_leave = r_leave + gaodu1          #200kmLEO
    v_park = math.sqrt(mu_leave/r_park_leave)


#print(v_park)
    v_leave_leave = math.sqrt(v_infinite_leave * v_infinite_leave +2 * mu_leave/r_park_leave)      #从停泊轨道加完速后的速度


#print(v_leave_leave)

    dv1 = v_leave_leave - v_park



     #print('dv1',dv1,'km/s')        #单位km/s

    position, velocity = kernel[0,b].compute_and_differentiate(arrive_date)  #4是火星barycenter 


    velocity_per_second = velocity / 86400.0

     #print('arrive vx vy vx',velocity_per_second)  #单位km/s


    v_infinite_arrive = np.linalg.norm(v2-velocity_per_second)


    r_park_arrive = r_arrive + gaodu2        #200km火星停泊轨道

    a_arrive = - mu_arrive / (v_infinite_arrive * v_infinite_arrive)

    v_arrive_arrive = math.sqrt(2 * mu_arrive / r_park_arrive - mu_arrive / a_arrive)

    v_escape = math.sqrt(2 * mu_arrive / r_park_arrive)           #火星逃逸速度

    dv2 = v_arrive_arrive - v_escape



     #print('dv2',dv2,'km/s')             #单位km/s

    total_dv = dv1 + dv2

     #print('total_dv',total_dv,'km/s')

    return dv1, dv2, total_dv






