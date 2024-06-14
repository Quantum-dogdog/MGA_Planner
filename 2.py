from izzo import izzo2015
import math
import numpy as np

from jplephem.spk import SPK


##算的是天问1号的total_dv结果是4.6km/s,与nasa的trajectory browser计算结果4.57km/s十分接近，完成交叉验证



np.set_printoptions(precision=2)  #小数点后两位


kernel = SPK.open('de405.bsp')



mu_venus = 3.24859e5
mu_earth = 3.986004418e5  # [km ** 3 / s ** 2]
mu_mars = 4.28284e4
mu_jupiter = 1.26713e8
mu_saturn = 3.79406e7
mu_uranus = 5.79456e6
mu_neptune = 6.83653e6


mu_sun = 1.32712440018e11  # [km ** 3 / s ** 2]
mu_moon = 4.90280e3



r_venus = 6051.8
r_earth = 6378.14 
r_moon = 1737.4    #平均球径，其他的都是赤道球径，单位km
r_mars = 3396.19
r_jupiter = 71492
r_saturn = 60268
r_uranus = 25559
r_neptune = 24764

r_sun = 695700

fashe_date = 2459053.5             #天问1号发射日期
flytime = 202                      #days
arrive_date = fashe_date + flytime

position1 = kernel[0,3].compute(fashe_date)
position2 = kernel[0,4].compute(arrive_date)

r1 = np.array(position1)  # [km]  #太阳到地球矢量
r2 = np.array(position2)  # [km]  #太阳到火星矢量



#r1 = np.array([7.59e+07,-1.19e+08,-5.18e+07])  # [km]  #太阳到地球矢量
#r2 = np.array([1.42e+07,2.13e+08,9.71e+07])  # [km]  #太阳到火星矢量

tof = flytime * 86400  # [s]

v1, v2 =izzo2015(mu_sun, r1, r2, tof, prograde=True, low_path=True)



##print(v1) # km/s 日心坐标系速度
##print(v2) # km/s
print()



position, velocity = kernel[0,3].compute_and_differentiate(fashe_date)


velocity_per_second = velocity / 86400.0

#print('earth vx vy vx',velocity_per_second)  #单位km/s

v_infinite_earth = np.linalg.norm(v1-velocity_per_second)


r_park_earth = r_earth + 200          #200kmLEO
v_park = math.sqrt(mu_earth/r_park_earth)


#print(v_park)
v_leave_earth = math.sqrt(v_infinite_earth * v_infinite_earth +2 * mu_earth/r_park_earth)      #从停泊轨道加完速后的速度


#print(v_leave_earth)

dv1 = v_leave_earth - v_park



print('dv1',dv1,'km/s')        #单位km/s

position, velocity = kernel[0,4].compute_and_differentiate(arrive_date)  #4是火星barycenter 


velocity_per_second = velocity / 86400.0

#print('mars vx vy vx',velocity_per_second)  #单位km/s


v_infinite_mars = np.linalg.norm(v2-velocity_per_second)


r_park_mars = r_mars + 200        #200km火星停泊轨道

a_mars = - mu_mars / (v_infinite_mars * v_infinite_mars)

v_arrive_mars = math.sqrt(2 * mu_mars / r_park_mars - mu_mars / a_mars)

v_escape = math.sqrt(2 * mu_mars / r_park_mars)           #火星逃逸速度

dv2 = v_arrive_mars - v_escape



print('dv2',dv2,'km/s')             #单位km/s

total_dv = dv1 + dv2

print('total_dv',total_dv,'km/s')



