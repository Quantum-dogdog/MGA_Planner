from izzo import izzo2015

import numpy as np

from jplephem.spk import SPK

np.set_printoptions(precision=2)  #小数点后两位


kernel = SPK.open('de405.bsp')

mu_earth = 3.986004418e5  # [km ** 3 / s ** 2]
r1 = np.array([15945.34, 0.0, 0.0])  # [km]
r2 = np.array([12214.83899, 10249.46731, 0.0])  # [km]
tof = 76.0 * 60  # [s]

v1, v2 =izzo2015(mu_earth, r1, r2, tof, prograde=True, low_path=True)

#expected_v1 = np.array([2.058913, 2.915965, 0.0])  # [km / s]
#expected_v2 = np.array([-3.451565, 0.910315, 0.0])  # [km / s]

print(v1)
print(v2)
print()


mu_earth = 3.986004418e5  # [km ** 3 / s ** 2]
r1 = np.array([5000.0, 10000.0, 2100.0])  # [km]
r2 = np.array([-14600.0, 2500.0, 7000.0])  # [km]
tof = 3600  # [s]

v1, v2 =izzo2015(mu_earth, r1, r2, tof, prograde=True, low_path=True)

print(v1)
print(v2)
print()
mu_sun = 39.47692641  # [AU ** 3 / year ** 2]
r1 = np.array([0.159321004, 0.579266185, 0.052359607])  # [AU]
r2 = np.array([0.057594337, 0.605750797, 0.068345246])  # [AU]
tof = 0.010794065  # [year]

v1, v2 =izzo2015(mu_sun, r1, r2, tof)

#expected_v1 = np.array([-9.303603251, 3.018641330, 1.536362143])  # [AU / year]

print(v1)
print()

mu_sun = 1.32712440018e11  # [km ** 3 / s ** 2]
r1 = np.array([-1.10e+08, 8.93e+07, 3.87e+07])  # [km]  #地球2015.2.8
r2 = np.array([-1.14e+08, 1.95e+08, 9.23e+07])  # [km]  #火星2015.2.8+200days
tof = 300 * 86400  # [s]

v1, v2 =izzo2015(mu_sun, r1, r2, tof, prograde=True, low_path=True)



print(v1) # km/s 日心坐标系速度
print(v2) # km/s
print()

jd = 2457261.5


position, velocity = kernel[0,3].compute_and_differentiate(2457061.5)


velocity_per_second = velocity / 86400.0

print('earth vx vy vx',velocity_per_second)  #单位km/s

dv1=np.linalg.norm(v1-velocity_per_second)

print(dv1)

position, velocity = kernel[0,4].compute_and_differentiate(jd)  #4是火星barycenter 


velocity_per_second = velocity / 86400.0

print('mars vx vy vx',velocity_per_second)  #单位km/s

dv2=np.linalg.norm(v2-velocity_per_second)

print(dv2)
