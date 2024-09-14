import totaldv
import constant

mu_leave = constant.mu_earth
mu_arrive = constant.mu_venus
r_leave = constant.r_earth
r_arrive = constant.r_venus

fashe_date = 2459053.5             #天问1号发射日期
flytime = 202                      #days


dv1, dv2, dv = totaldv.totaldv(3, 4, constant.mu_earth, constant.mu_mars, constant.r_earth, constant.r_mars, fashe_date, flytime, 200, 200)


print(dv1)
print(dv2)
print(dv)




dv1, dv2, dv = totaldv.totaldv(3, 2, mu_leave, mu_arrive, r_leave, r_arrive, fashe_date, flytime, 200, 200)


print(dv1)
print(dv2)
print(dv)
