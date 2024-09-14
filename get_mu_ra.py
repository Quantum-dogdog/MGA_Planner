import math


def muandra(shuzi):

    mu = [1.32712440018e11,0,3.24859e5,3.986004418e5,4.28284e4,1.26713e8,3.79406e7,5.79456e6,6.83653e6,4.90280e3]
    ra = [695700,0,6051.8,6378.14,3396.19,71492,60268,25559,24764,1737.4]
    return mu[shuzi], ra[shuzi]
'''
0-sun
1-null   引力弹弓用不到水星，不写了
2-venus
3-earth
4-mars
5-jupiter
6-saturn
7-uranus
8-neptune
9-moon

mu是GM    ra是radius

单位分别是：[km ** 3 / s ** 2]  km

只有月球是平均球径，其他的都是赤道球径，单位km

'''
