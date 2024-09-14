from get_dv import fashedv,buhuodv,dsmdv,flybyhou
from izzo import izzo2015
import numpy as np
from get_r_v import randv
from get_mu_ra import muandra
from erti import Twobody
import threading
import time
'''

def timer():
          str1 = '\n温馨提示：\n如果你还没有看到结果，那么不用等了，你输入的这个序列参数不可行；如果你看到结果了，那就当我没说'
          panduan = 0
          start = time.time()
          time.sleep(3)
          end = time.time()
          if end - start >2:
             print(str1)
t2 = threading.Thread(target = timer)
t2.start()
'''
    








def mga_dsm(xulie,riqi,  v1, tian, xishu, gaodu, jiao):
    
    zifuchuan = 'km/s'
    
   
    
    if len(xulie) == 2:
       
          jieguo = mga_dsm_2(riqi,xulie, v1, tian, xishu, gaodu, jiao)
          jieguo = str(jieguo) + zifuchuan
       
    if len(xulie) == 3: 
       
          jieguo = mga_dsm_3(riqi,xulie, v1, tian, xishu, gaodu, jiao)
          jieguo = str(jieguo) + zifuchuan
       
    if len(xulie) == 4:
       
          jieguo = mga_dsm_4(riqi,xulie, v1, tian, xishu, gaodu, jiao)
          jieguo = str(jieguo) + zifuchuan
       
    if len(xulie) == 5:
       
          jieguo = mga_dsm_5(riqi,xulie, v1, tian, xishu, gaodu, jiao)
          jieguo = str(jieguo) + zifuchuan
       
    if len(xulie) == 1: 
          jieguo = '\n初始化失败!\n序列至少应包含1个出发星球，1个目标星球，即序列长度至少为2!'
    if len(xulie) > 5: 
          jieguo = '\n初始化失败!\n序列长度过长(>5)，本程序不予计算!' 
    
     
    return jieguo









def mga_dsm_2(riqi, xulie, v1, tian, xishu, gaodu, jiao):

    mu_sun,ra_sun = muandra(0)
    
    r,v = randv(xulie[0], riqi)  #发射星球的位置和速度
    mu ,ra = muandra(xulie[0])   #发射星球的mu和radius
    dv1 = fashedv(xulie[0],v1,v,gaodu[0])
    tb = Twobody()
    h, i, e, raan, aop, ta = tb.statetoelement(r, v1, mu_sun)
    a = tb.getsma(r, v1, mu_sun)
    ta_next = tb.taafterdt(a, e, ta, tian[0] * xishu[0], mu_sun)
    r2,v2 = tb.elementtostate(h, i, e, raan, aop, ta_next)   #开普勒传播一定时间后的位置和速度，在这里施加dsm1
    r3,v5 = randv(xulie[1], riqi + tian[0])
    tof = tian[0] * (1 - xishu[0]) * 86400
    v3,v4 = izzo2015(mu_sun, r2, r3, tof, prograde=True, low_path=True)
    dsm1 = dsmdv(v3, v2)
        

    dv2 = buhuodv(xulie[1],v4,v5,gaodu[1])
    zong_dv = dv1 + dsm1 + dv2
    
    return zong_dv



def mga_dsm_3(riqi, xulie, v1, tian, xishu, gaodu, jiao):

    mu_sun,ra_sun = muandra(0)
    
    r,v = randv(xulie[0], riqi)  #发射星球的位置和速度
    mu ,ra = muandra(xulie[0])   #发射星球的mu和radius
    dv1 = fashedv(xulie[0],v1,v,gaodu[0])
    tb = Twobody()
    h, i, e, raan, aop, ta = tb.statetoelement(r, v1, mu_sun)
    a = tb.getsma(r, v1, mu_sun)
    ta_next = tb.taafterdt(a, e, ta, tian[0] * xishu[0], mu_sun)
    r2,v2 = tb.elementtostate(h, i, e, raan, aop, ta_next)   #开普勒传播一定时间后的位置和速度，在这里施加dsm1
    r3,v5 = randv(xulie[1], riqi + tian[0])
    tof = tian[0] * (1 - xishu[0]) * 86400
    v3,v4 = izzo2015(mu_sun, r2, r3, tof, prograde=True, low_path=True)
    dsm1 = dsmdv(v3, v2)
    r_leg1,v_leg1 = randv(xulie[1], riqi + tian[0])    
    v_out = flybyhou(xulie[1],v4,v_leg1,gaodu[1],jiao[0])
    h1, i1, e1, raan1, aop1, ta1 = tb.statetoelement(r_leg1, v_out, mu_sun)
    a1 = tb.getsma(r_leg1, v_out, mu_sun)    #半长轴
    
    ta1_next = tb.taafterdt(a1, e1, ta1, tian[1] * xishu[1], mu_sun)
    r6,v6 = tb.elementtostate(h1, i1, e1, raan1, aop1, ta1_next) #开普勒传播一定时间后的位置和速度，在这里施加dsm2
    r_end,v7 = randv(xulie[2], riqi + tian[0] + tian[1])
    tof1 = tian[1] * (1 - xishu[1]) * 86400
    v8,v9 = izzo2015(mu_sun, r6, r_end, tof1, prograde=True, low_path=True)
    dsm2 = dsmdv(v8, v6)
    dv2 = buhuodv(xulie[2],v9,v7,gaodu[2])
    zong_dv = dv1 + dsm1 + dsm2 + dv2
    
    return zong_dv

def mga_dsm_4(riqi, xulie, v1, tian, xishu, gaodu, jiao):

    mu_sun,ra_sun = muandra(0)
    
    r,v = randv(xulie[0], riqi)  #发射星球的位置和速度
    mu ,ra = muandra(xulie[0])   #发射星球的mu和radius
    dv1 = fashedv(xulie[0],v1,v,gaodu[0])
    tb = Twobody()
    h, i, e, raan, aop, ta = tb.statetoelement(r, v1, mu_sun)
    a = tb.getsma(r, v1, mu_sun)
    ta_next = tb.taafterdt(a, e, ta, tian[0] * xishu[0], mu_sun)
    r2,v2 = tb.elementtostate(h, i, e, raan, aop, ta_next)   #开普勒传播一定时间后的位置和速度，在这里施加dsm1
    r3,v5 = randv(xulie[1], riqi + tian[0])
    tof = tian[0] * (1 - xishu[0]) * 86400
    v3,v4 = izzo2015(mu_sun, r2, r3, tof, prograde=True, low_path=True)
    dsm1 = dsmdv(v3, v2)
    r_leg1,v_leg1 = randv(xulie[1], riqi + tian[0])    
    v_out = flybyhou(xulie[1],v4,v_leg1,gaodu[1],jiao[0])
    h1, i1, e1, raan1, aop1, ta1 = tb.statetoelement(r_leg1, v_out, mu_sun)
    a1 = tb.getsma(r_leg1, v_out, mu_sun)    #半长轴
    
    ta1_next = tb.taafterdt(a1, e1, ta1, tian[1] * xishu[1], mu_sun)
    r6,v6 = tb.elementtostate(h1, i1, e1, raan1, aop1, ta1_next) #开普勒传播一定时间后的位置和速度，在这里施加dsm2
    r_end,v7 = randv(xulie[2], riqi + tian[0] + tian[1])
    tof1 = tian[1] * (1 - xishu[1]) * 86400
    v8,v9 = izzo2015(mu_sun, r6, r_end, tof1, prograde=True, low_path=True)
    dsm2 = dsmdv(v8, v6)

    r_leg2,v_leg2 = randv(xulie[2], riqi + tian[0]+ tian[1])    
    v_out2 = flybyhou(xulie[2],v9,v_leg2,gaodu[2],jiao[1])
    h2, i2, e2, raan2, aop2, ta2 = tb.statetoelement(r_leg2, v_out2, mu_sun)
    a2 = tb.getsma(r_leg2, v_out2, mu_sun)    #半长轴
    
    ta2_next = tb.taafterdt(a2, e2, ta2, tian[2] * xishu[2], mu_sun)
    r10,v10 = tb.elementtostate(h2, i2, e2, raan2, aop2, ta2_next) #开普勒传播一定时间后的位置和速度，在这里施加dsm3
    r_end,v11 = randv(xulie[3], riqi + tian[0] + tian[1]+ tian[2])
    tof2 = tian[2] * (1 - xishu[2]) * 86400
    v12,v13 = izzo2015(mu_sun, r10, r_end, tof2, prograde=True, low_path=True)
    dsm3 = dsmdv(v12, v10)

 
    dv2 = buhuodv(xulie[3],v13,v11,gaodu[3])
   
    zong_dv = dv1 + dsm1 + dsm2 + dsm3 + dv2
    
    return zong_dv


def mga_dsm_5(riqi, xulie, v1, tian, xishu, gaodu, jiao):
    
    # 尝试执行计算，如果超时则抛出异常
    try:
        
            # 假设这里是你想要执行的耗时操作   
            mu_sun,ra_sun = muandra(0)
    
            r,v = randv(xulie[0], riqi)  #发射星球的位置和速度
            mu ,ra = muandra(xulie[0])   #发射星球的mu和radius
            dv1 = fashedv(xulie[0],v1,v,gaodu[0])
            tb = Twobody()
            h, i, e, raan, aop, ta = tb.statetoelement(r, v1, mu_sun)
            a = tb.getsma(r, v1, mu_sun)
            ta_next = tb.taafterdt(a, e, ta, tian[0] * xishu[0], mu_sun)
            r2,v2 = tb.elementtostate(h, i, e, raan, aop, ta_next)   #开普勒传播一定时间后的位置和速度，在这里施加dsm1
            r3,v5 = randv(xulie[1], riqi + tian[0])
            tof = tian[0] * (1 - xishu[0]) * 86400
            v3,v4 = izzo2015(mu_sun, r2, r3, tof, prograde=True, low_path=True)
            dsm1 = dsmdv(v3, v2)
            r_leg1,v_leg1 = randv(xulie[1], riqi + tian[0])    
            v_out = flybyhou(xulie[1],v4,v_leg1,gaodu[1],jiao[0])
            h1, i1, e1, raan1, aop1, ta1 = tb.statetoelement(r_leg1, v_out, mu_sun)
            a1 = tb.getsma(r_leg1, v_out, mu_sun)    #半长轴
    
            ta1_next = tb.taafterdt(a1, e1, ta1, tian[1] * xishu[1], mu_sun)
            r6,v6 = tb.elementtostate(h1, i1, e1, raan1, aop1, ta1_next) #开普勒传播一定时间后的位置和速度，在这里施加dsm2
            r_end,v7 = randv(xulie[2], riqi + tian[0] + tian[1])
            tof1 = tian[1] * (1 - xishu[1]) * 86400
            v8,v9 = izzo2015(mu_sun, r6, r_end, tof1, prograde=True, low_path=True)
            dsm2 = dsmdv(v8, v6)

            r_leg2,v_leg2 = randv(xulie[2], riqi + tian[0]+ tian[1])    
            v_out2 = flybyhou(xulie[2],v9,v_leg2,gaodu[2],jiao[1])
            h2, i2, e2, raan2, aop2, ta2 = tb.statetoelement(r_leg2, v_out2, mu_sun)
            a2 = tb.getsma(r_leg2, v_out2, mu_sun)    #半长轴
    
            ta2_next = tb.taafterdt(a2, e2, ta2, tian[2] * xishu[2], mu_sun)
            r10,v10 = tb.elementtostate(h2, i2, e2, raan2, aop2, ta2_next) #开普勒传播一定时间后的位置和速度，在这里施加dsm3
            r_end,v11 = randv(xulie[3], riqi + tian[0] + tian[1]+ tian[2])
            tof2 = tian[2] * (1 - xishu[2]) * 86400
            v12,v13 = izzo2015(mu_sun, r10, r_end, tof2, prograde=True, low_path=True)
            dsm3 = dsmdv(v12, v10)


            r_leg3,v_leg3 = randv(xulie[3], riqi + tian[0]+ tian[1]+ tian[2])    
            v_out3 = flybyhou(xulie[3],v13,v_leg3,gaodu[3],jiao[2])
            h3, i3, e3, raan3, aop3, ta3 = tb.statetoelement(r_leg3, v_out3, mu_sun)
            a3 = tb.getsma(r_leg3, v_out3, mu_sun)    #半长轴
    
            ta3_next = tb.taafterdt(a3, e3, ta3, tian[3] * xishu[3], mu_sun)
            r14,v14 = tb.elementtostate(h3, i3, e3, raan3, aop3, ta3_next) #开普勒传播一定时间后的位置和速度，在这里施加dsm4
            r_end,v15 = randv(xulie[4], riqi + tian[0] + tian[1]+ tian[2]+ tian[3])
            tof3 = tian[3] * (1 - xishu[3]) * 86400
            v16,v17 = izzo2015(mu_sun, r14, r_end, tof3, prograde=True, low_path=True)
            dsm4 = dsmdv(v16, v14)


            dv2 = buhuodv(xulie[4],v17,v15,gaodu[4])
   
            zong_dv = dv1 + dsm1 + dsm2 + dsm3 + dsm4 + dv2

    except:
        return 999      #有待后续修正，因为现在它到不了抛出错误的这一步就一直卡着不出结果，风扇狂转，猜测是函数里边的问题，就是当序列不可行的时候，它不会不算
    return zong_dv






 

