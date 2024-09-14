import random
from get_dv import fashedv,buhuodv,dsmdv,flybyhou
from izzo import izzo2015
import numpy as np
from get_r_v import randv
from get_mu_ra import muandra
from erti import Twobody


'''
tian = [202,500]        #每条腿间飞行总时间
xishu = [0.5,0.43]       #深空机动发生时刻占整条腿飞行时间的系数
gaodu = [200,200,2000]    #轨道高度距离星球表面
jiao = [1.5]    #飞越b平面上矢量b与T_bar的夹角
agent = [202,500,0.5,0.43,200,200,2000,1.5]
'''





# 目标函数，这里使用第1、3、5元素的和的负值，以实现最小化
def objective_function(agent):
    riqi = 2459053.5      #发射日期

    xulie = [3,4,5]   #引力弹弓总序列

    v1 = [28.6,12.5,7.5]

    
    
    mu_sun,ra_sun = muandra(0)
    
    r,v = randv(xulie[0], riqi)  #发射星球的位置和速度
    mu ,ra = muandra(xulie[0])   #发射星球的mu和radius
    dv1 = fashedv(xulie[0],v1,v,agent[4])
    tb = Twobody()
    h, i, e, raan, aop, ta = tb.statetoelement(r, v1, mu_sun)
    a = tb.getsma(r, v1, mu_sun)
    ta_next = tb.taafterdt(a, e, ta, agent[0] * agent[2], mu_sun)
    r2,v2 = tb.elementtostate(h, i, e, raan, aop, ta_next)   #开普勒传播一定时间后的位置和速度，在这里施加dsm1
    r3,v5 = randv(xulie[1], riqi + agent[0])
    tof = agent[0] * (1 - agent[2]) * 86400
    v3,v4 = izzo2015(mu_sun, r2, r3, tof, prograde=True, low_path=True)
    dsm1 = dsmdv(v3, v2)
    r_leg1,v_leg1 = randv(xulie[1], riqi + agent[0])    
    v_out = flybyhou(xulie[1],v4,v_leg1,agent[5],agent[-1])
    h1, i1, e1, raan1, aop1, ta1 = tb.statetoelement(r_leg1, v_out, mu_sun)
    a1 = tb.getsma(r_leg1, v_out, mu_sun)    #半长轴
    
    ta1_next = tb.taafterdt(a1, e1, ta1, agent[1] * agent[3], mu_sun)
    r6,v6 = tb.elementtostate(h1, i1, e1, raan1, aop1, ta1_next) #开普勒传播一定时间后的位置和速度，在这里施加dsm2
    r_end,v7 = randv(xulie[2], riqi + agent[0] + agent[1])
    tof1 = agent[1] * (1 - agent[3]) * 86400
    v8,v9 = izzo2015(mu_sun, r6, r_end, tof1, prograde=True, low_path=True)
    dsm2 = dsmdv(v8, v6)
    dv2 = buhuodv(xulie[2],v9,v7,agent[6])
    zong_dv = dv1 + dsm1 + dsm2 + dv2

    return zong_dv


# 初始化种群，确保每个元素都是两位小数
def initialize_population(pop_size, ranges):
    return [[round(random.uniform(ranges[i][0], ranges[i][1]), 2) for i in range(len(ranges))] for _ in range(pop_size)]

# 差分进化算法主程序
def differential_evolution(pop_size, F, CR, generations, ranges):
    # 初始化种群
    population = initialize_population(pop_size, ranges)
    dim = len(ranges)  # 问题维度

    for generation in range(generations):
        print(f"第 {generation} 代最佳适应度: {objective_function(population[0])}")
        print(f"第 {generation} 代最佳序列: {population[0]}")

        # 生成差分向量
        idxs = [idx for idx in range(pop_size) if idx != 0]
        a, b, c = random.sample(idxs, 3)
        mutant = []
        for j in range(dim):
            mutation = round(population[a][j] + F * (population[b][j] - population[c][j]), 2)
            # 保证变异后的值在定义域内
            mutant.append(max(ranges[j][0], min(ranges[j][1], mutation)))

        # 交叉操作
        trial = []
        for j in range(dim):
            if random.random() < CR:
                trial.append(mutant[j])
            else:
                trial.append(population[0][j])

        # 选择操作
        trial_fitness = objective_function(trial)
        current_fitness = objective_function(population[0])
        if trial_fitness < current_fitness:
            population[0] = trial

    # 找出最佳解
    best_fitness = float('inf')
    best_agent = None
    for agent in population:
        fitness = objective_function(agent)
        if fitness < best_fitness:
            best_fitness = fitness
            best_agent = agent

    return best_agent

# 参数设置
pop_size = 100  # 种群大小
F = 0.5  # 差分权重
CR = 0.9  # 交叉概率
generations = 100  # 迭代次数

# 定义取值范围
ranges = [
    (100, 600),  # 前2个元素
    (100, 600),
    (0.1, 0.9),  # 3-4个元素
    (0.1, 0.9),
    (200, 2000), # 5-6个元素
    (200, 2000),
    (200, 2000),
    (0.01, 3.14) # 最后一个元素
]

# 运行差分进化算法
best_agent = differential_evolution(pop_size, F, CR, generations, ranges)
print("最佳序列:", best_agent)
