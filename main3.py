import random
import math
import json        #python自带的，用来把结果写成step6能吃的nodes.json，不是新装的包
import datetime    #python自带的，把儒略日翻译成几月几号
from get_dv import fashedv, buhuodv     #这一版中间一次发动机都不打，所以dsmdv、flybyhou用不上
from izzo import izzo2015
import numpy as np
from get_r_v import randv
from get_mu_ra import muandra


'''
这一版是跟 github 上 book_scripts_for_GMAT_NASA/step8_纯MGA借力优化.py 学的，
跟 main2.py 的最大差别是：**一次发动机都不打**（纯借力、无深空机动），
染色体只剩 L+1 个数：

    第0个        发射日期
    接下来 L 个  每条腿飞几天

为什么能这么少？以前 main2.py 里"飞越多高""B面转到哪一圈"都当变量去搜，
其实这两样根本不用搜：前后两段路的方向一定下来，
要转多少角（δ）、贴多高（rp）就全定了，
    rp = μ / v_inf^2 * ( 1 / sin(δ/2) - 1 )
它是**算出来的输出量**，不是决策变量。（step8 的原话：rp 是输出量。）

代价是：前后两段各自用 Lambert 解出来以后，
进超速和出超速的**模长一般对不上**，转角落不进那颗星球能给的范围内，
这种"对不上多少"就叫缺陷 Δv_fb：
    δmax = 2*asin( 1 / (1 + rp_min*v_inf^2/μ) )        #这颗星球最多能转这么多
    把 vin 在 (vin,vout) 平面里朝 vout 转 min(δ, δmax)，模长不变，
    转完还差多少 = |vout - 转过的vin| = 缺陷
缺陷当罚函数加进目标里（权重很大），压到 0 就是真能弹道打出来的轨道。

目标函数（step8 的写法，MGA 惯例，只按超速算，不含发射损失和捕获）：
    J = 权重*|v_inf出发| + 权重*|v_inf到达| + W * Σ 缺陷^2
所以 J 不是总 Δv，别拿它跟 main2.py 的总 dv 直接比。
本文件额外用他 repo 里的 fashedv / buhuodv 把"从停泊轨道发射"和"木星捕获"
这两笔也报出来，方便跟 main2.py 对照。

还有一个坑（step8 注释里写得很清楚）：两点 Lambert 在转移角≈0° 时是退化的
（谐振返回那种绕一圈再撞回同一颗星球，过同一点、同 a 的椭圆有无穷多条），
所以这种腿 Lambert 怎么解缺陷都压不到 0。这一版先不管谐振段，
分支（绕几圈、顺行逆行、高路低路）用"束搜索"全试一遍，见下面 jie_fenzhi。
'''


#捕获以后要进的那条轨道，远点是星球半径的多少倍。任务要求，不是设计变量。
yuandian = 20

#停泊轨道高度（发射和捕获都用它）
jiatuo = 200
buhuo_gaodu = 200

#飞越缺陷的罚权重，step8 默认100，扫窗口的时候用2000。跟它一样，越大越先把弹道压对。
W_feiyue = 2000.0

#出发/到达超速在目标里的权重
quan_xian_bu = 1.0
quan_xian_da = 1.0

#每条腿最多允许绕太阳几圈（0=不绕，1=绕一圈……）。圈数不进染色体，直接全试。
M_zui_da = 1

#分支组合太多时，每条腿只留着最便宜的这么几种往下接（束宽）
liu_shao = 20

#火箭能给的最大超速 km/s（C3 上限的另一种写法）。写 None 就是不管。
vinf_shang = None

#每颗星球最低飞越高度 km，抄 step8 的 DEFAULT_RP_ALT
gaodu_xian = {2:300.0, 3:200.0, 4:200.0, 5:5000.0, 6:5000.0}

#发射日期上下界（儒略日）。写 None 就用 fanwei 里那默认的500天窗口。
riqi_qj = None

#每条腿飞行时间的上下界（天），想照着 step8 的算例写死就填，None 用默认。
tof_qj = None

zucheng = []      #最后一次算的明细，给打印和导出用
riqi_biao = []    #各节点日期
fenzhi_lie = []   #每条腿最后选中的是哪个分支
ke_xing = True    #这一次算出来的弹道可行不可行（缺陷<1e-3 且 rp 都够高）


#把一个矢量绕一根轴转一个角度，罗德里格公式，就三行，不用引别人的包
def xuanzhuan(zhou, jia, shi):
    z = zhou / np.linalg.norm(zhou)
    return (shi * math.cos(jia)
            + np.cross(z, shi) * math.sin(jia)
            + z * np.dot(z, shi) * (1 - math.cos(jia)))


#同一条腿的所有分支：绕几圈 x 顺行逆行 x 高路低路，全解出来列成一排
#（main2.py 是把圈数塞进染色体再取整，那是在糊弄差分进化；圈数是离散的，直接穷举就行）
def izzo_fenzhi(mu, r1, r2, tof_sec, M_ding):
    pai = []
    for M in range(0, M_ding + 1):
        if M == 0:
            fen = [True]                    #不绕圈只有一个解，高低路是一回事
        else:
            fen = [True, False]             #绕圈的有"高路/低路"两个分支
        for shun in [True, False]:
            for low in fen:
                try:
                    v1, v2 = izzo2015(mu, r1, r2, tof_sec, M=M, prograde=shun, low_path=low)
                except Exception:           #这个分支解不出来，跳过
                    continue
                if math.isnan(v1[0]) or math.isnan(v2[0]):
                    continue
                biao = 'M%d%s%s' % (M, '顺' if shun else '逆', '低' if low else '高')
                pai.append((v1, v2, biao))
    return pai


#一次零脉冲飞越：把 vin 尽量转到 vout 那个方向，还差多少就是缺陷
def feiyue_qing(vin, vout, mu_p, rp_min):
    chu = dict()
    mo_in = float(np.linalg.norm(vin))
    mo_out = float(np.linalg.norm(vout))
    if mo_in < 1e-9 or mo_out < 1e-9:
        return None
    a = vin / mo_in
    b = vout / mo_out
    c = float(np.dot(a, b))
    if c > 1.0:
        c = 1.0
    if c < -1.0:
        c = -1.0
    jia = math.acos(c)                                    #要转的角 δ
    v_inf = 0.5 * (mo_in + mo_out)                        #进和出的超速差不多，取平均
    jia_max = 2 * math.asin(min(1.0, 1.0 / (1.0 + rp_min * v_inf * v_inf / mu_p)))
    zhou = np.cross(a, b)
    if jia < 1e-12 or np.linalg.norm(zhou) < 1e-12:
        zhuanguo = vin                                    #本来就在一条线上，没得转
    else:
        zhuanguo = xuanzhuan(zhou, min(jia, jia_max), vin)
    que = float(np.linalg.norm(vout - zhuanguo))          #缺陷，单位 km/s
    sin_ban = math.sin(jia / 2)
    if sin_ban > 1e-9:
        rp_yao = mu_p / (v_inf * v_inf) * (1.0 / sin_ban - 1.0)   #这个δ唯一需要的近心半径
    else:
        rp_yao = 1e12
    #注意 rp_yao 必须用**真实转角**反解，不能用 min(δ,δmax) 截断，
    #否则转角不够的时候看不出违反约束（step8 注释里专门提了这一点）
    chu['mo_in'] = mo_in
    chu['mo_out'] = mo_out
    chu['v_inf'] = v_inf
    chu['jia'] = jia
    chu['jia_max'] = jia_max
    chu['que'] = que
    chu['rp'] = rp_yao
    return chu


def objective_function(agent, xulie):

    global zucheng, riqi_biao, fenzhi_lie, ke_xing

    L = len(xulie) - 1                      #几条腿
    riqi = agent[0]
    tian = []
    for k in range(L):
        tian.append(agent[1 + k])
    for k in range(L + 1):
        if k == 0:
            riqi_biao = [riqi]
        else:
            riqi_biao.append(riqi_biao[k - 1] + tian[k - 1])

    mu_sun, ra_sun = muandra(0)

    #--------- 先把每条腿在每个节点历元上的星球位置和所有分支都解出来 ---------
    r_lie = []
    v_lie = []
    for k in range(L + 1):
        r, v = randv(xulie[k], riqi_biao[k])
        r_lie.append(r)
        v_lie.append(v)

    fen_lie = []
    for k in range(L):
        pai = izzo_fenzhi(mu_sun, r_lie[k], r_lie[k + 1], tian[k] * 86400, M_zui_da)
        if len(pai) == 0:
            return 99999                    #这条腿一个分支都解不出来，这个日期不要
        fen_lie.append(pai)

    #--------- 逐条腿往下接：每条腿的分支全试，只留着最便宜的 liu_shao 种 ---------
    #状态就是一个三元组 (累计罚, 编号, 过程)，过程 = [各段的(v出发,v到达), 分支名...]
    qidong = []
    hao = 0
    for j in range(len(fen_lie[0])):
        v1, v2, biao = fen_lie[0][j]
        fa = float(np.linalg.norm(v1 - v_lie[0]))         #出发超速，最后再乘权重
        hao = hao + 1
        qidong.append((0.0, hao, [(v1, v2)], [biao], fa))
    for k in range(1, L):
        mu_p, ra_p = muandra(xulie[k])
        if xulie[k] in gaodu_xian:
            rp_min = ra_p + gaodu_xian[xulie[k]]
        else:
            rp_min = ra_p + 200.0
        xinshi = []
        for zhuang in qidong:
            v_qian = zhuang[2][-1][1]                     #上一段的到达速度
            vin_in = v_qian - v_lie[k]
            for j in range(len(fen_lie[k])):
                v1, v2, biao = fen_lie[k][j]
                vin_out = v1 - v_lie[k]
                ft = feiyue_qing(vin_in, vin_out, mu_p, rp_min)
                if ft is None:
                    continue
                hao = hao + 1
                lei = zhuang[0] + W_feiyue * ft['que'] * ft['que']
                xinshi.append((lei, hao, zhuang[2] + [(v1, v2)], zhuang[3] + [biao], zhuang[4]))
        if len(xinshi) == 0:
            return 99999
        xinshi.sort()                                     #按累计罚从小到大
        qidong = xinshi[0:liu_shao]

    #--------- 束里剩的这几条，补齐出发/到达超速，挑总分最小的那个 ---------
    zui_J = 1e30
    zui = None
    for zhuang in qidong:
        legs = zhuang[2]
        da = float(np.linalg.norm(legs[-1][1] - v_lie[L]))
        J = quan_xian_bu * zhuang[4] + quan_xian_da * da + zhuang[0]
        if vinf_shang is not None and zhuang[4] > vinf_shang:
            J = J + 100.0 * (zhuang[4] - vinf_shang) ** 2      #火箭打不上去，罚（step8同款）
        if J < zui_J:
            zui_J = J
            zui = (zhuang, da, J)
    if zui is None:
        return 99999

    zhuang, v_inf_da, J = zui
    legs = zhuang[2]
    fenzhi_lie = zhuang[3]

    #--------- 把明细一条一条记下来，顺便看这个弹道可不可行 ---------
    zucheng = []
    zucheng.append(float(np.linalg.norm(legs[0][0] - v_lie[0])))      #出发超速
    zucheng.append(v_inf_da)                                           #到达超速
    da_lei = []
    que_zui_da = 0.0
    for k in range(1, L):
        mu_p, ra_p = muandra(xulie[k])
        if xulie[k] in gaodu_xian:
            rp_min = ra_p + gaodu_xian[xulie[k]]
        else:
            rp_min = ra_p + 200.0
        vin_in = legs[k - 1][1] - v_lie[k]
        vin_out = legs[k][0] - v_lie[k]
        ft = feiyue_qing(vin_in, vin_out, mu_p, rp_min)
        da_lei.append((ft['rp'], ft['rp'] - ra_p, math.degrees(ft['jia']),
                       math.degrees(ft['jia_max']), ft['que']))
        if ft['que'] > que_zui_da:
            que_zui_da = ft['que']
        if ft['rp'] < rp_min:
            que_zui_da = 1.0                                          #贴太低了，也不算可行
    zucheng.append(da_lei)
    ke_xing = que_zui_da < 1e-3

    #--------- 顺手算一下"真实要花多少dv"：发射损失 + 到达捕获，中间一次机动都没有 ---------
    dv_fa = fashedv(xulie[0], legs[0][0], v_lie[0], jiatuo)
    dv_da = buhuodv(xulie[L], legs[-1][1], v_lie[L], buhuo_gaodu, yuandian)
    zucheng.append(dv_fa)
    zucheng.append(dv_da)

    if math.isnan(J):
        return 99999
    return J


#把算不出来的情况都变成一个大数，交给机器自己挑
def chengben(agent, xulie):
    try:
        return objective_function(agent, xulie)
    except Exception:
        return 99999


#初始化种群，确保每个元素都是两位小数
#zhong_zi：step8 里叫 --x-seed，把已经算出来的染色体塞进初始种群第一个，
#这样"先不加火箭上限搜一遍，再塞回去加上限搜一遍"就有个像样的起点了。不用就写None。
def initialize_population(pop_size, ranges, zhong_zi=None):
    renqun = [[round(random.uniform(ranges[i][0], ranges[i][1]), 2) for i in range(len(ranges))] for _ in range(pop_size)]
    if zhong_zi is not None:
        for j in range(len(zhong_zi)):
            lo, hi = ranges[j]
            zhi = round(zhong_zi[j], 2)
            if zhi < lo:
                zhi = lo
            if zhi > hi:
                zhi = hi
            renqun[0][j] = zhi
    return renqun


#差分进化算法主程序（跟 main2.py 一模一样，只是变量少了）
def differential_evolution(pop_size, F, CR, generations, ranges, xulie, zhong_zi=None):
    population = initialize_population(pop_size, ranges, zhong_zi)
    dim = len(ranges)

    fitness = []
    for agent in population:
        fitness.append(chengben(agent, xulie))

    for generation in range(generations):

        for k in range(pop_size):

            idxs = [idx for idx in range(pop_size) if idx != k]
            a, b, c = random.sample(idxs, 3)
            must = random.randint(0, dim - 1)

            mutant = []
            for j in range(dim):
                mutation = round(population[a][j] + F * (population[b][j] - population[c][j]), 2)
                mutant.append(max(ranges[j][0], min(ranges[j][1], mutation)))

            trial = []
            for j in range(dim):
                if random.random() < CR or j == must:
                    trial.append(mutant[j])
                else:
                    trial.append(population[k][j])

            trial_fitness = chengben(trial, xulie)
            if trial_fitness < fitness[k]:
                population[k] = trial
                fitness[k] = trial_fitness

        best_k = 0
        for k in range(pop_size):
            if fitness[k] < fitness[best_k]:
                best_k = k
        if generations <= 20 or generation % 20 == 0:
            print(f"第 {generation} 代最佳适应度: {fitness[best_k]:.4f}")

    best_fitness = float('inf')
    best_agent = None
    for agent in population:
        fitness_now = chengben(agent, xulie)
        if fitness_now < best_fitness:
            best_fitness = fitness_now
            best_agent = agent

    return best_agent, best_fitness


#step8 里 DE 完之后还要用 scipy 的 SLSQP 精修，把缺陷压到机器精度。
#不能引新包，就自己写个土办法：每个变量±动一小步，谁让J变小就留着，
#一整轮都没改善就把步长砍一半，砍到很小就停。这叫模式搜索，只有十几行。
def jing_xiu(agent, xulie, bu_fa):
    qu_jian = fanwei(xulie)
    zui = chengben(agent, xulie)
    xian = []
    for k in range(len(bu_fa)):
        xian.append(bu_fa[k])
    for lun_shu in range(0, 2000):                 #最多试2000轮，免得在噪声里打转
        gai_le = False
        for k in range(len(agent)):
            for qian in [1, -1]:
                shi = []
                for j in range(len(agent)):
                    shi.append(agent[j])
                shi[k] = shi[k] + qian * xian[k]
                lo, hi = qu_jian[k]                     #别跑出界
                if shi[k] < lo:
                    shi[k] = lo
                if shi[k] > hi:
                    shi[k] = hi
                zhi = chengben(shi, xulie)
                #只认比现在好1e-8以上的，再小的差别是数值噪声，别理它
                if zhi < zui - 0.00000001:
                    zui = zhi
                    agent = shi
                    gai_le = True
        zui_xian = 0.0
        if gai_le == False:
            for k in range(len(xian)):
                xian[k] = xian[k] / 2.0
                if xian[k] > zui_xian:
                    zui_xian = xian[k]
            if zui_xian < 0.000001:                     #步长砍到1e-6天(约0.1秒)就收工
                break
    return agent, zui


#这个版本只有 L+1 个变量，区间就是日期和各段天数
def fanwei(xulie):
    L = len(xulie) - 1
    if riqi_qj is None:
        ranges = [(2458553.5, 2459553.5)]    #发射日期，默认前后各500天
    else:
        ranges = [riqi_qj]
    for k in range(1, L + 1):
        if tof_qj is not None and k - 1 < len(tof_qj):
            ranges.append(tof_qj[k - 1])
            continue
        mubiao = xulie[k]
        if mubiao == 2:
            ranges.append((100, 400))        #到金星
        elif mubiao == 3:
            ranges.append((20, 700))         #回地球
        elif mubiao == 4:
            ranges.append((100, 600))        #到火星
        elif mubiao == 5:
            ranges.append((400, 2600))       #到木星
        else:
            ranges.append((100, 600))
    return ranges


def riqi_zh(jd):
    d = datetime.date(1970, 1, 1) + datetime.timedelta(days=int(jd - 2440587.5))
    return str(d)


names = {2: 'venus', 3: 'earth', 4: 'mars', 5: 'jupiter', 6: 'saturn'}
zhwen = {2: '金星', 3: '地球', 4: '火星', 5: '木星', 6: '土星'}


def print_detail(agent, zong, xulie):
    objective_function(agent, xulie)
    L = len(xulie) - 1
    hang = ''
    for k in range(len(xulie)):
        if k == 0:
            hang = '%s %s' % (zhwen[xulie[k]], riqi_zh(riqi_biao[0]))
        else:
            hang = hang + ' -> %s %s' % (zhwen[xulie[k]], riqi_zh(riqi_biao[k]))
    print(hang)
    print('  分支: %s   一共 %.0f 天' % (','.join(fenzhi_lie), riqi_biao[-1] - riqi_biao[0]))
    print('  出发 v_inf %.3f   到达 v_inf %.3f   J %.4f   %s'
          % (zucheng[0], zucheng[1], zong, '可行(零脉冲)' if ke_xing else '不可行'))
    if vinf_shang is not None and zucheng[0] > vinf_shang:
        print('    注意：出发超速 %.3f 超过火箭上限 %.2f（C3 %.1f），这一发打不上去，J里已经罚过了'
              % (zucheng[0], vinf_shang, zucheng[0] * zucheng[0]))
    for k in range(len(zucheng[2])):
        rp, alt, jia, jia_max, que = zucheng[2][k]
        print('    飞越%d %-4s  rp %9.1f km (高 %8.1f km)  δ %6.2f° / 最大 %6.2f°  缺陷 %8.3f m/s'
              % (k + 1, zhwen[xulie[k + 1]], rp, alt, jia, jia_max, que * 1000))
    print('  发射损失 %.3f + 捕获 %.3f = 真实总dv %.3f km/s（中间没有机动）'
          % (zucheng[-2], zucheng[-1], zucheng[-2] + zucheng[-1]))


#写成 step6/step7/step8 通用的 nodes JSON（单位 km、km/s、MJD），这样能直接丢进GMAT链去验
def export_nodes(agent, xulie, path):
    objective_function(agent, xulie)
    L = len(xulie) - 1
    mu_sun, ra_sun = muandra(0)
    r_lie = []
    v_lie = []
    for k in range(L + 1):
        r, v = randv(xulie[k], riqi_biao[k])
        r_lie.append(r)
        v_lie.append(v)

    #把选中的分支重算一遍，拿到每段的速度矢量（跟objective_function里一模一样）
    fen_lie = []
    for k in range(L):
        fen_lie.append(izzo_fenzhi(mu_sun, r_lie[k], r_lie[k + 1], tian_of(agent, L)[k] * 86400, M_zui_da))
    legs = []
    for k in range(L):
        for (v1, v2, biao) in fen_lie[k]:
            if biao == fenzhi_lie[k]:
                legs.append((v1, v2))
                break

    nodes = dict()
    nodes['departure'] = dict(body=names[xulie[0]],
                              epoch_mjd=float(riqi_biao[0] - 2400000.5),
                              vinf_km_s=[float(q) for q in (legs[0][0] - v_lie[0])])
    fei = []
    for k in range(1, L):
        mu_p, ra_p = muandra(xulie[k])
        vin_in = legs[k - 1][1] - v_lie[k]
        vin_out = legs[k][0] - v_lie[k]
        rp_min = ra_p + gaodu_xian.get(xulie[k], 200.0)
        ft = feiyue_qing(vin_in, vin_out, mu_p, rp_min)
        fei.append(dict(body=names[xulie[k]],
                        epoch_mjd=float(riqi_biao[k] - 2400000.5),
                        vinf_in_km_s=[float(q) for q in vin_in],
                        vinf_out_km_s=[float(q) for q in vin_out],
                        rp_km=float(ft['rp']),
                        alt_km=float(ft['rp'] - ra_p),
                        turn_deg=float(math.degrees(ft['jia'])),
                        turn_max_deg=float(math.degrees(ft['jia_max'])),
                        defect_m_s=float(ft['que'] * 1000)))
    nodes['flybys'] = fei
    nodes['dsms'] = [None] * L
    nodes['arrival'] = dict(epoch_mjd=float(riqi_biao[L] - 2400000.5),
                            target=names[xulie[L]])
    nodes['meta'] = dict(source='MGA_Planner main3 (纯MGA, 零脉冲借力, 无DSM)',
                         t0_mjd2000=float(riqi_biao[0] - 2400000.5 - 51544.5),
                         tof_days=[float(t) for t in tian_of(agent, L)],
                         objective=float(chengben(agent, xulie)),
                         vinf_dep_km_s=float(zucheng[0]),
                         vinf_arr_km_s=float(zucheng[1]),
                         branch_ids=[str(q) for q in fenzhi_lie],
                         max_flyby_defect_m_s=float(max([f['defect_m_s'] for f in fei], default=0.0)),
                         dv_launch_km_s=float(zucheng[-2]),
                         dv_capture_km_s=float(zucheng[-1]),
                         feasible=bool(ke_xing),
                         x=[float(q) for q in agent])
    with open(path, 'w', encoding='utf-8') as fh:
        json.dump(nodes, fh, indent=1, ensure_ascii=False)
    return nodes, path


def tian_of(agent, L):
    tian = []
    for k in range(L):
        tian.append(agent[1 + k])
    return tian


#精修的第一批步子：日期半天半天试，各段天数两天两天试
def chu_shi_bu_fa(L):
    bu_fa = [0.5]
    for k in range(L):
        bu_fa.append(2.0)
    return bu_fa


#照着 step8 的 --scan：发射窗口一片一片切，每片单独跑一轮短DE+精修，一行一个结果写进csv
#（好窗口不是猜出来的，是扫出来的；扫的时候窗口要给宽一点，别只给20天）
def sao_chuang(xulie, kai, zhi, bu, kuan, pop, gen, lu_jing):
    global riqi_qj
    hang = []
    riqi = kai
    while riqi <= zhi:
        riqi_qj = (riqi - kuan / 2.0, riqi + kuan / 2.0)
        random.seed(1)
        a, v = differential_evolution(pop, F, CR, gen, fanwei(xulie), xulie)
        a, v = jing_xiu(a, xulie, chu_shi_bu_fa(len(xulie) - 1))
        chao = ''
        for k in range(len(zucheng[2])):
            chao = chao + ',%.0f' % zucheng[2][k][1]      #每次飞越的高度km
        hang.append('%.2f,%.6f,%.4f,%.4f,%d%s,%.0f'
                    % (riqi, v, zucheng[0], zucheng[1], 1 if ke_xing else 0,
                       chao, riqi_biao[-1] - riqi_biao[0]))
        print('  窗口 %s  J=%.4f  出发v_inf=%.3f  到达v_inf=%.3f  %s'
              % (riqi_zh(riqi), v, zucheng[0], zucheng[1],
                 '可行' if ke_xing else '不可行'))
        riqi = riqi + bu
    wen = open(lu_jing, 'w', encoding='utf-8')
    wen.write('riqi,J,vinf_bu,vinf_da,kexing,feiyue_gaodu_km,zong_tian\n')
    for q in hang:
        wen.write(q + '\n')
    wen.close()
    print('扫完 %d 个窗口，已写 %s' % (len(hang), lu_jing))


#================== 参数设置 ==================
xulie = [3, 2, 3, 3, 5]     #地球-金星-地球-地球-木星（Galileo那种VEEGA，4条腿3次借力）

#照着 step8 的 eveej 算例给窗口：1989-06-01 ~ 1989-12-31
riqi_qj = (2447699.5 - 140.0, 2447699.5 + 140.0)
tof_qj = [(60, 220), (180, 420), (600, 900), (800, 1500)]

pop_size = 60
F = 0.5
CR = 0.9
generations = 150
chong_qi = 2
sao = False             #True = 顺便把发射窗口扫一遍（很慢，但是好窗口就是这么找的）

agent_best, best_fitness = None, 99999999.0
for lun in range(chong_qi):
    random.seed(lun + 1)
    agent_lun, dv_lun = differential_evolution(pop_size, F, CR, generations,
                                               fanwei(xulie), xulie)
    #DE 找到的点缺陷还有几十上百 m/s，再用土办法精修一下，把缺陷压到1 m/s以下
    agent_lun, dv_lun = jing_xiu(agent_lun, xulie, chu_shi_bu_fa(len(xulie) - 1))
    print()
    print('----- 第%d轮  J = %.6f（精修后）' % (lun + 1, dv_lun))
    print_detail(agent_lun, dv_lun, xulie)
    if dv_lun < best_fitness:
        agent_best, best_fitness = agent_lun, dv_lun

print()
print('最佳agent:', agent_best)
print('J =', round(best_fitness, 4), '(MGA惯例，只算两个超速加缺陷罚，不是总dv)')
print_detail(agent_best, best_fitness, xulie)
nodes, lu_jing = export_nodes(agent_best, xulie, 'main3_nodes.json')
print('已写出 %s（可以直接喂给 step6 的 --nodes）' % lu_jing)

#想扫窗口就把这行的False改成True（会先跑上面那一轮，再一片一片扫）
#参数：起、止（儒略日）、步长天、窗宽天、种群、代数
if sao:
    print()
    print('----- 开始扫发射窗口 -----')
    sao_chuang(xulie, 2447360.5, 2449447.5, 90.0, 150.0, 60, 60, 'main3_sao.csv')
