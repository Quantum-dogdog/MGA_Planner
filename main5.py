# -*- coding: utf-8 -*-
"""
main4.py —— 纯 MGA 借力搜索：序列枚举 + 引力弹弓自动搜任意共振腿

跟 main3.py 是同一个老流程，没拆模块、没搞类、没搞嵌套，全部平铺在一个文件里，
变量名还是那套拼音。想改参数就改【第 0 块】，然后从下往上看，一个函数是一件事。

跟 main3.py 比，多了两样东西（都是从「行星发射窗口计算器」学来的）：

一、序列不用手写一行了
    xulie 可以写编号 [3,2,3,3,5]，也可以写紧凑序列码 'EEJ'、'VVEJ'、'EVVEJ'，
    水星出发就写 'HJN'。想看某起点到某目标一共多少条路线，把 lie_yi_xia 改成
    True 就打印穷举结果（≤5 段）。序列码和 TOF 区间都是自动生成的。

二、同体腿自动搜共振借力
    两点 Lambert 在转移角≈0° 时退化：绕一圈再撞回同一颗星球（E→E、V→V…）
    的腿，缺陷怎么都压不到 0，所以老 main3.py 的结果里没有共振借力。
    这一版只要发现序列里 names[i] == names[i+1]（**任意**天体，**自动**识别，
    不用你声明），就给这条腿挂一个「圈数基因」，并且在【共振轮】里改走打靶：
      * 第一段就是同体腿（EEJ 的 E→E）：用「自由出发返回」，
        τ = n·行星周期，航天器周期 τ/k 反推半长轴和速度大小；
      * 中间的同体腿（EVEEJ 的 E→E）：用 (δ, ψ, τ) 三个量做 LM 打靶，
        通用变量开普勒传播，把「末端位置 - 行星位置」的三维残差压到最小。
    搜索分两轮，都是老流程那套差分进化 + 精修，只是成本函数换了一个：
      第一轮「普通轮」：纯 Lambert，跟老 main3.py 一模一样；
      第二轮「共振轮」：同体腿打靶，并且把 (n,k) 共振族种子塞进初始种群。
    最后两轮结果一起打印，共振解会标出「命中误差 / 飞越缺陷」，不冒充可行解。

还有一个顺手的改动：星历内置了。老 main3.py 要同目录有 10MB 的 de405.bsp，
没有就 import 失败。这里内置了一份 JPL 低精度行星星历（纯 numpy 数值表），
有 de405.bsp 就照旧用 SPK（精度不变），没有就自动用内置的，拷走这一个 py 就能跑。

文件目录（从上往下）：
    第 0 块  参数：改这里
    第 1 块  常量表 + 内置星历：muandra / randv / xingqi / ban_zhou
    第 2 块  序列枚举：jie_ma / lie_xulie / tongti_tui / she_zhi / tui_fanwei / fanwei
    第 3 块  老几何：xuanzhuan / izzo_fenzhi / feiyue_qing
    第 4 块  共振打靶：chuanbo / fangxiang / datui / ziyou_fanhui
    第 5 块  普通轮目标函数：objective_function（老 main3.py 原样）
    第 6 块  共振轮目标函数：ping_gongzhen
    第 7 块  搜索：chengben / initialize_population / differential_evolution /
                  jing_xiu / chu_shi_bu_fa / sao_chuang
    第 8 块  打印与导出：riqi_zh / print_detail / print_gongzhen / export_nodes
    第 9 块  Lambert 求解器（老仓库 izzo.py 原样搬来，只用到 izzo2015，不用细看）
    第 10 块 开跑（改完参数直接 python main4.py）
"""

import random
import math
import json
import datetime
import sys
import time              # 第 9 块那个 Lambert 求解器内部要计时，别删

import numpy as np
from numpy import cross, pi
from numpy.linalg import norm


# =========================================================================
# 第 0 块  参数：想改就改这里
# =========================================================================

# 序列：写编号列表 [3,2,3,3,5]，或写紧凑序列码 'EEJ' / 'VVEJ' / 'EVVEJ' / 'HJN'
# 也可以在命令行直接给： python main4.py EEJ
xulie = [3, 2, 3, 3, 5]        # 地球-金星-地球-地球-木星（伽利略 VEEGA）
qidian = 3                      # 出发天体编号：1水 2金 3地 4火 5木 6土 7天 8海 10冥
mubiao = 5                      # 目标天体编号（只在穷举时用）

# 把这一行改成 True，就只打印 qidian → mubiao 的全部 ≤zui_da_tui 段路线，不跑搜索
lie_yi_xia = False
zui_da_tui = 5                  # 最多几段

# 共振开关：同体腿自动改打靶，不用声明是哪条腿
gongzhen = True
k_zui_da = 3                    # 同体腿最多试几圈（1:1 / 2:1 / 3:1 共振返回）

# 下面这些跟老 main3.py 一个意思，值也一样
yuandian = 20                   # 捕获后椭圆远点 = 星球半径的多少倍
jiatuo = 200                    # 停泊轨道高度 km
buhuo_gaodu = 200
W_feiyue = 2000.0               # 飞越缺陷罚权重
quan_xian_bu = 1.0              # 出发超速权重
quan_xian_da = 1.0              # 到达超速权重
M_zui_da = 1                    # 非同体腿最多绕太阳几圈（圈数不进染色体，全试）
liu_shao = 20                   # 束搜索：每条腿只留最便宜的这么多种往下接
vinf_shang = None               # 火箭最大出发超速 km/s，None = 不管
gaodu_xian = {1: 200.0, 2: 300.0, 3: 200.0, 4: 200.0, 5: 5000.0,
              6: 5000.0, 7: 5000.0, 8: 5000.0, 10: 3000.0}

riqi_qj = (2447699.5 - 140.0, 2447699.5 + 140.0)   # 发射日期窗口（儒略日）
# 想复现 1989 VEEGA 参考解就留着这四段；换序列就改成 None，让 tui_fanwei 自动给
tof_qj = [(60, 220), (180, 420), (600, 900), (800, 1500)]

pop_size = 40                   # 普通轮 DE 种群
F = 0.5
CR = 0.9
generations = 80                # 想复现 1989 参考解（J≈13.7）调到 150 以上
chong_qi = 1                    # 普通轮重跑几次（换随机种子）
sui_ji_zhong = 1

gong_pop = 32                   # 共振轮 DE 种群（打靶比 Lambert 贵，给小点）
gong_gen = 80
sao = False                     # True = 顺便扫发射窗口（很慢）


# =========================================================================
# 第 1 块  常量表 + 内置星历
# =========================================================================

# 编号表（老 get_mu_ra.py 那套，1 号水星和 10 号冥王星是这次补的，2~9 没动）：
#   0太阳 1水星 2金星 3地球 4火星 5木星 6土星 7天王星 8海王星 9月球 10冥王星
mu_lie = [1.32712440018e11, 2.2032e4, 3.24859e5, 3.986004418e5,
          4.28284e4, 1.26713e8, 3.79406e7, 5.79456e6, 6.83653e6,
          4.90280e3, 8.696e2]
ra_lie = [695700, 2440.0, 6051.8, 6378.14, 3396.19, 71492,
          60268, 25559, 24764, 1737.4, 1188.3]


def muandra(shuzi):
    """返回 (mu km³/s², 赤道半径 km)。跟老 get_mu_ra.py 一模一样。"""
    i = int(shuzi)
    if i < 0 or i >= len(mu_lie):
        raise ValueError('未知天体编号: %s' % shuzi)
    return mu_lie[i], ra_lie[i]


AU = 1.495978707e8              # km
DAY = 86400.0                   # s
MU_SUN = 1.32712440018e11       # km³/s²
DU = math.pi / 180.0

# JPL 低精度行星星历根数表（a/e/i/L/peri/node + 每儒略世纪的变率）。
# 拟合区间 1800-2050，参考系日心黄道 J2000，地球取地月质心根数。
yuansu = {
    'mercury': dict(el=[0.38709927, 0.20563593, 7.00497902, 252.25032350, 77.45779628, 48.33076593],
                    dot=[0.00000037, 0.00001906, -0.00594749, 149472.67411175, 0.16047689, -0.12534081]),
    'venus':   dict(el=[0.72333566, 0.00677672, 3.39467605, 181.97909950, 131.60246718, 76.67984255],
                    dot=[0.00000390, -0.00004107, -0.00078890, 58517.81538729, 0.00268329, -0.27769418]),
    'earth':   dict(el=[1.00000261, 0.01671123, -0.00001531, 100.46457166, 102.93768193, 0.0],
                    dot=[0.00000562, -0.00004392, -0.01294668, 35999.37244981, 0.32327364, 0.0]),
    'mars':    dict(el=[1.52371034, 0.09339410, 1.84969142, -4.55343205, -23.94362959, 49.55953891],
                    dot=[0.00001847, 0.00007882, -0.00813131, 19140.30268499, 0.44441088, -0.29257343]),
    'jupiter': dict(el=[5.20288700, 0.04838624, 1.30439695, 34.39644051, 14.72847983, 100.47390909],
                    dot=[-0.00011607, -0.00013253, -0.00183714, 3034.74612775, 0.21252668, 0.20469106]),
    'saturn':  dict(el=[9.53667594, 0.05386179, 2.48599187, 49.95424423, 92.59887831, 113.66242448],
                    dot=[-0.00125060, -0.00050991, 0.00193609, 1222.49362201, -0.41897216, -0.28867794]),
    'uranus':  dict(el=[19.18916464, 0.04725744, 0.77263783, 313.23810451, 170.95427630, 74.01692503],
                    dot=[-0.00196176, -0.00004397, -0.00242939, 428.48202785, 0.40805281, 0.04240589]),
    'neptune': dict(el=[30.06992276, 0.00859048, 1.77004347, -55.12002969, 44.96476227, 131.78422574],
                    dot=[0.00026291, 0.00005105, 0.00035372, 218.45945325, -0.32241464, -0.00508664]),
    # 冥王星用旧版 Standish 根数，仅供外行星窗口初筛
    'pluto':   dict(el=[39.48211675, 0.24882730, 17.14001206, 238.92903833, 224.06891629, 110.30393684],
                    dot=[-0.00031596, 0.00005170, 0.00004818, 145.20780515, -0.04062942, -0.01183482]),
}

ming_zi = {1: 'mercury', 2: 'venus', 3: 'earth', 4: 'mars', 5: 'jupiter',
           6: 'saturn', 7: 'uranus', 8: 'neptune', 10: 'pluto'}
zhwen = {1: '水星', 2: '金星', 3: '地球', 4: '火星', 5: '木星',
         6: '土星', 7: '天王星', 8: '海王星', 10: '冥王星'}
zhduan = {1: '水', 2: '金', 3: '地', 4: '火', 5: '木',
          6: '土', 7: '天', 8: '海', 10: '冥'}
names_en = {1: 'mercury', 2: 'venus', 3: 'earth', 4: 'mars', 5: 'jupiter',
            6: 'saturn', 7: 'uranus', 8: 'neptune', 10: 'pluto'}


def gui_yi_jiao(x):
    t = x % (2.0 * math.pi)
    return t + 2.0 * math.pi if t < 0 else t


def ping_jin_zhen(M, e):
    """平近点角 -> 真近点角，牛顿迭代解开普勒方程。"""
    E = M + e * math.sin(M)
    for _ in range(40):
        f = E - e * math.sin(E) - M
        fp = 1.0 - e * math.cos(E)
        if abs(fp) < 1e-15:
            break
        d = f / fp
        E -= d
        if abs(d) < 1e-13:
            break
    return 2.0 * math.atan2(math.sqrt(1.0 + e) * math.sin(E / 2.0),
                            math.sqrt(1.0 - e) * math.cos(E / 2.0))


def xuan_z(v, a):
    c, s = math.cos(a), math.sin(a)
    return np.array([c * v[0] - s * v[1], s * v[0] + c * v[1], v[2]])


def xuan_x(v, a):
    c, s = math.cos(a), math.sin(a)
    return np.array([v[0], c * v[1] - s * v[2], s * v[1] + c * v[2]])


def gen_shu_zhuang(a, e, inc, Om, w, nu, mu):
    """六根数 -> 位置速度矢量（轨道面里算好再转三次）。"""
    p = a * (1.0 - e * e)
    cnu, snu = math.cos(nu), math.sin(nu)
    r = p / (1.0 + e * cnu)
    rpf = np.array([r * cnu, r * snu, 0.0])
    vpf = np.array([-math.sqrt(mu / p) * snu,
                    math.sqrt(mu / p) * (e + cnu), 0.0])
    return (xuan_z(xuan_x(xuan_z(rpf, w), inc), Om),
            xuan_z(xuan_x(xuan_z(vpf, w), inc), Om))


def xingli(shuzi, riqi):
    """内置星历：某天体在儒略日 riqi 的日心位置、速度（km、km/s）。"""
    key = ming_zi.get(int(shuzi))
    if key is None:
        raise ValueError('不支持的天体编号: %s' % shuzi)
    d = yuansu[key]
    t = float(riqi) - 2451545.0                 # JD -> MJD2000
    if not (-73048.0 < t < 44560.0):            # 1800 ~ 2122
        raise ValueError('内置星历有效期 1800-2122，当前 MJD2000=%.1f' % t)
    dt = t / 36525.0
    a = (d['el'][0] + d['dot'][0] * dt) * AU
    e = d['el'][1] + d['dot'][1] * dt
    inc = (d['el'][2] + d['dot'][2] * dt) * DU
    L = (d['el'][3] + d['dot'][3] * dt) * DU
    W = (d['el'][4] + d['dot'][4] * dt) * DU
    Om = (d['el'][5] + d['dot'][5] * dt) * DU
    w = W - Om
    M = gui_yi_jiao(L - W)
    nu = ping_jin_zhen(M, e)
    return gen_shu_zhuang(a, e, inc, Om, w, nu, MU_SUN)


# 有 de405.bsp 就照旧用 jplephem，精度不变；没有就自动回退到上面的内置星历。
try:
    from jplephem.spk import SPK
    _kernel = SPK.open('de405.bsp')
except Exception:
    _kernel = None


def randv(shuzi, riqi):
    """老 main3.py 里就是这个名字，单位 km、km/s，时间儒略日。"""
    if _kernel is not None:
        r, v = _kernel[0, 9 if int(shuzi) == 10 else int(shuzi)].compute_and_differentiate(riqi)
        return r, v / 86400.0                    # jplephem 给的是 km/day
    return xingli(shuzi, riqi)


def xingqi(r, v):
    """由日心状态算轨道周期（天），双曲轨道返回 None。共振族要用。"""
    eps = float(np.dot(v, v)) / 2.0 - MU_SUN / float(np.linalg.norm(r))
    if eps >= 0:
        return None
    a = -MU_SUN / (2.0 * eps)
    return 2.0 * math.pi * math.sqrt(a ** 3 / MU_SUN) / DAY


def ban_zhou(shuzi):
    """近似半长轴（AU），自动给 TOF 区间时用。"""
    return yuansu[ming_zi[int(shuzi)]]['el'][0]


# =========================================================================
# 第 2 块  序列枚举
# =========================================================================

dai_ma = {'H': 1, 'V': 2, 'E': 3, 'M': 4, 'J': 5,
          'S': 6, 'U': 7, 'N': 8, 'P': 10}
ma_biao = {v: k for k, v in dai_ma.items()}
shun_xu = [1, 2, 3, 4, 5, 6, 7, 8, 10]
pai_xu = {1: 0, 2: 1, 3: 2, 4: 3, 5: 4, 6: 5, 7: 6, 8: 7, 10: 8}

# 手工整理的常见腿 TOF 搜索盒（天）。只定搜索范围，不是解，刻意给宽。
tui_fw_biao = {
    (3, 3): (300, 1150), (3, 2): (60, 220), (2, 3): (40, 420),
    (3, 5): (500, 1800), (3, 6): (1600, 4000), (3, 7): (3000, 6000),
    (3, 8): (3500, 6500), (3, 10): (4000, 7000), (3, 1): (60, 320),
    (3, 4): (140, 460),
    (2, 5): (200, 900), (2, 6): (400, 1600), (2, 7): (600, 2400),
    (2, 8): (800, 3000), (2, 10): (800, 3000), (2, 1): (30, 220),
    (2, 4): (60, 320), (2, 2): (180, 900),
    (1, 3): (60, 420), (1, 2): (30, 220), (1, 4): (80, 420),
    (1, 5): (200, 1000), (1, 1): (80, 360),
    (4, 3): (140, 460), (4, 2): (60, 320), (4, 1): (80, 420),
    (4, 5): (280, 1700), (4, 4): (400, 2200),
    (5, 6): (700, 2600), (5, 7): (1400, 3800), (5, 8): (2800, 5200),
    (5, 10): (2200, 4600), (5, 4): (280, 1200), (5, 1): (200, 1000),
    (5, 5): (800, 5000),
    (6, 7): (1200, 4000), (6, 8): (2000, 5000), (6, 10): (1800, 4500),
    (6, 5): (700, 2600), (6, 6): (1500, 8000),
    (7, 6): (1200, 4000), (7, 7): (3000, 16000), (7, 8): (2000, 6000),
    (8, 7): (2000, 6000), (8, 8): (4000, 20000), (8, 10): (1000, 3500),
}

# 下面这几个全局量由 she_zhi() 一次算好，后面所有函数都直接拿来用
names = []          # 当前序列的天体编号列表
L = 0               # 几条腿
tui_fw_lie = []     # 每条腿的 TOF 搜索盒
k_tui = []          # 哪些腿是同体腿（自动识别，不用声明）
k_wei = {}          # 同体腿 -> 圈数基因在染色体里的下标
k_shang = {}        # 同体腿 -> 最多几圈
k_lie_chang = 0     # 染色体尾部有几个圈数基因


def jie_ma(xu, qi=3):
    """把 [3,2,3,3,5] 或 'EEJ' / 'VVEJ' / 'HJN' 变成编号列表。

    地球出发保留伽利略式简写：VVEJ 等价于 EVVEJ；其它起点必须写全。
    """
    if isinstance(xu, (list, tuple)):
        return [int(q) for q in xu]
    s = ''.join(c for c in str(xu).upper() if c.isalpha())
    if not s:
        return None
    if qi == 3 and s[0] != 'E':
        s = 'E' + s
    out = []
    for c in s:
        if c not in dai_ma:
            return None
        out.append(dai_ma[c])
    if len(out) < 2 or out[0] != qi or len(out) > 6:
        return None
    return out


def lie_xulie(mb, zui_da=5, qi=3):
    """穷举 qi -> mb 的全部 ≤zui_da 段路线，返回紧凑序列码列表。

    中间借力天体不超出「起点/目标中较外侧轨道再往内一格」，起点本身始终
    允许同体返回。地球出发、5 段时：木星 341 / 土星 781 / 天王星 1555 /
    海王星 2801 / 水星 121 条，跟参考项目一个不差。
    """
    mb, qi = int(mb), int(qi)
    out = []
    if mb not in pai_xu or mb == qi:
        return out
    cap = min(5, max(1, int(zui_da)))
    zui_wai = pai_xu[4] if mb == 1 else max(pai_xu[qi], pai_xu[mb]) - 1
    ke_xuan = []
    for b in shun_xu:
        if b == mb:
            continue
        if b != qi and pai_xu[b] > zui_wai:
            continue
        ke_xuan.append(b)

    def zou(ming_dan):
        if len(ming_dan) < cap - 1:
            for b in ke_xuan:
                zou(ming_dan + [b])
        if len(ming_dan) <= cap - 1:
            out.append(''.join(ma_biao[n] for n in [qi] + ming_dan + [mb]))

    zou([])
    return out


def tongti_tui(ming_dan):
    """所有同体腿的腿号。这就是共振轮要打靶的腿，自动识别，不硬编码。"""
    return [i for i in range(len(ming_dan) - 1)
            if ming_dan[i] == ming_dan[i + 1]]


def tui_fanwei(a, b):
    """一条腿的 TOF 搜索盒（天）。手工表没有就按半长轴估。"""
    v = tui_fw_biao.get((a, b))
    if v:
        return [float(v[0]), float(v[1])]
    ra, rb = ban_zhou(a), ban_zhou(b)
    if a == b:                                   # 同体腿要覆盖 1:1 / 2:1 / 3:1 共振
        P = 365.25 * ra ** 1.5
        return [max(40.0, round(0.6 * P)), round(3.2 * P)]
    half = 182.625 * ((ra + rb) * 0.5) ** 1.5
    return [max(30.0, round(0.55 * half)), round(2.5 * half)]


def she_zhi(ming_dan):
    """把当前序列设成全局量：TOF 区间 + 同体腿自动挂圈数基因。"""
    global names, L, tui_fw_lie, k_tui, k_wei, k_shang, k_lie_chang
    names = list(ming_dan)
    L = len(names) - 1
    tui_fw_lie = [tui_fanwei(names[i], names[i + 1]) for i in range(L)]
    k_tui = tongti_tui(names)                   # 关键：同体腿自动挂基因
    k_wei, k_shang = {}, {}
    for m, i in enumerate(k_tui):
        k_wei[i] = 1 + L + m                    # 染色体 = [t0, T1..TL, k1..km]
        k_shang[i] = k_zui_da
    k_lie_chang = len(k_tui)


def jie_k(agent):
    """把染色体尾部的圈数基因取整，返回每条腿的圈数（非同体腿为 0）。"""
    k_lie = [0] * L
    for i in k_tui:
        k_lie[i] = min(k_shang[i], max(1, int(round(agent[k_wei[i]]))))
    return k_lie


def fanwei(xu, yao_k=False):
    """染色体每一维的上下界。普通轮 L+1 维（跟老版一样），共振轮多几个圈数基因。"""
    n = len(xu) - 1
    if riqi_qj is None:
        ranges = [(2458553.5, 2459553.5)]
    else:
        ranges = [riqi_qj]
    for k in range(n):
        # tof_qj 是给默认 VEEGA 那四条腿写死的；换序列腿数对不上就自动忽略，
        # 免得 EEJ 被 (60,220) 这种短区间卡死（共振腿要一年以上）
        if tof_qj is not None and len(tof_qj) == n and k < len(tof_qj):
            ranges.append(tof_qj[k])
        else:
            ranges.append(tui_fw_lie[k])
    if yao_k:
        for i in k_tui:
            ranges.append((0.0, float(k_shang[i])))
    return ranges


def shuo_xulie(ming_dan):
    return ''.join(ma_biao.get(n, '?') for n in ming_dan) + ' · ' + \
           '-'.join(zhduan.get(n, str(n)) for n in ming_dan)


# =========================================================================
# 第 3 块  老几何（这几个函数跟 main3.py 一模一样）
# =========================================================================

def xuanzhuan(zhou, jia, shi):
    """把矢量绕一根轴转一个角，罗德里格公式，就三行。"""
    z = zhou / np.linalg.norm(zhou)
    return (shi * math.cos(jia)
            + np.cross(z, shi) * math.sin(jia)
            + z * np.dot(z, shi) * (1 - math.cos(jia)))


def izzo_fenzhi(mu, r1, r2, tof_sec, M_ding):
    """同一条腿的所有分支：绕几圈 × 顺行逆行 × 高路低路，全解出来列成一排。"""
    pai = []
    for M in range(0, M_ding + 1):
        fen = [True] if M == 0 else [True, False]
        for shun in [True, False]:
            for low in fen:
                try:
                    v1, v2 = izzo2015(mu, r1, r2, tof_sec, M=M,
                                      prograde=shun, low_path=low)
                except Exception:
                    continue
                if math.isnan(v1[0]) or math.isnan(v2[0]):
                    continue
                biao = 'M%d%s%s' % (M, '顺' if shun else '逆', '低' if low else '高')
                pai.append((v1, v2, biao))
    return pai


def feiyue_qing(vin, vout, mu_p, rp_min):
    """一次零脉冲飞越：把 vin 尽量转到 vout 那个方向，还差多少就是缺陷。"""
    chu = dict()
    mo_in = float(np.linalg.norm(vin))
    mo_out = float(np.linalg.norm(vout))
    if mo_in < 1e-9 or mo_out < 1e-9:
        return None
    a = vin / mo_in
    b = vout / mo_out
    c = float(np.dot(a, b))
    c = max(-1.0, min(1.0, c))
    jia = math.acos(c)                                    # 要转的角 δ
    v_inf = 0.5 * (mo_in + mo_out)
    jia_max = 2 * math.asin(min(1.0, 1.0 / (1.0 + rp_min * v_inf * v_inf / mu_p)))
    zhou = np.cross(a, b)
    if jia < 1e-12 or np.linalg.norm(zhou) < 1e-12:
        zhuanguo = vin
    else:
        zhuanguo = xuanzhuan(zhou, min(jia, jia_max), vin)
    que = float(np.linalg.norm(vout - zhuanguo))          # 缺陷 km/s
    sin_ban = math.sin(jia / 2)
    if sin_ban > 1e-9:
        rp_yao = mu_p / (v_inf * v_inf) * (1.0 / sin_ban - 1.0)
    else:
        rp_yao = 1e12
    # rp_yao 必须用真实转角反解，不能用 min(δ,δmax) 截断，否则看不出贴太低
    chu['mo_in'] = mo_in
    chu['mo_out'] = mo_out
    chu['v_inf'] = v_inf
    chu['jia'] = jia
    chu['jia_max'] = jia_max
    chu['que'] = que
    chu['rp'] = rp_yao
    return chu


# =========================================================================
# 第 4 块  共振打靶（这一段是这次新加的，全是平铺的小函数）
# =========================================================================

def stumpff_c(z):
    if z > 1e-6:
        s = math.sqrt(z)
        return (1.0 - math.cos(s)) / z
    if z < -1e-6:
        s = math.sqrt(-z)
        return (math.cosh(s) - 1.0) / (-z)
    return 0.5 - z / 24.0 + z * z / 720.0 - z ** 3 / 40320.0


def stumpff_s(z):
    if z > 1e-6:
        s = math.sqrt(z)
        return (s - math.sin(s)) / (z * s)
    if z < -1e-6:
        s = math.sqrt(-z)
        return (math.sinh(s) - s) / (s * -z)
    return 1.0 / 6.0 - z / 120.0 + z * z / 5040.0 - z ** 3 / 362880.0


def chuanbo(r0, v0, dt):
    """通用变量开普勒传播：从 (r0,v0) 走 dt 秒，返回 (r,v)。"""
    r0 = np.asarray(r0, float)
    v0 = np.asarray(v0, float)
    rn = float(np.linalg.norm(r0))
    vn = float(np.linalg.norm(v0))
    rv = float(np.dot(r0, v0))
    alpha = 2.0 / rn - vn * vn / MU_SUN
    if alpha > 1e-10:
        x = math.sqrt(MU_SUN) * dt * alpha
    elif alpha < -1e-10:
        sg = 1.0 if dt >= 0 else -1.0
        arg = (-2.0 * MU_SUN * alpha * dt) / (
            rv + sg * math.sqrt(-MU_SUN / alpha) * (1.0 - rn * alpha))
        if arg <= 0:
            return None
        x = sg * math.sqrt(-1.0 / alpha) * math.log(arg)
    else:
        x = math.sqrt(MU_SUN) * dt / rn
    for _ in range(80):
        z = alpha * x * x
        C, S = stumpff_c(z), stumpff_s(z)
        F = ((rv / math.sqrt(MU_SUN)) * x * x * C
             + (1.0 - alpha * rn) * x ** 3 * S + rn * x - math.sqrt(MU_SUN) * dt)
        dF = ((rv / math.sqrt(MU_SUN)) * x * (1.0 - z * S)
              + (1.0 - alpha * rn) * x * x * C + rn)
        if not math.isfinite(dF) or abs(dF) < 1e-14:
            break
        dx = F / dF
        x -= dx
        if abs(dx) < 1e-10 * max(1.0, abs(x)):
            break
    z = alpha * x * x
    C, S = stumpff_c(z), stumpff_s(z)
    f = 1.0 - (x * x / rn) * C
    g = dt - (x ** 3 / math.sqrt(MU_SUN)) * S
    r = f * r0 + g * v0
    rnn = float(np.linalg.norm(r))
    if rnn < 1e-9:
        return None
    fd = (math.sqrt(MU_SUN) / (rnn * rn)) * x * (z * S - 1.0)
    gd = 1.0 - (x * x / rnn) * C
    return r, fd * r0 + gd * v0


def fangxiang(a_hat, delta, psi):
    """和入射渐近线方向 a_hat 成 delta 度、绕它转 psi 度的单位方向（出射渐近线）。"""
    a = np.asarray(a_hat, float)
    a = a / np.linalg.norm(a)
    tmp = np.array([0.0, 0.0, 1.0]) if abs(a[2]) < 0.9 else np.array([1.0, 0.0, 0.0])
    p1 = np.cross(a, tmp)
    p1 = p1 / np.linalg.norm(p1)
    p2 = np.cross(a, p1)
    d = math.radians(delta)
    q = math.radians(psi)
    return a * math.cos(d) + (p1 * math.cos(q) + p2 * math.sin(q)) * math.sin(d)


nuan = {}       # 热启动：上一条染色体在同一 (腿, k) 上打中的解，几次迭代就能收敛


def datui(tui_hao, r0, v0, vin, ep0, T_seed, k, mubiao, mu_p, rp_min):
    """一条中段同体返回腿的 (δ, ψ, τ) 打靶。

    航天器从 v_planet + |v∞|·u(δ,ψ) 出发，用通用变量传播 τ 天，
    Levenberg-Marquardt 把「末端位置 - 行星位置」的三维残差压到最小。
    τ 是输出量不是约束；靠 |k_eff - k| < 0.35 挑共振族（k 是染色体给的圈数）。
    """
    nv = float(np.linalg.norm(vin))
    if not (0.05 < nv < 100.0):
        return None
    sd = 1.0 / (1.0 + rp_min * nv * nv / mu_p)
    if not (0.0 < sd < 1.0):
        return None
    dmax = math.degrees(2.0 * math.asin(sd))        # 这颗星球最多能转的角度
    if dmax <= 1e-3:
        return None
    a_hat = vin / nv
    lo, hi = T_seed * 0.8, T_seed * 1.2
    dlo, dhi = 0.5, min(170.0, dmax)
    if hi <= lo or dhi <= dlo:
        return None

    def zhi(q):
        """归一化变量 -> (δ度, ψ度, τ天)"""
        return [dlo + q[0] * (dhi - dlo),
                -720.0 + q[1] * 1440.0,
                lo + q[2] * (hi - lo)]

    def gui(z):
        return [max(0.0, min(1.0, (z[0] - dlo) / (dhi - dlo))),
                max(0.0, min(1.0, (z[1] + 720.0) / 1440.0)),
                max(0.0, min(1.0, (z[2] - lo) / (hi - lo)))]

    def can(z):
        """给定三个量，算出末端位置残差"""
        u = fangxiang(a_hat, z[0], z[1])
        vsc = v0 + u * nv
        pr = chuanbo(r0, vsc, z[2] * DAY)
        if pr is None:
            return None
        rb, _ = randv(mubiao, ep0 + z[2])
        return dict(res=pr[0] - rb, vsc=vsc, varr=pr[1], z=z, u=u)

    def jie(g):
        if g is None:
            return None
        D = float(np.linalg.norm(g['res']))
        eps = float(np.dot(g['vsc'], g['vsc'])) / 2.0 - MU_SUN / float(np.linalg.norm(r0))
        a_km = -MU_SUN / (2.0 * eps) if eps < 0 else float('nan')
        P = 2.0 * math.pi * math.sqrt(a_km ** 3 / MU_SUN) / DAY if a_km > 0 else float('nan')
        k_eff = g['z'][2] / P if P > 0 else 1e9
        return dict(D=D, v1=g['vsc'], v2=g['varr'], tof=g['z'][2], k_eff=k_eff)

    def lm(q0, iters):
        q = np.array(q0, float)
        g = can(zhi(q))
        if g is None:
            return None
        nr = float(np.linalg.norm(g['res']))
        lam, stall = 1.0, 0
        for _ in range(iters):
            A = np.zeros((3, 3))
            h = 1e-6
            for j in range(3):
                qp = q.copy()
                qp[j] = min(1.0, qp[j] + h) if qp[j] + h <= 1.0 else 1.0 - h
                gm = can(zhi(qp))
                if gm is not None:
                    A[:, j] = (gm['res'] - g['res']) / h
            try:
                dq = np.linalg.solve(A.T @ A + lam * np.eye(3), -A.T @ g['res'])
            except np.linalg.LinAlgError:
                break
            mx = float(np.max(np.abs(dq)))
            step = 0.5 / mx if mx > 0.5 else 1.0
            qn = np.clip(q + dq * step, 0.0, 1.0)
            gn = can(zhi(qn))
            if gn is None:
                lam *= 6.0
                stall += 1
                if lam > 1e8 or stall > 10:
                    break
                continue
            nn = float(np.linalg.norm(gn['res']))
            if nn < nr:
                q, g, nr = qn, gn, nn
                lam, stall = max(1e-9, lam / 4.0), 0
            else:
                lam *= 6.0
                stall += 1
                if lam > 1e8 or stall > 10:
                    break
            if nn < 1e-4:
                break
        return jie(g)

    # 热启动：同一个 (腿, k) 上一次的解通常就在新窗口旁边
    wen = (tui_hao, k)
    if wen in nuan:
        w = nuan[wen]
        if lo <= w['tof'] <= hi:
            fr = lm(gui([w['delta'], w['psi'], w['tof']]), 20)
            if fr and fr['D'] < 2e3 and abs(fr['k_eff'] - k) < 0.35:
                return fr

    # 冷启动：先用下一段 Lambert 的出发渐近线给几个像样的初值，再补一圈网格
    qishi = []
    rb0, _ = randv(mubiao, ep0 + T_seed)
    for (v1l, _, _) in izzo_fenzhi(MU_SUN, r0, rb0, T_seed * DAY, k):
        du = v1l - v0
        dn = float(np.linalg.norm(du))
        if dn <= 1e-3 or abs(dn - nv) > nv * 0.05:
            continue
        u = du / dn
        dd = math.degrees(math.acos(max(-1.0, min(1.0, float(np.dot(a_hat, u))))))
        tmp = np.array([0.0, 0.0, 1.0]) if abs(a_hat[2]) < 0.9 else np.array([1.0, 0.0, 0.0])
        p1 = np.cross(a_hat, tmp)
        p1 = p1 / np.linalg.norm(p1)
        p2 = np.cross(a_hat, p1)
        psi = math.degrees(math.atan2(float(np.dot(u, p2)), float(np.dot(u, p1))))
        if dlo <= dd <= dhi:
            qishi.append([dd, psi, T_seed])
    for fd in (0.2, 0.45, 0.7, 0.95):
        for psi in (0.0, 90.0, 180.0, 270.0):
            qishi.append([min(fd * dmax, 170.0), psi, T_seed])

    zui = None
    for z0 in qishi:
        fr = lm(gui(z0), 40)
        if fr is None or fr['D'] > 1e5 or abs(fr['k_eff'] - k) > 0.35:
            continue
        pai = (fr['D'], abs(fr['k_eff'] - k), z0)
        if zui is None or pai[:2] < zui[0]:
            zui = (pai[:2], fr, z0)
    if zui is None:
        return None
    nuan[wen] = dict(delta=zui[2][0], psi=zui[2][1], tof=zui[1]['tof'])
    return zui[1]


def ziyou_fanhui(t0, tian, k_lie, r0, v0, xu, rp_min):
    """首段自由出发返回（EEJ 这种第一段就是同体腿）。

    没有入射渐近线可打靶，就令 τ = n·行星周期，航天器周期 τ/k 反推半长轴 a，
    再用 vis-viva 定出速度大小 |vdep|；方向从下一段 Lambert 的出射渐近线里挑。
    返回 (vdep, τ, 下一段分支号)。分支号要带下去，否则两条腿各说各话。
    """
    P0 = xingqi(r0, v0)
    if not P0 or P0 <= 0:
        return None
    n0 = min(max(int(round(tian[0] / P0)), 1), k_zui_da)
    kk = max(1, int(k_lie[0]))
    tau = n0 * P0
    ep1 = t0 + tau
    r1, v1 = randv(xu[1], ep1)
    rr = float(np.linalg.norm(r1))
    a = (MU_SUN * (tau * DAY / (2.0 * math.pi * kk)) ** 2) ** (1.0 / 3.0)
    vmag = math.sqrt(max(0.0, MU_SUN * (2.0 / rr - 1.0 / a)))
    if not (vmag > 0 and math.isfinite(vmag)):
        return None
    mu_p, _ = muandra(xu[1])

    sols = []
    if L >= 2:
        ep2 = ep1 + tian[1]
        r2, _ = randv(xu[2], ep2)
        sols = izzo_fenzhi(MU_SUN, r1, r2, tian[1] * DAY, M_zui_da + k_lie[1])
    if not sols:
        nv1 = float(np.linalg.norm(v1))
        u = v1 / nv1 if nv1 > 0 else np.array([1.0, 0.0, 0.0])
        return (u * vmag, tau, None)

    u0 = v0 / float(np.linalg.norm(v0))
    v1n2 = float(np.dot(v1, v1))
    zui = None
    for c, (v1l, _, _) in enumerate(sols):
        vout = v1l - v1
        R = float(np.linalg.norm(vout))
        na = float(np.linalg.norm(v1l))
        if R < 1e-9:
            continue
        u = vout / R
        # 三个自然方向 + 一个「让 vin 和 vout 共线」的解析方向
        cands = []
        if na > 0:
            cands.append(v1l / na * vmag)
        cands.append(u * vmag)
        cands.append(u0 * vmag)
        b = float(np.dot(v1, u))
        disc = b * b + vmag * vmag - v1n2
        if disc >= 0:
            s = -b + math.sqrt(disc)
            if s > 0:
                cands.append(v1 + s * u)
        for vdep in cands:
            ft = feiyue_qing(vdep - v1, vout, mu_p, rp_min)
            if ft is None:
                continue
            sc = float(np.linalg.norm(vdep - v0)) + 100.0 * ft['que'] ** 2
            if zui is None or sc < zui[0]:
                zui = (sc, vdep, c)
    if zui is None:
        return None
    return (zui[1], tau, zui[2])


# =========================================================================
# 第 5 块  普通轮目标函数（老 main3.py 原样，纯 Lambert）
# =========================================================================

zucheng = []
riqi_biao = []
fenzhi_lie = []
ke_xing = True


def objective_function(agent, xu):
    global zucheng, riqi_biao, fenzhi_lie, ke_xing

    n = len(xu) - 1
    riqi = agent[0]
    tian = [agent[1 + k] for k in range(n)]
    riqi_biao = [riqi]
    for k in range(n):
        riqi_biao.append(riqi_biao[k] + tian[k])

    mu_sun, ra_sun = muandra(0)

    # 先把每条腿在每个节点历元上的星球位置和所有分支都解出来
    r_lie, v_lie = [], []
    for k in range(n + 1):
        r, v = randv(xu[k], riqi_biao[k])
        r_lie.append(r)
        v_lie.append(v)

    fen_lie = []
    for k in range(n):
        pai = izzo_fenzhi(mu_sun, r_lie[k], r_lie[k + 1],
                          tian[k] * 86400, M_zui_da)
        if len(pai) == 0:
            return 99999
        fen_lie.append(pai)

    # 逐条腿往下接：每条腿的分支全试，只留着最便宜的 liu_shao 种
    # 状态 = (累计罚, 编号, 过程[(v出发,v到达)...], 分支名..., 出发超速)
    qidong = []
    hao = 0
    for j in range(len(fen_lie[0])):
        v1, v2, biao = fen_lie[0][j]
        fa = float(np.linalg.norm(v1 - v_lie[0]))
        hao += 1
        qidong.append((0.0, hao, [(v1, v2)], [biao], fa))
    for k in range(1, n):
        mu_p, ra_p = muandra(xu[k])
        rp_min = ra_p + gaodu_xian.get(xu[k], 200.0)
        xinshi = []
        for zhuang in qidong:
            v_qian = zhuang[2][-1][1]
            vin_in = v_qian - v_lie[k]
            for j in range(len(fen_lie[k])):
                v1, v2, biao = fen_lie[k][j]
                vin_out = v1 - v_lie[k]
                ft = feiyue_qing(vin_in, vin_out, mu_p, rp_min)
                if ft is None:
                    continue
                hao += 1
                lei = zhuang[0] + W_feiyue * ft['que'] * ft['que']
                xinshi.append((lei, hao, zhuang[2] + [(v1, v2)],
                               zhuang[3] + [biao], zhuang[4]))
        if len(xinshi) == 0:
            return 99999
        xinshi.sort()
        qidong = xinshi[0:liu_shao]

    # 束里剩的这几条，补齐出发/到达超速，挑总分最小的那个
    zui_J = 1e30
    zui = None
    for zhuang in qidong:
        legs = zhuang[2]
        da = float(np.linalg.norm(legs[-1][1] - v_lie[n]))
        J = quan_xian_bu * zhuang[4] + quan_xian_da * da + zhuang[0]
        if vinf_shang is not None and zhuang[4] > vinf_shang:
            J += 100.0 * (zhuang[4] - vinf_shang) ** 2
        if J < zui_J:
            zui_J = J
            zui = (zhuang, da, J)
    if zui is None:
        return 99999

    zhuang, v_inf_da, J = zui
    legs = zhuang[2]
    fenzhi_lie = zhuang[3]

    # 把明细一条一条记下来，顺便看这个弹道可不可行
    zucheng = []
    zucheng.append(float(np.linalg.norm(legs[0][0] - v_lie[0])))
    zucheng.append(v_inf_da)
    da_lei = []
    que_zui_da = 0.0
    for k in range(1, n):
        mu_p, ra_p = muandra(xu[k])
        rp_min = ra_p + gaodu_xian.get(xu[k], 200.0)
        vin_in = legs[k - 1][1] - v_lie[k]
        vin_out = legs[k][0] - v_lie[k]
        ft = feiyue_qing(vin_in, vin_out, mu_p, rp_min)
        da_lei.append((ft['rp'], ft['rp'] - ra_p, math.degrees(ft['jia']),
                       math.degrees(ft['jia_max']), ft['que']))
        if ft['que'] > que_zui_da:
            que_zui_da = ft['que']
        if ft['rp'] < rp_min:
            que_zui_da = 1.0
    zucheng.append(da_lei)
    ke_xing = que_zui_da < 1e-3

    # 顺手算一下真实要花多少 dv：发射损失 + 到达捕获，中间一次机动都没有
    dv_fa = fashedv(xu[0], legs[0][0], v_lie[0], jiatuo)
    dv_da = buhuodv(xu[n], legs[-1][1], v_lie[n], buhuo_gaodu, yuandian)
    zucheng.append(dv_fa)
    zucheng.append(dv_da)

    if math.isnan(J):
        return 99999
    return J


# 老 get_dv.py 里就这两个函数用到，20 行，一起搬进来，省得再引一个文件
def fashedv(shuzi, v_need, v_now, gaodu):
    v_inf_out = np.linalg.norm(v_need - v_now)
    mu, ra = muandra(shuzi)
    r_park = ra + gaodu
    v_park = math.sqrt(mu / r_park)
    v_after = math.sqrt(v_inf_out * v_inf_out + 2 * mu / r_park)
    return v_after - v_park


def buhuodv(shuzi, v_arrive, v_now, gaodu, yuandian=0):
    mu, ra = muandra(shuzi)
    r_park = ra + gaodu
    v_inf_in = np.linalg.norm(v_arrive - v_now)
    a_arrive = -mu / (v_inf_in * v_inf_in) if v_inf_in > 0 else -1e12
    v_inpe = math.sqrt(2 * mu / r_park - mu / a_arrive)
    r_yuan = r_park if yuandian == 0 else ra * yuandian
    a_zhu = (r_park + r_yuan) / 2
    v_zhuan = math.sqrt(2 * mu / r_park - mu / a_zhu)
    return v_inpe - v_zhuan


# =========================================================================
# 第 6 块  共振轮目标函数
# =========================================================================

gong_mingxi = None


def ping_gongzhen(agent, xu):
    """共振轮：同体腿自动打靶，其它腿 Lambert，逐段往前推，返回 (J, 明细)。

    跟第 5 块的束搜索不同，这里每段只挑「让飞越缺陷最小」的那一个分支往下接
    （贪心），因为打靶腿依赖上一段的到达速度，没法先全解出来再组合。
    普通轮还是老流程全套束搜索，两条路互不干扰。
    """
    global gong_mingxi
    n = len(xu) - 1
    tian = [agent[1 + k] for k in range(n)]
    k_lie = jie_k(agent)

    ep = [agent[0]]
    r0, v0 = randv(xu[0], ep[0])
    r_lie, v_lie = [r0], [v0]
    legs, fenzhi, hit_lie = [], [], []
    fa = 0.0                        # 飞越缺陷平方和
    que_zui = 0.0
    fei_lie = []
    qiang = {}                      # 首段自由返回绑定的下一段分支号

    for i in range(n):
        ep0 = ep[-1]
        r_p, v_p = r_lie[-1], v_lie[-1]
        mu_p, ra_p = muandra(xu[i]) if i > 0 else muandra(xu[0])
        rp_min = ra_p + gaodu_xian.get(xu[i], 200.0) if i > 0 else \
            muandra(xu[1])[1] + gaodu_xian.get(xu[1], 200.0)

        if i in k_tui and gongzhen:
            if i == 0:
                # 第一段就是同体腿：自由出发返回
                jie = ziyou_fanhui(ep0, tian, k_lie, r_p, v_p, xu, rp_min)
                if jie is None:
                    return 99999, None
                vdep, tau, c = jie
                v1, v2 = vdep, vdep
                tof = tau
                fenzhi.append('自由返回')
                if c is not None:
                    qiang[1] = c
                r1, _ = randv(xu[1], ep0 + tau)
                hit_lie.append(float(np.linalg.norm(r1 - r_p)))
            else:
                # 中间同体腿：(δ,ψ,τ) 打靶
                vin = legs[-1][1] - v_p
                rec = datui(i, r_p, v_p, vin, ep0, tian[i], k_lie[i],
                            xu[i + 1], mu_p, rp_min)
                if rec is None:
                    return 99999, None
                v1, v2 = rec['v1'], rec['v2']
                tof = rec['tof']
                fenzhi.append('打靶')
                hit_lie.append(rec['D'])
                ft = feiyue_qing(vin, v1 - v_p, mu_p, rp_min)
                if ft is None:
                    return 99999, None
            legs.append((v1, v2))
            ep.append(ep0 + tof)
        else:
            ep.append(ep0 + tian[i])
            r2, _ = randv(xu[i + 1], ep[-1])
            pai = izzo_fenzhi(MU_SUN, r_p, r2, tian[i] * DAY,
                              M_zui_da + k_lie[i])
            if not pai:
                return 99999, None
            if i in qiang:
                pai = [pai[min(qiang[i], len(pai) - 1)]]
            # 挑让飞越缺陷最小的一支（首段再看出发超速，末段再看到达超速）
            zui_pai = None
            for (a1, a2, biao) in pai:
                score = 0.0
                if i > 0:
                    ft = feiyue_qing(legs[-1][1] - v_p, a1 - v_p, mu_p, rp_min)
                    if ft is None:
                        continue
                    score += W_feiyue * ft['que'] ** 2
                if i == 0:
                    score += float(np.linalg.norm(a1 - v_p))
                if i == n - 1:
                    score += float(np.linalg.norm(a2 - randv(xu[n], ep[-1])[1]))
                if zui_pai is None or score < zui_pai[0]:
                    zui_pai = (score, a1, a2, biao, ft if i > 0 else None)
            if zui_pai is None:
                return 99999, None
            v1, v2, biao = zui_pai[1], zui_pai[2], zui_pai[3]
            legs.append((v1, v2))
            fenzhi.append(biao)

        r_end, v_end = randv(xu[i + 1], ep[-1])
        r_lie.append(r_end)
        v_lie.append(v_end)

    # 每个飞越点再算一遍明细（打印用）
    for k in range(1, n):
        mu_p, ra_p = muandra(xu[k])
        rp_min = ra_p + gaodu_xian.get(xu[k], 200.0)
        ft = feiyue_qing(legs[k - 1][1] - v_lie[k], legs[k][0] - v_lie[k],
                         mu_p, rp_min)
        if ft is None:
            return 99999, None
        fei_lie.append((ft['rp'], ft['rp'] - ra_p, math.degrees(ft['jia']),
                        math.degrees(ft['jia_max']), ft['que']))
        # 罚项统一用「最终 legs 重算出来」的缺陷，跟上面打印的明细是同一个量。
        # 早先贪心选分支时用的是临时值，DE 压的和打印的对不上，会漏掉打靶腿
        # 之后那个节点（只保证位置命中、不保证速度接续）的大残差。
        fa += ft['que'] ** 2
        que_zui = max(que_zui, ft['que'])

    vd = float(np.linalg.norm(legs[0][0] - v_lie[0]))
    va = float(np.linalg.norm(legs[-1][1] - v_lie[-1]))
    max_hit = max(hit_lie) if hit_lie else 0.0
    hit_fa = 10.0 * sum((d / 1e4) ** 2 for d in hit_lie)
    J = quan_xian_bu * vd + quan_xian_da * va + W_feiyue * fa + hit_fa
    if vinf_shang is not None and vd > vinf_shang:
        J += 100.0 * (vd - vinf_shang) ** 2
    if not math.isfinite(J):
        return 99999, None

    dv_fa = fashedv(xu[0], legs[0][0], v_lie[0], jiatuo)
    dv_da = buhuodv(xu[n], legs[-1][1], v_lie[n], buhuo_gaodu, yuandian)
    gong_mingxi = dict(ep=ep[:], tian=[ep[q + 1] - ep[q] for q in range(n)],
                       vd=vd, va=va, max_defect=que_zui, max_hit=max_hit,
                       fenzhi=fenzhi[:], fei_lie=fei_lie, k_lie=k_lie,
                       dv_fa=dv_fa, dv_da=dv_da,
                       feasible=(que_zui < 1e-3 and max_hit < 1e4))
    return J, gong_mingxi


# =========================================================================
# 第 7 块  搜索（差分进化 + 精修，跟老 main3.py 一模一样，只多一个 moshi）
# =========================================================================

def chengben(agent, xu, moshi='putong'):
    """算不出来的情况都变成一个大数，交给机器自己挑。moshi='gongzhen' 走共振轮。"""
    try:
        if moshi == 'gongzhen':
            return ping_gongzhen(agent, xu)[0]
        return objective_function(agent, xu)
    except Exception:
        return 99999


def initialize_population(pop_size, ranges, zhong_zi=None):
    """初始化种群，两位数。zhong_zi 可以是一条染色体，也可以是好几条（共振种子）。"""
    renqun = [[round(random.uniform(ranges[i][0], ranges[i][1]), 2)
               for i in range(len(ranges))] for _ in range(pop_size)]
    if zhong_zi is not None:
        if len(zhong_zi) > 0 and isinstance(zhong_zi[0], (list, tuple)):
            duo = zhong_zi
        else:
            duo = [zhong_zi]
        for m, zz in enumerate(duo):
            if m >= pop_size:
                break
            for j in range(min(len(zz), len(ranges))):
                lo, hi = ranges[j]
                renqun[m][j] = round(min(hi, max(lo, zz[j])), 2)
    return renqun


def differential_evolution(pop_size, F, CR, generations, ranges, xu,
                           zhong_zi=None, moshi='putong'):
    population = initialize_population(pop_size, ranges, zhong_zi)
    dim = len(ranges)
    fitness = [chengben(a, xu, moshi) for a in population]

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
            trial_fitness = chengben(trial, xu, moshi)
            if trial_fitness < fitness[k]:
                population[k] = trial
                fitness[k] = trial_fitness
        best_k = min(range(pop_size), key=lambda q: fitness[q])
        if generations <= 20 or generation % 20 == 0:
            print('  %s 第 %d 代最佳 J: %.4f'
                  % ('共振轮' if moshi == 'gongzhen' else '普通轮',
                     generation, fitness[best_k]))

    best_fitness = float('inf')
    best_agent = None
    for agent in population:
        fitness_now = chengben(agent, xu, moshi)
        if fitness_now < best_fitness:
            best_fitness = fitness_now
            best_agent = agent
    return best_agent, best_fitness


def jing_xiu(agent, xu, bu_fa, moshi='putong'):
    """每个变量 ±动一小步，谁让 J 变小就留着；一整轮没改善就把步长砍一半。"""
    qu_jian = fanwei(xu, yao_k=(moshi == 'gongzhen'))
    zui = chengben(agent, xu, moshi)
    xian = list(bu_fa)
    for _ in range(0, 2000):
        gai_le = False
        for k in range(len(agent)):
            if k >= len(xian) or xian[k] <= 0:
                continue
            for qian in [1, -1]:
                shi = list(agent)
                shi[k] = shi[k] + qian * xian[k]
                lo, hi = qu_jian[k]
                shi[k] = min(hi, max(lo, shi[k]))
                zhi = chengben(shi, xu, moshi)
                if zhi < zui - 1e-8:
                    zui = zhi
                    agent = shi
                    gai_le = True
        if not gai_le:
            xian = [q / 2.0 for q in xian]
            if max(xian) < 1e-6:
                break
    return agent, zui


def chu_shi_bu_fa(n, yao_k=False):
    bu_fa = [0.5] + [2.0] * n
    if yao_k:
        bu_fa += [0.0] * k_lie_chang      # 圈数是整数，精修不动它
    return bu_fa


def gongzhen_zhongzi(xu, ranges):
    """把若干 (n,k) 共振族种子摆进初始种群，DE 一起步就在正确的盆地里。"""
    out = []
    n = len(xu) - 1
    t0 = 0.5 * (ranges[0][0] + ranges[0][1])
    for i in k_tui:
        r, v = randv(xu[i], t0)
        P = xingqi(r, v)
        if not P:
            continue
        for nn in range(1, k_zui_da + 1):
            for kk in range(1, k_zui_da + 1):
                x = [0.0] * len(ranges)
                x[0] = t0
                for q in range(n):
                    lo, hi = ranges[1 + q]
                    x[1 + q] = 0.5 * (lo + hi)
                x[1 + i] = nn * P
                if x[1 + i] < ranges[1 + i][0] or x[1 + i] > ranges[1 + i][1]:
                    continue
                x[k_wei[i]] = kk
                out.append(x)
    return out[:12]


def sao_chuang(xu, kai, zhi, bu, kuan, pop, gen, lu_jing):
    """照着老版 --scan：发射窗口一片一片切，每片跑一轮短 DE + 精修，写进 csv。"""
    global riqi_qj
    hang = []
    riqi = kai
    while riqi <= zhi:
        riqi_qj = (riqi - kuan / 2.0, riqi + kuan / 2.0)
        random.seed(1)
        a, v = differential_evolution(pop, F, CR, gen, fanwei(xu), xu)
        a, v = jing_xiu(a, xu, chu_shi_bu_fa(len(xu) - 1))
        chao = ''
        for k in range(len(zucheng[2])):
            chao += ',%.0f' % zucheng[2][k][1]
        hang.append('%.2f,%.6f,%.4f,%.4f,%d%s,%.0f'
                    % (riqi, v, zucheng[0], zucheng[1], 1 if ke_xing else 0,
                       chao, riqi_biao[-1] - riqi_biao[0]))
        print('  窗口 %s  J=%.4f  %s'
              % (riqi_zh(riqi), v, '可行' if ke_xing else '不可行'))
        riqi += bu
    with open(lu_jing, 'w', encoding='utf-8') as wen:
        wen.write('riqi,J,vinf_bu,vinf_da,kexing,feiyue_gaodu_km,zong_tian\n')
        for q in hang:
            wen.write(q + '\n')
    print('扫完 %d 个窗口，已写 %s' % (len(hang), lu_jing))


# =========================================================================
# 第 8 块  打印与导出
# =========================================================================

def riqi_zh(jd):
    d = datetime.date(1970, 1, 1) + datetime.timedelta(days=int(jd - 2440587.5))
    return str(d)


def print_detail(agent, zong, xu):
    objective_function(agent, xu)
    n = len(xu) - 1
    hang = '%s %s' % (zhwen.get(xu[0], xu[0]), riqi_zh(riqi_biao[0]))
    for k in range(1, n + 1):
        hang += ' -> %s %s' % (zhwen.get(xu[k], xu[k]), riqi_zh(riqi_biao[k]))
    print(hang)
    print('  分支: %s   一共 %.0f 天'
          % (','.join(fenzhi_lie), riqi_biao[-1] - riqi_biao[0]))
    print('  出发 v_inf %.3f   到达 v_inf %.3f   J %.4f   %s'
          % (zucheng[0], zucheng[1], zong,
             '可行(零脉冲)' if ke_xing else '不可行'))
    for k in range(len(zucheng[2])):
        rp, alt, jia, jia_max, que = zucheng[2][k]
        print('    飞越%d %-4s  rp %9.1f km (高 %8.1f km)  δ %6.2f° / 最大 %6.2f°  缺陷 %8.3f m/s'
              % (k + 1, zhwen.get(xu[k + 1], xu[k + 1]), rp, alt,
                 jia, jia_max, que * 1000))
    print('  发射损失 %.3f + 捕获 %.3f = 真实总dv %.3f km/s（中间没有机动）'
          % (zucheng[-2], zucheng[-1], zucheng[-2] + zucheng[-1]))


def print_gongzhen(agent, J, xu):
    """共振轮结果单独打印：明确标出命中误差和缺陷，不冒充可行解。"""
    ping_gongzhen(agent, xu)
    d = gong_mingxi
    if d is None:
        print('  共振轮没算出可用结果')
        return
    hang = '%s %s' % (zhwen.get(xu[0], xu[0]), riqi_zh(d['ep'][0]))
    for k in range(1, len(d['ep'])):
        hang += ' -> %s %s' % (zhwen.get(xu[k], xu[k]), riqi_zh(d['ep'][k]))
    print(hang)
    print('  分支: %s   一共 %.0f 天   圈数基因 %s'
          % (','.join(d['fenzhi']), d['ep'][-1] - d['ep'][0],
             ','.join(str(q) for q in d['k_lie'])))
    print('  出发 v_inf %.3f   到达 v_inf %.3f   J %.4f   %s'
          % (d['vd'], d['va'], J,
             '可行' if d['feasible'] else '近似解(带残差，别当可行轨道用)'))
    print('    最大飞越缺陷 %.1f m/s   最大命中误差 %.0f km'
          % (d['max_defect'] * 1000, d['max_hit']))
    for k, (rp, alt, jia, jia_max, que) in enumerate(d['fei_lie']):
        print('    飞越%d %-4s  rp %9.1f km (高 %8.1f km)  δ %6.2f°  缺陷 %8.3f m/s'
              % (k + 1, zhwen.get(xu[k + 1], xu[k + 1]), rp, alt, jia, que * 1000))
    print('  发射损失 %.3f + 捕获 %.3f = 真实总dv %.3f km/s（中间没有机动）'
          % (d['dv_fa'], d['dv_da'], d['dv_fa'] + d['dv_da']))
    print('  注意：这一轮 J 里含缺陷罚，别跟普通轮的 J 直接比大小；'
          '它列出来是为了让你看见 Lambert 看不见的共振链')


def tian_of(agent, n):
    return [agent[1 + k] for k in range(n)]


def export_nodes(agent, xu, path):
    """写成 step6/step7/step8 通用的 nodes JSON（纯 Lambert 普通轮最优解）。"""
    objective_function(agent, xu)
    n = len(xu) - 1
    mu_sun, ra_sun = muandra(0)
    r_lie, v_lie = [], []
    for k in range(n + 1):
        r, v = randv(xu[k], riqi_biao[k])
        r_lie.append(r)
        v_lie.append(v)
    fen_lie = []
    for k in range(n):
        fen_lie.append(izzo_fenzhi(mu_sun, r_lie[k], r_lie[k + 1],
                                   tian_of(agent, n)[k] * 86400, M_zui_da))
    legs = []
    for k in range(n):
        for (v1, v2, biao) in fen_lie[k]:
            if biao == fenzhi_lie[k]:
                legs.append((v1, v2))
                break
    nodes = dict()
    nodes['departure'] = dict(body=names_en.get(xu[0], str(xu[0])),
                              epoch_mjd=float(riqi_biao[0] - 2400000.5),
                              vinf_km_s=[float(q) for q in (legs[0][0] - v_lie[0])])
    fei = []
    for k in range(1, n):
        mu_p, ra_p = muandra(xu[k])
        vin_in = legs[k - 1][1] - v_lie[k]
        vin_out = legs[k][0] - v_lie[k]
        rp_min = ra_p + gaodu_xian.get(xu[k], 200.0)
        ft = feiyue_qing(vin_in, vin_out, mu_p, rp_min)
        fei.append(dict(body=names_en.get(xu[k], str(xu[k])),
                        epoch_mjd=float(riqi_biao[k] - 2400000.5),
                        vinf_in_km_s=[float(q) for q in vin_in],
                        vinf_out_km_s=[float(q) for q in vin_out],
                        rp_km=float(ft['rp']),
                        alt_km=float(ft['rp'] - ra_p),
                        turn_deg=float(math.degrees(ft['jia'])),
                        turn_max_deg=float(math.degrees(ft['jia_max'])),
                        defect_m_s=float(ft['que'] * 1000)))
    nodes['flybys'] = fei
    nodes['dsms'] = [None] * n
    nodes['arrival'] = dict(epoch_mjd=float(riqi_biao[n] - 2400000.5),
                            target=names_en.get(xu[n], str(xu[n])))
    nodes['meta'] = dict(source='MGA_Planner main4 (纯MGA, 零脉冲借力, 无DSM)',
                         t0_mjd2000=float(riqi_biao[0] - 2400000.5 - 51544.5),
                         tof_days=[float(t) for t in tian_of(agent, n)],
                         objective=float(chengben(agent, xu)),
                         vinf_dep_km_s=float(zucheng[0]),
                         vinf_arr_km_s=float(zucheng[1]),
                         branch_ids=[str(q) for q in fenzhi_lie],
                         max_flyby_defect_m_s=float(max(
                             [f['defect_m_s'] for f in fei], default=0.0)),
                         dv_launch_km_s=float(zucheng[-2]),
                         dv_capture_km_s=float(zucheng[-1]),
                         feasible=bool(ke_xing),
                         x=[float(q) for q in agent])
    with open(path, 'w', encoding='utf-8') as fh:
        json.dump(nodes, fh, indent=1, ensure_ascii=False)
    return nodes, path


# =========================================================================
# 第 10 块 开跑（跟老 main3.py 一样，改完参数直接 python main4.py）
# =========================================================================

# =========================================================================
# 第 9 块  Lambert 求解器
# -------------------------------------------------------------------------
# 这一段是老仓库 izzo.py 原样搬过来的（Izzo 2015 算法），只用到 izzo2015
# 这一个入口。它是纯数值活儿，跟借力搜索没关系，不用看、也不用改。
# 搬进来是为了让 main4.py 真的只有这一个文件：拷走它就能跑。
# =========================================================================

def hypergeometric_f(z, tol):
    Sj = 1.0
    Cj = 1.0
    err = 1.0
    Cj1 = 0.0
    Sj1 = 0.0
    j = 0
    while err > tol:
        Cj1 = Cj * (3.0 + j) * (1.0 + j) / (2.5 + j) * z / (j + 1)
        Sj1 = Sj + Cj1
        err = abs(Cj1)
        Sj = Sj1
        Cj = Cj1
        j += 1
    return Sj


def izzo2015(
    mu,
    r1,
    r2,
    tof,
    M=0,
    prograde=True,
    low_path=True,
    maxiter=35,
    atol=1e-5,
    rtol=1e-7,
    full_output=False,
):
    r"""
    Solves Lambert problem using Izzo's devised algorithm.

    Parameters
    ----------
    mu: float
        Gravitational parameter, equivalent to :math:`GM` of attractor body.
    r1: numpy.array
        Initial position vector.
    r2: numpy.array
        Final position vector.
    M: int
        Number of revolutions. Must be equal or greater than 0 value.
    prograde: bool
        If `True`, specifies prograde motion. Otherwise, retrograde motion is imposed.
    low_path: bool
        If two solutions are available, it selects between high or low path.
    maxiter: int
        Maximum number of iterations.
    atol: float
        Absolute tolerance.
    rtol: float
        Relative tolerance.
    full_output: bool
        If True, the number of iterations is also returned.

    Returns
    -------
    v1: numpy.array
        Initial velocity vector.
    v2: numpy.array
        Final velocity vector.
    numiter: list
        Number of iterations.

    Notes
    -----
    This is the algorithm devised by Dario Izzo[1] in 2015. It inherits from
    the one developed by Lancaster[2] during the 60s, following the universal
    formulae approach. It is one of the most modern solvers, being a complete
    Lambert's problem solver (zero and Multiple-revolution solutions). It shows
    high performance and robustness while requiring no more than four iterations
    to reach a solution.

    All credits of the implementation go to Juan Luis Cano Rodríguez and the
    poliastro development team, from which this routine inherits. Some changes
    were made to adapt it to `lamberthub` API. In addition, the hypergeometric
    function is the one from SciPy.

    Copyright (c) 2012-2021 Juan Luis Cano Rodríguez and the poliastro development team

    References
    ----------
    [1] Izzo, D. (2015). Revisiting Lambert’s problem. Celestial Mechanics
           and Dynamical Astronomy, 121(1), 1-15.

    [2] Lancaster, E. R., & Blanchard, R. C. (1969). A unified form of
           Lambert's theorem (Vol. 5368). National Aeronautics and Space
           Administration.

    """

   

    # Chord
    c = r2 - r1
    c_norm, r1_norm, r2_norm = norm(c), norm(r1), norm(r2)

    # Semiperimeter
    s = (r1_norm + r2_norm + c_norm) * 0.5

    # Versors
    i_r1, i_r2 = r1 / r1_norm, r2 / r2_norm
    i_h = cross(i_r1, i_r2)
    i_h = i_h / norm(i_h)

    # Geometry of the problem
    ll = np.sqrt(1 - min(1.0, c_norm / s))

    # Compute the fundamental tangential directions
    if i_h[2] < 0:
        ll = -ll
        i_t1, i_t2 = cross(i_r1, i_h), cross(i_r2, i_h)
    else:
        i_t1, i_t2 = cross(i_h, i_r1), cross(i_h, i_r2)

    # Correct transfer angle parameter and tangential vectors regarding orbit's
    # inclination
    ll, i_t1, i_t2 = (
        (-ll, -i_t1, -i_t2) if prograde is False else (ll, i_t1, i_t2)
    )

    # Non dimensional time of flight
    T = np.sqrt(2 * mu / s**3) * tof

    # Find solutions and filter them
    x, y, numiter, tpi = _find_xy(ll, T, M, maxiter, atol, rtol, low_path)

    # Reconstruct
    gamma = np.sqrt(mu * s / 2)
    rho = (r1_norm - r2_norm) / c_norm
    sigma = np.sqrt(1 - rho**2)

    # Compute the radial and tangential components at initial and final
    # position vectors
    V_r1, V_r2, V_t1, V_t2 = _reconstruct(
        x, y, r1_norm, r2_norm, ll, gamma, rho, sigma
    )

    # Solve for the initial and final velocity
    v1 = V_r1 * (r1 / r1_norm) + V_t1 * i_t1
    v2 = V_r2 * (r2 / r2_norm) + V_t2 * i_t2

    return (v1, v2, numiter, tpi) if full_output is True else (v1, v2)


def _reconstruct(x, y, r1, r2, ll, gamma, rho, sigma):
    """Reconstruct solution velocity vectors."""
    V_r1 = gamma * ((ll * y - x) - rho * (ll * y + x)) / r1
    V_r2 = -gamma * ((ll * y - x) + rho * (ll * y + x)) / r2
    V_t1 = gamma * sigma * (y + ll * x) / r1
    V_t2 = gamma * sigma * (y + ll * x) / r2
    return [V_r1, V_r2, V_t1, V_t2]


def _find_xy(ll, T, M, maxiter, atol, rtol, low_path):
    """Computes all x, y for given number of revolutions."""
    # For abs(ll) == 1 the derivative is not continuous
    assert abs(ll) < 1

    M_max = np.floor(T / pi)
    T_00 = np.arccos(ll) + ll * np.sqrt(1 - ll**2)  # T_xM

    # Refine maximum number of revolutions if necessary
    if T < T_00 + M_max * pi and M_max > 0:
        _, T_min = _compute_T_min(ll, M_max, maxiter, atol, rtol)
        if T < T_min:
            M_max -= 1

    # Check if a feasible solution exist for the given number of revolutions
    # This departs from the original paper in that we do not compute all solutions
    if M > M_max:
        raise ValueError("No feasible solution, try lower M!")

    # Initial guess
    x_0 = _initial_guess(T, ll, M, low_path)

    # Start Householder iterations from x_0 and find x, y
    x, numiter, tpi = _householder(x_0, T, ll, M, atol, rtol, maxiter)
    y = _compute_y(x, ll)

    return x, y, numiter, tpi


def _compute_y(x, ll):
    """Computes y."""
    return np.sqrt(1 - ll**2 * (1 - x**2))


def _compute_psi(x, y, ll):
    """Computes psi.

    "The auxiliary angle psi is computed using Eq.(17) by the appropriate
    inverse function"

    """
    if -1 <= x < 1:
        # Elliptic motion
        # Use arc cosine to avoid numerical errors
        return np.arccos(x * y + ll * (1 - x**2))
    elif x > 1:
        # Hyperbolic motion
        # The hyperbolic sine is bijective
        return np.arcsinh((y - x * ll) * np.sqrt(x**2 - 1))
    else:
        # Parabolic motion
        return 0.0


def _tof_equation(x, T0, ll, M):
    """Time of flight equation."""
    return _tof_equation_y(x, _compute_y(x, ll), T0, ll, M)


def _tof_equation_y(x, y, T0, ll, M):
    """Time of flight equation with externally computated y."""
    if M == 0 and np.sqrt(0.6) < x < np.sqrt(1.4):
        eta = y - ll * x
        S_1 = (1 - ll - x * eta) * 0.5
        Q = 4 / 3 * hypergeometric_f(S_1, 1e-11)
        T_ = (eta**3 * Q + 4 * ll * eta) * 0.5
    else:
        psi = _compute_psi(x, y, ll)
        T_ = np.divide(
            np.divide(psi + M * pi, np.sqrt(np.abs(1 - x**2))) - x + ll * y,
            (1 - x**2),
        )

    return T_ - T0


def _tof_equation_p(x, y, T, ll):
    # TODO: What about derivatives when x approaches 1?
    return (3 * T * x - 2 + 2 * ll**3 * x / y) / (1 - x**2)


def _tof_equation_p2(x, y, T, dT, ll):
    return (3 * T + 5 * x * dT + 2 * (1 - ll**2) * ll**3 / y**3) / (
        1 - x**2
    )


def _tof_equation_p3(x, y, _, dT, ddT, ll):
    return (
        7 * x * ddT + 8 * dT - 6 * (1 - ll**2) * ll**5 * x / y**5
    ) / (1 - x**2)


def _compute_T_min(ll, M, maxiter, atol, rtol):
    """Compute minimum T."""
    if ll == 1:
        x_T_min = 0.0
        T_min = _tof_equation(x_T_min, 0.0, ll, M)
    else:
        if M == 0:
            x_T_min = np.inf
            T_min = 0.0
        else:
            # Set x_i > 0 to avoid problems at ll = -1
            x_i = 0.1
            T_i = _tof_equation(x_i, 0.0, ll, M)
            x_T_min = _halley(x_i, T_i, ll, atol, rtol, maxiter)
            T_min = _tof_equation(x_T_min, 0.0, ll, M)

    return [x_T_min, T_min]


def _initial_guess(T, ll, M, low_path):
    """Initial guess."""
    if M == 0:
        # Single revolution
        T_0 = np.arccos(ll) + ll * np.sqrt(1 - ll**2) + M * pi  # Equation 19
        T_1 = 2 * (1 - ll**3) / 3  # Equation 21
        if T >= T_0:
            x_0 = (T_0 / T) ** (2 / 3) - 1
        elif T < T_1:
            x_0 = 5 / 2 * T_1 / T * (T_1 - T) / (1 - ll**5) + 1
        else:
            # This is the real condition, which is not exactly equivalent
            # elif T_1 < T < T_0
            x_0 = (T_0 / T) ** (np.log2(T_1 / T_0)) - 1

        return x_0
    else:
        # Multiple revolution
        x_0l = (((M * pi + pi) / (8 * T)) ** (2 / 3) - 1) / (
            ((M * pi + pi) / (8 * T)) ** (2 / 3) + 1
        )
        x_0r = (((8 * T) / (M * pi)) ** (2 / 3) - 1) / (
            ((8 * T) / (M * pi)) ** (2 / 3) + 1
        )

        # Filter out the solution
        x_0 = (
            np.max([x_0l, x_0r]) if low_path is True else np.min([x_0l, x_0r])
        )

        return x_0


def _halley(p0, T0, ll, atol, rtol, maxiter):
    """Find a minimum of time of flight equation using the Halley method.

    Note
    ----
    This function is private because it assumes a calling convention specific to
    this module and is not really reusable.

    """
    for ii in range(1, maxiter + 1):
        y = _compute_y(p0, ll)
        fder = _tof_equation_p(p0, y, T0, ll)
        fder2 = _tof_equation_p2(p0, y, T0, fder, ll)
        if fder2 == 0:
            raise RuntimeError("Derivative was zero")
        fder3 = _tof_equation_p3(p0, y, T0, fder, fder2, ll)

        # Halley step (cubic)
        p = p0 - 2 * fder * fder2 / (2 * fder2**2 - fder * fder3)

        if abs(p - p0) < rtol * np.abs(p0) + atol:
            return p
        p0 = p

    raise RuntimeError("Failed to converge")


def _householder(p0, T0, ll, M, atol, rtol, maxiter):
    """Find a zero of time of flight equation using the Householder method.

    Note
    ----
    This function is private because it assumes a calling convention specific to
    this module and is not really reusable.

    """

    # The clock starts together with the iteration
    tic = time.perf_counter()
    for numiter in range(1, maxiter + 1):
        y = _compute_y(p0, ll)
        fval = _tof_equation_y(p0, y, T0, ll, M)
        T = fval + T0
        fder = _tof_equation_p(p0, y, T, ll)
        fder2 = _tof_equation_p2(p0, y, T, fder, ll)
        fder3 = _tof_equation_p3(p0, y, T, fder, fder2, ll)

        # Householder step (quartic)
        p = p0 - fval * (
            (fder**2 - fval * fder2 / 2)
            / (fder * (fder**2 - fval * fder2) + fder3 * fval**2 / 6)
        )

        if abs(p - p0) < rtol * np.abs(p0) + atol:
            # Stop the clock and compute the time per iteration
            tac = time.perf_counter()
            tpi = (tac - tic) / numiter

            return p, numiter, tpi
        p0 = p

    raise RuntimeError("Failed to converge")


if __name__ == '__main__':

    # 命令行给个序列码就临时换序列： python main4.py EEJ
    if len(sys.argv) > 1:
        xulie = sys.argv[1]

    if lie_yi_xia:
        ma_lie = lie_xulie(mubiao, zui_da_tui, qidian)
        print('%s(%d) → %s(%d)：%d 条 ≤%d 段路线'
              % (zhwen[qidian], qidian, zhwen[mubiao], mubiao,
                 len(ma_lie), zui_da_tui))
        print('前 20 条：', '  '.join(ma_lie[:20]))
        sys.exit(0)

    names = jie_ma(xulie, qidian)
    if not names:
        print('序列写错了：%s' % xulie)
        sys.exit(1)
    she_zhi(names)

    print('=' * 72)
    print('序列 %s   同体腿（自动识别，会搜共振）: %s'
          % (shuo_xulie(names),
             ','.join('%d(%s→%s)' % (i, zhwen[names[i]], zhwen[names[i + 1]])
                      for i in k_tui) if k_tui else '无'))
    print('星历: %s' % ('de405.bsp (jplephem SPK)' if _kernel is not None
                       else '内置 JPL 低精度星历（无外部文件）'))
    print('=' * 72)

    # ---------- 第一轮：普通轮，纯 Lambert，老流程原样 ----------
    agent_best, best_J = None, 1e30
    for lun in range(chong_qi):
        random.seed(sui_ji_zhong + lun)
        agent_lun, dv_lun = differential_evolution(
            pop_size, F, CR, generations, fanwei(names), names)
        agent_lun, dv_lun = jing_xiu(agent_lun, names,
                                     chu_shi_bu_fa(len(names) - 1))
        print('----- 普通轮 第%d轮  J = %.6f（精修后）' % (lun + 1, dv_lun))
        if dv_lun < best_J:
            agent_best, best_J = agent_lun, dv_lun
    print()
    if best_J < 99999:
        print('【普通轮 · 纯 Lambert】J = %.4f' % best_J)
        print_detail(agent_best, best_J, names)
        nodes, lu_jing = export_nodes(agent_best, names, 'main4_nodes.json')
        print('已写出 %s' % lu_jing)
    else:
        print('【普通轮】这个窗口里纯 Lambert 一条解都没算出来，跳过打印和导出')
        print('（换个发射窗口、放宽 TOF，或者把 M_zui_da 调大再试）')

    # ---------- 第二轮：共振轮，同体腿自动打靶 ----------
    if gongzhen and k_tui:
        print()
        print('----- 共振轮：同体腿 %s 改打靶，按 (n,k) 共振族播种 -----'
              % ','.join(str(i) for i in k_tui))
        ranges_g = fanwei(names, yao_k=True)
        zhong = gongzhen_zhongzi(names, ranges_g)
        random.seed(sui_ji_zhong + 99)
        agent_g, J_g = differential_evolution(
            gong_pop, F, CR, gong_gen, ranges_g, names,
            zhong_zi=zhong, moshi='gongzhen')
        agent_g, J_g = jing_xiu(agent_g, names,
                                chu_shi_bu_fa(len(names) - 1, yao_k=True),
                                moshi='gongzhen')
        print()
        print('【共振轮 · 引力弹弓自动搜共振腿】J = %.4f' % J_g)
        print_gongzhen(agent_g, J_g, names)

    if sao:
        print()
        print('----- 开始扫发射窗口 -----')
        sao_chuang(names, 2447360.5, 2449447.5, 90.0, 150.0, 40, 40,
                   'main4_sao.csv')
