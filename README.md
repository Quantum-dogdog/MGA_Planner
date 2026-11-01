# MGA_Planner

本项目全名Multiple Gravity Assist Planner（引力助推规划器），又名SlingShot Trajectory Optimizer（SSTO，引力弹弓轨迹优化器），目前属于可用的初级版本（26-11-01）。



## 项目优点

与网络上已有的 mga planner 相比，优点是：

1. **纯 python**，代码逻辑简单易懂，不和 C++ 掺和，不用下载 matlab。
2. **免部署**，no need setup，no need build，no need install。下载到本地之后，在 cmd 里 cd 到解压后的文件夹内 >> `python main3.py`，或者用 python 自带的 IDLE 打开 main3.py。需要注意的是，de405.bsp 这个 10M 的星历表必须跟 py 文件在同一个文件夹里，因为代码里我用的是相对路径，换个目录打开就报找不到文件。
3. **只需要 numpy 这一个第三方库**，其余的 jplephem 我直接把源码放进来了，所以不会像别的项目一样因为某个依赖库更新了、废弃了就 run error。


## 这个仓库能干什么

给它一条借力序列，比如地球-金星-地球-地球-木星，再给一个发射日期窗口，它帮你搜出什么时候出发、每段飞多少天，让"出发超速 + 到达超速 + 一次小于200m/s的dv"尽量小，最后写成 `main3_nodes.json`，可以直接喂给我的[book_scripts_for_GMAT_NASA](https://github.com/Quantum-dogdog/book_scripts_for_GMAT_NASA)的step6。

现在有三个入口：

| 文件 | 是什么 | 建议 |
| --- | --- | --- |
| `main.py` | 最早那版，序列长度写死成 3，中间还会打深空机动 | 留着对照，不推荐跑 |
| `main2.py` | MGA-1DSM，中间允许打一次发动机，先 DE 粗搜再精修 | 想带 DSM 就跑这个 |
| `main3.py` | 纯借力，全程一次发动机都不打，只靠星球甩和最后捕获 | **当前主版本** |

注意：`main2.py` 和 `main3.py` 文件末尾我直接就写运行代码了，没写 `if __name__ == '__main__'`，所以**别去 import 它**，一 import 就开始跑优化，想当库用你自己加个守卫。

## 运行

```bash
pip install numpy
cd 解压后的目录
python main3.py
```

## 天体编号表



```text
0-太阳   1-空（水星借力用不到，没写）   2-金星   3-地球   4-火星
5-木星   6-土星   7-天王星   8-海王星
```

所以 `xulie = [3, 2, 3, 3, 5]` 读作地球→金星→地球→地球→木星。

单位约定：距离 km，速度 km/s，日期儒略日 JD（导出的 json 里是 MJD），main3 里的时间单位是天。

## 文件介绍

**1. `jplephem/`** 由 Brandon Rhodes 大佬开发，用于他的 skyfield 项目，这里我们依据 MIT license 对其进行引用，这样就不用再安装 SPICE 了。

**2. `de405.bsp`** JPL 行星星历表，10M，记得自己下载。

**3. `constant.py`** 定义了一系列常量，它们代表不同天体的标准引力参数（mu）和半径（r）。注意同样一套数字在 `get_mu_ra.py` 里我又抄了一份，两处数值现在是一样的。

**4. `get_mu_ra.py`** `muandra(shuzi)` 返回 (mu, ra)，按上面那张编号表取。mu 单位 km³/s²，ra 单位 km。除月球是平均半径，其他都是赤道半径。

**5. `get_r_v.py`** `randv(shuzi, riqi)` 返回某天体在日心坐标系下的位置和速度。jplephem 给的速度是 km/day，我除了 86400 换成 km/s。

**6. `erti.py`** 天体力学和轨道计算相关，它是一个类，用法是 `tb = Twobody()` 然后 `tb.xxx`：

1. `statetoelement`：位置和速度向量转轨道六要素（半长轴、轨道倾角、离心率、升交点赤经、近心点幅角、真近点角）。
2. `elementtostate`：轨道六要素转回位置和速度向量。
3. `getsma`：计算半长轴。
4. `eafromta`：真近点角转偏近点角。
5. `taafterdt`：给定时间间隔后的真近点角，t1 单位是天。

另外还有模块级的 `newtonRootSolve` 用来解开普勒方程。注意：


- `newtonRootSolve` 原来没有任何保护，遇到斜率是 0、或者数值溢出算出 nan 的死轨道会一直转不出来，整个优化就卡在那。现在加了斜率为 0、出 nan、迭代超过 100 次这三道保险，一律返回 nan，由上层把这个解丢掉。容差也从 1e-15 放回 1e-12，1e-15 在 double 里经常是取不到的，只是逼着循环跑到黑。



**7. `get_dv.py`** 一系列计算航天器进行不同类型机动时所需速度变化量（Δv）的函数：

1. `fashedv(shuzi, v_need, v_now, gaodu)`：从停泊轨道出发，加速到所需超速要多少 dv，考虑引力常数 mu、轨道高度 gaodu 和星球半径 ra。
2. `buhuodv(shuzi, v_arrive, v_now, gaodu, yuandian=0)`：捕获要减多少速度。原来那版是拿"减速到逃逸速度"当捕获，现在加了个 `yuandian` 参数，减速目标是近拱点在停泊轨道、远拱点在 yuandian 倍星球半径的那条捕获椭圆；写 0 就是直接进圆轨道。
3. `dsmdv(v_after, v_before)`：两个速度之间的差异，深空机动用。
4. `flybyhou(target_body, v_in, v_now, altitude, incidence_angle)`：给定飞越高度和 B 面转角算出口速度。文件里有两份 flybyhou，旧那份被 `'''` 包着不生效，生效的是英文参数这份。另外 **main3 用不到它**（main3 中间不打发动机），用的是自己写的 `feiyue_qing`。

**8. `izzo.py`** 用于解决 Lambert 问题，取自 pykep，依据 GPL-3.0 license 对其进行了微调，使其不再依赖外部 hypergeometric_f 函数。

**9. `ulambert.py` / `bijiao_izzo_ulambert.py`** 另一套多圈 Lambert 的写法，以及它跟 izzo.py 的对拍。

**10. `main.py`** 实现了一个差分进化算法，用于优化空间任务中的飞行总时间、深空机动时刻的比例系数、轨道高度、弹弓角度等参数。
上一版我在这里写的"如果你想写飞行序列为 4 或 5 的代码，你需要改 objective_function、agent、ranges"——**现在**，main2 和 main3 都是按 `len(xulie)` 循环写的，改序列只改一行 `xulie` 就行，别的不用动。

**11. `mgadsm.py`** 根据输入的序列长度调用 mga_dsm_2/3/4/5 算总 dv。你可以看到我写的代码是低端代码，远不如 KSP-MGA-Planner 中抽象、高等。

**12. `sequence.py`** 定义行星际飞行序列。

**13. `spice_kernel_reader.py`** 用于读取和解析 SPICE 核心文件，是 jplephem 原项目里的 tmp54.py。

**14. `main3_nodes.json`** main3 跑完自动写出的结果，字段见下面"输出格式"。

**15. `test_*.py`、`totaldv.py`、`tianwen1_total_dv.py`** 一些打印对拍的脚本，主要是拿天问一号的数据比一比。

## main3.py 为什么变量这么少

染色体只有 L+1 个数：出发日期 1 个，加上每段飞行时间 L 个（L 就是借力次数 = `len(xulie)-1`）。

以前 main.py 把飞越高度 gaodu、B 面转角 jiao、机动时刻系数 xishu 也当变量去搜，其实前后两段的方向一定下来，要转多少角 δ 就定了，贴着多高 rp 也就定了：

```text
rp        = mu_p / v_inf^2 * (1 / sin(delta/2) - 1)      # rp 是输出量，不是变量
delta_max = 2 * asin(1 / (1 + rp_min * v_inf^2 / mu_p))  # 这颗星球最多能给多少角
```

代价是：前后两段各自用 Lambert 解出来以后，进超速和出超速的**模长一般对不上**，转角落不进这颗星球能给的范围内。这个"对不上多少"我叫它飞越缺陷，当罚函数加进目标：

```text
J = 1.0 * vinf_dep + 1.0 * vinf_arr + 2000 * sum(缺陷^2)
```

所以 **J 不是总 dv**，别拿它跟 main.py / mgadsm.py 的总 dv 直接比。main3 另外用 fashedv 和 buhuodv 把"从停泊轨道发射"和"末端捕获"这两笔单独报出来，方便跟 main2 对照。缺陷小于 1e-3 km/s 记为可行（ke_xing）。

其他参数顺便写了：

| 名字 | 意思 | 现在的值 |
| --- | --- | --- |
| `M_zui_da` | 每条腿最多绕太阳几圈。圈数不进染色体，直接穷举（M × 顺行逆行 × 高路低路全解一遍，分支名标成 `M0顺低` 这种） | 1 |
| `liu_shao` | 束宽。分支组合是指数级的，每条腿只留最便宜的这么几个往下接 | 20 |
| `gaodu_xian` | 每颗星球最低飞越高度 | 金星 300，地球/火星 200，木星/土星 5000 km |
| `jiatuo` / `buhuo_gaodu` | 停泊轨道高度，发射和捕获都用它 | 200 km |
| `yuandian` | 捕获椭圆远拱点是星球半径的多少倍，任务要求不是设计变量 | 20 |
| `vinf_shang` | 火箭能给的超速上限（C3 上限的另一种写法），写 None 就是不管 | None |
| DE 参数 | 自己的差分进化 | pop=60, F=0.5, CR=0.9, gen=150, 重起 2 次 |
| `jing_xiu` | DE 完再精修一下。step8 用的是 scipy 的 SLSQP，我不能引新包，就自己写土办法：每个变量 ± 动一小步，谁能把 J 变小就留着，一整轮没改善就把步长砍一半，这叫模式搜索，十几行 | — |
| `sao_chuang` | 把发射日期窗口按天往前扫，每个起点各跑一次。`sao=True` 打开，很慢，但好窗口就是这么找的 | False |

## 输出格式（main3_nodes.json）

```text
departure : body, epoch_mjd, vinf_km_s（矢量）
flybys[]  : body, epoch_mjd, vinf_in_km_s, vinf_out_km_s,
            rp_km, alt_km, turn_deg（实际转角）, turn_max_deg（星球能给的上限）,
            defect_m_s（缺陷，单位 m/s）
dsms[]    : 纯 MGA 全为 null，占位给下游用
arrival   : epoch_mjd, target
meta      : J、两个超速、每段天数、每腿选中的分支号（如 "M0顺低"）
```

## main3 哪里对了

1. **rp 和转角不再是变量了。** 老代码搜 gaodu 和 jiao，等于把已经被物理定死的东西再扔进优化器搜一遍，费算力，还容易搜出前后速度对不上的假解。维度从 8 个降到 L+1 个。
2. **圈数不再塞进连续变量再取整了。** 取整一下函数值就跳变，那是在糊弄差分进化。现在是老老实实穷举。
3. **不可行这件事被显式建模了。** 老代码只加 dv，一条根本飞不成的弹道也能给出很低的 dv。main3 有罚项和可行门，J 小同时要求缺陷小。这是从"算 dv"到"算弹道"的区别。
4. **分支组合用束搜索压住了**，所以它实际能跑完，不是纸面正确。
5. **结果变成可复查的接口。** 每个节点带速度矢量、高度、实际转角和上限、缺陷。main.py 那版只 print 一个数。

## 现在跑出来的是什么

序列 E-V-E-E-J，窗口 1989-02-02 到 1989-11-09：

```text
出发      1989-06-24   vinf = 7.766 km/s
金星借力  1989-12-17   高度 32222 km，转 9.50°（上限 41.26°），缺陷 0.051 m/s
地球借力1 1990-09-29   高度  1939 km，转 42.76°（上限 49.72°），缺陷 0.168 m/s
地球借力2 ...          高度  1725 km，转 43.12°（上限 49.31°），缺陷 0.316 m/s
到达木星  1996-03-28   vinf = 5.817 km/s
四段飞行时间 176 / 286 / 823 / 1185 天，总 6.76 年
J = 13.583，所有腿选中分支都是 M0顺低
```



## 还没解决的问题

1. 两点 Lambert 在转移角接近 0 度的时候是退化的：过同一点、同 a 的双曲线有无穷多条，这种腿缺陷怎么都压不到 0，所以答案里没有共振借力（还没从[book_scripts_for_GMAT_NASA](https://github.com/Quantum-dogdog/book_scripts_for_GMAT_NASA)copy过来）。




该段介绍由智谱清言拟稿、Qoder 接手开发，我润色。


# 开发背景

2023年9月份，在探索使用krpc控制ksp（坎巴拉太空计划）中的火箭的过程中，我对行星际轨道转移、轨迹优化产生了浓厚的兴趣，经过搜索了解，我认识到行星际轨道转移、轨迹优化用专业术语来讲，是一个mga1dsm问题和一个全局优化问题。这是总体大局观上的认识。

在过去的一年里，我对网络上，尤其是github上有关引力弹弓（引力辅助）项目进行了广泛的学习，最终决定采用python自己造一下轮子，这样才算是“从知道到达懂得”（此句来自北大2011年微电影《女生日记》，女主江小夏）。

下面是我搜集到的项目列表，在这里我按照功能全面度、文档丰富度对这些项目进行纯主观排序，供大家学习参考：

1.[pykep](https://github.com/esa/pykep)
由达里奥·伊佐（Dario Izzo）大佬及其团队开发的pykep，是我心目中毫无疑问的No.1.大佬本人还是兰伯特问题最高效求解算法的发明者，本项目中使用的就是izzo算法。
该项目功能十分全面，但是需要安装多个第三方python库，以及c++运行的一些依赖，使用pykep优化mga1dsm问题需要结合pygmo（也是大佬团队开发的），在windows系统中安装基本上都能成功，但是有一些例子已经dead。

2.[tudat-space](https://github.com/tudat-team/tudatpy)
由荷兰代尔夫特大学开发，同样需要安装多个第三方python库，以及c++运行的一些依赖，优化mga1dsm问题也需要结合pygmo，但是它在windows系统中安装不上，好在它部署在Binder上，你可以在线使用。

3.[ksptot](https://github.com/Arrowstar/ksptot)
由Arrowstar开发，与上述两个项目相比，大佬个人开发此项目且文档详实、赏心悦目，但是它需要你安装MATLAB Compiler Runtime R2022a才能用（需要占用你20GB的空间），若是你想用于太阳系内行星际转移的优化，还需要自己设置星球的参数（毕竟正如你所见到的，它是ksp的tot）。

4.[KSP-MGA-Planner](https://github.com/Krafpy/KSP-MGA-Planner)
由Krafpy开发，同样是个人开发者，该项目为网页部署，若你想要离线运行，需要下载安装node.js，然后npm安装http-server.它可以用于太阳系内行星际转移的优化，但是星球数据是用的截取的固定的六根数，不是实时读取spice kernal,采用的是遗传算法来优化轨迹。

5.[astra](https://github.com/andreabellome/astra)
由Andrea Bellome开发，需要你安装整个MATLAB以及C语言的编译器，它采用的是动态规划来优化轨迹。

6.[Aerospace Trajectory Optimization](https://sourceforge.net/projects/aero-trajectory-optimization/files/)
由cdeaglejr开发,github上仅上传了pdf,代码在sourceforge和dropbox上,有MATLAB、fortran多个版本，例子是优化的emj飞行序列。

7.[AUTOMATE](https://github.com/HadrienAFSA/AUTOMATE)
同样是一个MATLAB项目，它用的是最经典的Tisserand图解法，所以相对来讲，在算法方面不那么高大上。

8.[EMTG](https://github.com/nasa/EMTG)
由Jacob Englander领衔开发，是一个属于NASA的项目，用于优化低推力轨迹，环境配置超麻烦，还需要snopt，花钱。

9.[SkippingStone](https://github.com/rodyo/FEX-SkippingStone)
是MATLAB项目。

10.[MGALT-STOpS](https://github.com/jpcaldwell01/PI_MGALT_STOpS)
原本是MATLAB项目，向python迁移中，所以是个混合项目,用于优化低推力轨迹。

11.[MOLTO-IT](https://github.com/uc3m-aerospace/MOLTO-IT)
是MATLAB项目，用于优化低推力轨迹。

12.[TrajOpt_DE](https://github.com/padmanabhapsimha/TrajOpt_DE)
是C++、MATLAB混合项目。

13.[IGATO](https://github.com/tingspain/IGATO)
是C++项目，引用的有pykep,无文档。

14.[LTGA](https://github.com/Lpyshmily/LTGA)
是C、C++混合项目，用于优化低推力轨迹。

15.[Wayfinder](https://github.com/Muetdhiver-lab/Wayfinder)
引用pykep&pygmo，适用于ksp.

16.[GAToptimization](https://github.com/tomginsberg/GAToptimization)
超级混合项目，现实生活中不重复造轮子，到处引用的样子。

17.[Optimum_Interplanetary_Trajectory](https://github.com/AdamHibberd/Optimum_Interplanetary_Trajectory)
是MATLAB项目，需要SPICE、NOMAD，没有文档，默认你会安装SPICE,NOMAD是一个和snopt差不多的东西，花钱.

18.[hodographic-shaping-method-python](https://github.com/lstubbig/hodographic-shaping-method-python)
使用全息图优化低推力轨迹，尚不能看懂。

19.[Lunar-Gravity-Assist-Patched-Conic-Simulation](https://github.com/MatthewSlater12/Lunar-Gravity-Assist-Patched-Conic-Simulation)
地月系统内引力弹弓。

20.[Lunar-Gravity-Assist](https://github.com/swtnikita50/Lunar-Gravity-Assist)
月球弹弓。

21.[Gravity_Assist_Planner](https://github.com/Yourself1011/Gravity_Assist_Planner)
需要安装Processing这个小众软件，可视化演示中既没有Assist,也没有Planner，认定为toy。

22.[gravity-assist-flyby-optimizer](https://github.com/itchono/gravity-assist-flyby-optimizer)
还没有开发到轨迹优化的部分，属于烂尾项目。

23.[pythonflyby](https://github.com/kepler-69c/pythonFlyby)
效果演示项目，代码简单，认定为toy.

24.[Trajectory_Optimization_Python](https://github.com/camstillo/Trajectory_Optimization_Python)
需要GEKKO这个python包，并不含引力弹弓。

25.[gravity assist](https://github.com/search?q=gravity%20assist&type=repositories)
github上大量重名的项目，认定为平平无奇的张三。

