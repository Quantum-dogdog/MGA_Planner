# MGA_Planner

本项目全名Multiple Gravity Assist Planner（引力助推规划器），又名SlingShot Trajectory Optimizer（SSTO，引力弹弓轨迹优化器），目前属于简陋的发行版本。


# 项目优点

本项目与网络上已有的mga planner相比，优点是：

  （1）纯python，代码逻辑简单易懂，不和C++掺和，不用下载matlab

  （2）免部署，no need setup，no need build，no need install，下载到本地之后，在cmd里cd到解压后的文件夹内 >> python main.py运行或者用python自带的IDLE打开main.py

  （3）零依赖，除默认python库(指numpy)外，无需安装多个第三方python库。不会像别的项目一样，因为依赖库的更新、废弃而run error

# 文件介绍

1.[jplephem](https://github.com/brandon-rhodes/python-jplephem)由Brandon Rhodes大佬开发，用于他的[skyfield](https://github.com/skyfielders/python-skyfield/)项目，这里我们依据MIT license对其进行引用，这样就不用再安装SPICE了。

2.constant.py定义了一系列常量，它们代表不同天体的标准引力参数（mu）和半径（r）等。

3.erti.py 文件包含了一些与天体力学和轨道计算相关的功能。它包含以下方法：

  （1）statetoelement: 将位置和速度向量转换为轨道六要素（半长轴、轨道倾角、离心率、升交点赤经、近心点幅角、真近点角）。
       
  （2）elementtostate: 将轨道六要素转换回位置和速度向量。
       
  （3）getsma: 计算半长轴。
       
  （4）taafterdt: 计算在给定时间间隔后的真近点角。

4.get_dv.py 文件包含了一系列计算航天器进行不同类型机动时所需速度变化量（Δv）的函数。以下是文件的主要组成部分：

  （1）fashedv 函数：计算从当前速度 v_now 加速到所需速度 v_need 所需的速度变化量。这个函数考虑了引力常数 mu、轨道高度 gaodu 和目标天体的半径 ra。
       
  （2）buhuodv 函数：计算从当前速度 v_now 减速到所需速度 v_arrive 所需的速度变化量。这个函数也考虑了引力常数 mu、轨道高度 gaodu 和目标天体的半径 ra。
       
  （3）dsmdv 函数：计算两个速度之间的差异。
       
  （4）flybyhou 函数：计算航天器在经过目标天体附近时所需的速度变化量。这个函数考虑了目标天体的引力常数 mu、半径 ra、飞越高度 altitude 和入射角 incidence_angle。

5.get_mu_ra.py 文件用于计算不同天体的引力常数（μ）和半径（r）。

6.get_r_v.py 文件用于计算特定日期和天体位置的距离（r）和速度（v）。

7.izzo.py 文件用于解决 Lambert 问题，取自[pykep](https://github.com/esa/pykep)，依据GPL-3.0 license对其进行了微调，使其不再依赖外部hypergeometric_f函数。

8.main.py 文件实现了一个差分进化算法，用于优化空间任务中的飞行总时间、深空机动时刻的比例系数、轨道高度、弹弓角度等参数，没有优化发射日期、发射序列、发射速度等参数，但是也可以了不是吗。
需要注意的是，目前main.py文件里我只写了飞行序列为3的代码，如果你想写飞行序列为4或5的代码，你需要改objective_function、agent、ranges。

9.mgadsm.py 文件可以根据输入的序列长度（表示空间任务的复杂性），调用不同的子函数（mga_dsm_2, mga_dsm_3, mga_dsm_4, mga_dsm_5: 这些子函数分别处理长度为2、3、4、5的序列）来计算总速度变化量（delta-v），你可以看到我写的代码是低端代码，远不如[KSP-MGA-Planner](https://github.com/Krafpy/KSP-MGA-Planner)中抽象、高等。

10.sequence.py 文件定义行星际飞行序列。

11.spice_kernal_reader.py 文件用于读取和解析 SPICE 核心文件，是[jplephem](https://github.com/brandon-rhodes/python-jplephem)原项目中的tmp54.py。

12.剩下的是一些例子，主要是用天问1号的数据对结果进行交叉验证。

该段介绍由[智谱清言](https://chatglm.cn/)拟稿，我润色。

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

