# MGA_Planner

本项目全名Multiple Gravity Assist Planner（引力助推规划器），又名SlingShot Trajectory Optimizer（SSTO，引力弹弓轨迹优化器），目前属于简陋的发行版本。


# 介绍

本项目与网络上已有的mga planner相比，优点是：

（1）纯python，代码逻辑简单易懂，不和C++掺和，不用下载matlab

（2）免部署，no need setup，no need build，no need install，下载到本地之后，直接>>python main.py运行

（3）零依赖，除默认python库(指numpy)外，无需安装多个第三方python库。不会像别的项目一样，因为依赖库的更新、废弃而run error


# 开发背景

2023年9月份，在探索使用krpc控制ksp（坎巴拉太空计划）中的火箭的过程中，我对行星际轨道转移、轨迹优化产生了浓厚的兴趣，经过搜索了解，我认识到行星际轨道转移、轨迹优化用专业术语来讲，是一个mga1dsm问题和一个全局优化问题。这是总体大局观上的认识。

在过去的一年里，我对网络上，尤其是github上有关引力弹弓（引力辅助）项目进行了广泛的学习，最终决定采用python自己造一下轮子，这样才算是从知道到达懂得。

下面是我搜集到的项目列表，在这里我按照功能全面度、文档丰富度对这些项目进行纯主观排序，供大家学习参考：

1.[pykep](https://github.com/esa/pykep)
由达里奥·伊佐（Dario Izzo）大佬及其团队开发的pykep，是我心目中毫无疑问的No.1.大佬本人还是兰伯特问题最高效求解算法的发明者，本项目中使用的就是izzo算法。
该项目功能十分全面，但是需要安装多个第三方python库，以及c++运行的一些依赖，使用pykep优化mga1dsm问题需要结合pygmo（也是大佬团队开发的），在windows系统中安装基本上都能成功。

2.[tudat-space](https://github.com/tudat-team/tudatpy)
由荷兰代尔夫特大学开发，同样需要安装多个第三方python库，以及c++运行的一些依赖，优化mga1dsm问题也需要结合pygmo，但是它在windows系统中安装不上，好在它部署在Binder上，你可以在线使用。

3.[ksptot](https://github.com/Arrowstar/ksptot)
由Arrowstar开发，与上述两个项目相比，大佬个人开发此项目且文档详实、赏心悦目，但是它需要你安装MATLAB Compiler Runtime R2022a才能用（需要占用你20GB的空间），若是你想用于太阳系内行星际转移的优化，还需要自己设置星球的参数（毕竟正如你所见到的，它是ksp的tot）。

4.[KSP-MGA-Planner](https://github.com/Krafpy/KSP-MGA-Planner)
由Krafpy开发，同样是个人开发者，该项目为网页部署，若你想要离线运行，需要下载安装node.js，然后npm安装http-server.它可以用于太阳系内行星际转移的优化，采用的是遗传算法来优化轨迹。

5.[astra](https://github.com/andreabellome/astra)
由Andrea Bellome开发，需要你安装整个MATLAB以及C语言的编译器，它采用的是动态规划来优化轨迹。

6.[AUTOMATE](https://github.com/HadrienAFSA/AUTOMATE)
同样是一个MATLAB项目，它用的是最经典的Tisserand图解法，所以相对来讲，在算法方面不那么高大上。

7.[EMTG](https://github.com/nasa/EMTG)
由Jacob Englander领衔开发，是一个属于NASA的项目，用于优化低推力轨迹，环境配置超麻烦，还需要snopt。

8.[SkippingStone](https://github.com/rodyo/FEX-SkippingStone)
MATLAB项目。

9.[IGATO](https://github.com/tingspain/IGATO)
是C++项目，引用的有pykep,无文档。

10.[Wayfinder](https://github.com/Muetdhiver-lab/Wayfinder)
引用pykep&pygmo，适用于ksp.

11.
