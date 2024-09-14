# MGA_Planner

本项目全名Multiple Gravity Assist Planner（引力助推规划器），又名SlingShot Trajectory Optimizer（SSTO，引力弹弓轨迹优化器），目前属于简陋的发行版本。


# 介绍

本项目与网络上已有的mga planner相比，优点是：

（1）纯python，代码逻辑简单易懂，不和C++掺和，不用下载matlab

（2）免部署，no need setup，no need build，no need install，下载到本地之后，直接>>python main.py运行

（3）零依赖，除默认python库(指numpy)外，无需安装多个第三方python库。不会像别的项目一样，因为依赖库的更新、废弃而run error
