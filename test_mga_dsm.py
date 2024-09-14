from mgadsm import mga_dsm

#riqi = 2461339.5
riqi = 2459053.5      #发射日期



#xulie = [3,4,5]

xulie = [3,4,5,6,7]   #引力弹弓总序列
#v1 = [28.66504442,14.52063624,7.55242545]
#v1 = [-20.29,24.5,12.49]
v1 = [28.6,12.5,7.5]
tian = [202,500,800,1600]        #每条腿间飞行总时间
xishu = [0.5,0.43,0.3,0.3]       #深空机动发生时刻占整条腿飞行时间的系数
gaodu = [200,200,2000,2000,2000]    #轨道高度距离星球表面
jiao = [1.5,2.1,0.1]    #飞越b平面上矢量b与T_bar的夹角



jieguo = mga_dsm(xulie,riqi, v1, tian, xishu, gaodu, jiao)


print('结果：',jieguo)

xulie = [3,4,5]
v1 = [28.6,12.5,7.5]
tian = [202,500]        #每条腿间飞行总时间
xishu = [0.5,0.43]       #深空机动发生时刻占整条腿飞行时间的系数
gaodu = [200,200,2000]    #轨道高度距离星球表面
jiao = [1.5] 

jieguo = mga_dsm(xulie,riqi, v1, tian, xishu, gaodu, jiao)


print('结果：',jieguo)



