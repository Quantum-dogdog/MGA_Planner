import numpy as np
from get_mu_ra import muandra
import math

def fashedv(shuzi,v_need,v_now,gaodu):

    v_inf_out = np.linalg.norm(v_need - v_now)   #结果是数
    mu ,ra = muandra(shuzi)
    r_park = ra + gaodu
    v_park = math.sqrt(mu / r_park)
    v_after = math.sqrt(v_inf_out * v_inf_out +2 * mu/r_park)
    dv = v_after - v_park
    return dv

def buhuodv(shuzi,v_arrive,v_now,gaodu):
    
    mu ,ra = muandra(shuzi)
    r_park = ra + gaodu
    v_inf_in = np.linalg.norm(v_arrive - v_now)
    a_arrive = - mu / (v_inf_in * v_inf_in)
    v_inpe = math.sqrt(2 * mu / r_park - mu / a_arrive)
    v_escape = math.sqrt(2 * mu / r_park)           #逃逸速度
    dv = v_inpe - v_escape
    return dv

def dsmdv(v_after,v_before):
    
    dv = np.linalg.norm(v_after - v_before)
    return dv


'''


def flybyhou(shuzi,v_in,v_now,gaodu,jiao):
    
    mu ,ra = muandra(shuzi)
    v_inf_in = v_in - v_now
    v1 = np.linalg.norm(v_inf_in)  #v1是数
    s_bar = v_inf_in / v1
    t = np.cross(v_now,v_inf_in)
    t_bar = t / np.linalg.norm(t) 
    cos_theta = np.cos(jiao)        #jiao是弧度
    sin_theta = np.sin(jiao)
    u = s_bar
    b_bar = np.array([[cos_theta + u[0] * u[0] * (1 - cos_theta), u[0] * u[1] * (1 - cos_theta) - u[2] * sin_theta, u[0] * u[2] * (1 - cos_theta) + u[1] * sin_theta],
                  [u[1] * u[0] * (1 - cos_theta) + u[2] * sin_theta, cos_theta + u[1] * u[1] * (1 - cos_theta), u[1] * u[2] * (1 - cos_theta) - u[0] * sin_theta],
                  [u[2] * u[0] * (1 - cos_theta) - u[1] * sin_theta, u[2] * u[1] * (1 - cos_theta) + u[0] * sin_theta, cos_theta + u[2] * u[2] * (1 - cos_theta)]])
    b_bar = np.dot(b_bar,t_bar)   #旋转矩阵在前，需要转的矢量在后
    n = np.cross(b_bar,s_bar)
    n_bar = n / np.linalg.norm(n)
    r_pe = ra + gaodu 
    jiao2 = 2*math.asin(1/(1+(r_pe*v1*v1/mu)))
    cos_theta = np.cos(jiao2)        #jiao是弧度
    sin_theta = np.sin(jiao2)
    u = n_bar
    v_out =  np.array([[cos_theta + u[0] * u[0] * (1 - cos_theta), u[0] * u[1] * (1 - cos_theta) - u[2] * sin_theta, u[0] * u[2] * (1 - cos_theta) + u[1] * sin_theta],
                  [u[1] * u[0] * (1 - cos_theta) + u[2] * sin_theta, cos_theta + u[1] * u[1] * (1 - cos_theta), u[1] * u[2] * (1 - cos_theta) - u[0] * sin_theta],
                  [u[2] * u[0] * (1 - cos_theta) - u[1] * sin_theta, u[2] * u[1] * (1 - cos_theta) + u[0] * sin_theta, cos_theta + u[2] * u[2] * (1 - cos_theta)]])
    v_inf_out = np.dot(v_out,v_in)
    v_out = v_inf_out + v_now
    return v_out



'''

def flybyhou(target_body, v_in, v_now, altitude, incidence_angle):
    # 获取目标天体的引力常数和半径
    mu, ra = muandra(target_body)
    
    # 计算相对速度
    v_inf_in = v_in - v_now
    v_inf_in_magnitude = np.linalg.norm(v_inf_in)  # 相对速度的大小
    
    # 单位化相对速度向量
    s_hat = v_inf_in / v_inf_in_magnitude
    
    # 计算垂直于相对速度和当前速度的向量
    t = np.cross(v_now, v_inf_in)
    t_hat = t / np.linalg.norm(t)
    
    # 计算旋转矩阵所需的三角函数值
    cos_theta = np.cos(incidence_angle)  # incidence_angle 应该是弧度
    sin_theta = np.sin(incidence_angle)
    
    # 构造旋转矩阵
    rotation_matrix = np.array([
        [cos_theta + s_hat[0]**2 * (1 - cos_theta), s_hat[0] * s_hat[1] * (1 - cos_theta) - s_hat[2] * sin_theta, s_hat[0] * s_hat[2] * (1 - cos_theta) + s_hat[1] * sin_theta],
        [s_hat[1] * s_hat[0] * (1 - cos_theta) + s_hat[2] * sin_theta, cos_theta + s_hat[1]**2 * (1 - cos_theta), s_hat[1] * s_hat[2] * (1 - cos_theta) - s_hat[0] * sin_theta],
        [s_hat[2] * s_hat[0] * (1 - cos_theta) - s_hat[1] * sin_theta, s_hat[2] * s_hat[1] * (1 - cos_theta) + s_hat[0] * sin_theta, cos_theta + s_hat[2]**2 * (1 - cos_theta)]
    ])
    
    # 旋转垂直向量
    b_hat = np.dot(rotation_matrix, t_hat)
    
    # 计算新的垂直向量
    n = np.cross(b_hat, s_hat)
    n_hat = n / np.linalg.norm(n)
    
    # 计算近心点距离
    r_periapsis = ra + altitude
    
    # 计算飞越角
    jiao2 = 2 * np.arcsin(1 / (1 + (r_periapsis * v_inf_in_magnitude**2 / mu)))
    
    # 使用飞越角构造旋转矩阵
    cos_phi = np.cos(jiao2)
    sin_phi = np.sin(jiao2)
    rotation_matrix_out = np.array([
        [cos_phi + n_hat[0]**2 * (1 - cos_phi), n_hat[0] * n_hat[1] * (1 - cos_phi) - n_hat[2] * sin_phi, n_hat[0] * n_hat[2] * (1 - cos_phi) + n_hat[1] * sin_phi],
        [n_hat[1] * n_hat[0] * (1 - cos_phi) + n_hat[2] * sin_phi, cos_phi + n_hat[1]**2 * (1 - cos_phi), n_hat[1] * n_hat[2] * (1 - cos_phi) - n_hat[0] * sin_phi],
        [n_hat[2] * n_hat[0] * (1 - cos_phi) - n_hat[1] * sin_phi, n_hat[2] * n_hat[1] * (1 - cos_phi) + n_hat[0] * sin_phi, cos_phi + n_hat[2]**2 * (1 - cos_phi)]
    ])
    
    # 计算飞越后的相对速度
    v_inf_out = np.dot(rotation_matrix_out, v_inf_in)
    
    # 计算飞越后的绝对速度
    v_out = v_inf_out + v_now
    
    return v_out

  
