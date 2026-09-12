"""
Test script for ulambert.py (universal variable Lambert solver).

Same test cases as test_izzo.py, so the two solvers can be compared
directly. Only depends on numpy.
"""

import numpy as np

from ulambert import universal_lambert

np.set_printoptions(precision=6)  # 小数点后6位，方便和期望值对比


def run_case(name, mu, r1, r2, tof, eta=0, N=0,
             expected_v1=None, expected_v2=None, tol=1e-4):
    """Run one Lambert case, print results and (optionally) check them."""
    print("=" * 60)
    print(name)

    v1, v2, xout = universal_lambert(r1, r2, tof, mu, eta, N)

    # 单解情形：取第一列
    v1 = v1[:, 0]
    v2 = v2[:, 0]

    print("v1 =", v1)
    print("v2 =", v2)
    print("x  =", xout)

    if expected_v1 is not None:
        err1 = np.linalg.norm(v1 - expected_v1)
        print("err(v1) =", err1)
        assert err1 < tol, f"v1 mismatch: {v1} vs {expected_v1}"
    if expected_v2 is not None:
        err2 = np.linalg.norm(v2 - expected_v2)
        print("err(v2) =", err2)
        assert err2 < tol, f"v2 mismatch: {v2} vs {expected_v2}"
    print("-" * 60)


# ----------------------------------------------------------------------
# Case 1: Curtis, "Orbital Mechanics for Engineering Students", Example 5.2
# 地球引力场，76 分钟转移（椭圆解）
# ----------------------------------------------------------------------
mu_earth = 3.986004418e5  # [km^3 / s^2]
r1 = np.array([15945.34, 0.0, 0.0])           # [km]
r2 = np.array([12214.83899, 10249.46731, 0.0])  # [km]
tof = 76.0 * 60                                 # [s]

run_case(
    "Case 1: Curtis Example 5.2 (elliptic)",
    mu_earth, r1, r2, tof,
    expected_v1=np.array([2.058913, 2.915965, 0.0]),   # [km/s]
    expected_v2=np.array([-3.451565, 0.910315, 0.0]),  # [km/s]
)

# ----------------------------------------------------------------------
# Case 2: Curtis Example 5.2 之后的第二个算例（椭圆解）
# ----------------------------------------------------------------------
mu_earth = 3.986004418e5  # [km^3 / s^2]
r1 = np.array([5000.0, 10000.0, 2100.0])       # [km]
r2 = np.array([-14600.0, 2500.0, 7000.0])      # [km]
tof = 3600                                     # [s]

run_case(
    "Case 2: 3D elliptic transfer, tof = 1h",
    mu_earth, r1, r2, tof,
    expected_v1=np.array([-5.9925, 1.9254, 3.2456]),   # [km/s]
    expected_v2=np.array([-3.3125, -4.1966, -0.38529]),  # [km/s]
)

# ----------------------------------------------------------------------
# Case 3: 日心坐标系算例（和 test_izzo.py 相同，单位 AU / year）
# ----------------------------------------------------------------------
mu_sun = 39.47692641  # [AU^3 / year^2]
r1 = np.array([0.159321004, 0.579266185, 0.052359607])  # [AU]
r2 = np.array([0.057594337, 0.605750797, 0.068345246])  # [AU]
tof = 0.010794065                                       # [year]

run_case(
    "Case 3: heliocentric transfer (AU, year)",
    mu_sun, r1, r2, tof,
    expected_v1=np.array([-9.303603251, 3.018641330, 1.536362143]),  # [AU/year]
)

# ----------------------------------------------------------------------
# Case 4: 地球 -> 火星，2015.2.8 出发，飞行 300 天（日心坐标系）
# 和 test_izzo.py 一样，可与 DE405 星历求出的实际速度做对比
# （需要 de405.bsp 和 jplephem，若没有可注释掉这部分）
# ----------------------------------------------------------------------
try:
    from jplephem.spk import SPK

    mu_sun = 1.32712440018e11  # [km^3 / s^2]
    r1 = np.array([-1.10e8, 8.93e7, 3.87e7])   # [km] 地球 2015.2.8
    r2 = np.array([-1.14e8, 1.95e8, 9.23e7])   # [km] 火星 2015.2.8 + 200 days
    tof = 300 * 86400                           # [s]

    v1, v2, xout = universal_lambert(r1, r2, tof, mu_sun, eta=0, N=0)
    v1, v2 = v1[:, 0], v2[:, 0]

    print("=" * 60)
    print("Case 4: Earth -> Mars, 300 day transfer")
    print("v1 =", v1, "km/s  (日心坐标系速度)")
    print("v2 =", v2, "km/s")
    print("-" * 60)

    kernel = SPK.open("de405.bsp")

    jd = 2457261.5
    position, velocity = kernel[0, 3].compute_and_differentiate(2457061.5)
    velocity_per_second = velocity / 86400.0
    print("earth vx vy vz", velocity_per_second)  # 单位 km/s
    dv1 = np.linalg.norm(v1 - velocity_per_second)
    print("dv1 =", dv1)

    position, velocity = kernel[0, 4].compute_and_differentiate(jd)  # 4 是火星 barycenter
    velocity_per_second = velocity / 86400.0
    print("mars  vx vy vz", velocity_per_second)  # 单位 km/s
    dv2 = np.linalg.norm(v2 - velocity_per_second)
    print("dv2 =", dv2)
except ImportError:
    print("=" * 60)
    print("Case 4 skipped: jplephem / de405.bsp not available")

print("=" * 60)
print("All tests passed!")
