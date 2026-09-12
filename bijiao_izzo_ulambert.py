"""
Side-by-side comparison of izzo2015 and universal_lambert.
"""
import numpy as np
from izzo import izzo2015
from ulambert import universal_lambert

np.set_printoptions(precision=6)


def pick_eta(r1, r2, prograde):
    """根据几何自动选择 eta，使其与 izzo 的 prograde 对应。

    短程转移的角动量方向 = r1 x r2：
      - 若 (r1 x r2)_z > 0，短程即顺行 -> prograde=True  对应 eta=0
      - 若 (r1 x r2)_z < 0，短程为逆行 -> prograde=True  对应 eta=1
    """
    hz = np.cross(r1, r2)[2]
    short_way_is_prograde = hz > 0
    return 0 if (short_way_is_prograde == prograde) else 1


def compare(name, mu, r1, r2, tof, prograde=True, low_path=True, N=0):
    r1 = np.asarray(r1, dtype=float)
    r2 = np.asarray(r2, dtype=float) 

    eta = pick_eta(r1, r2, prograde)

    v1i, v2i = izzo2015(mu, r1, r2, tof,
                        prograde=prograde, low_path=low_path, M=N)
    v1u, v2u, xout = universal_lambert(r1, r2, tof, mu, eta, N)
    v1u, v2u = v1u[:, 0], v2u[:, 0]

    dv1 = np.linalg.norm(v1i - v1u)
    dv2 = np.linalg.norm(v2i - v2u)

    print("=" * 62)
    print(f"{name}  (eta = {eta})")
    print("izzo      v1 =", v1i)
    print("ulambert  v1 =", v1u)
    print("|dv1| =", dv1)
    print("izzo      v2 =", v2i)
    print("ulambert  v2 =", v2u)
    print("|dv2| =", dv2)
    return dv1, dv2


mu_earth = 3.986004418e5   # km^3/s^2
mu_sun = 1.32712440018e11  # km^3/s^2

compare("Case 1: Curtis Ex 5.2", mu_earth,
        [15945.34, 0, 0], [12214.83899, 10249.46731, 0], 76 * 60)

compare("Case 2: 3D elliptic", mu_earth,
        [5000, 10000, 2100], [-14600, 2500, 7000], 3600)

compare("Case 3: heliocentric (km)", mu_sun,
        [-1.10e8, 8.93e7, 3.87e7], [-1.14e8, 1.95e8, 9.23e7], 0.010794065 * 365.25 * 86400,
        prograde=False)   # 注意：Case 3 在 test_izzo 里没指定 prograde，默认 prograde=True
                          # 但该几何短程即为顺行方向；如果对不上可切换 prograde 试试

compare("Case 4: Earth -> Mars, 300 d", mu_sun,
        [-1.10e8, 8.93e7, 3.87e7], [-1.14e8, 1.95e8, 9.23e7], 300 * 86400,
        prograde=True)    # 之前对不上的用例，自动映射 eta=1
