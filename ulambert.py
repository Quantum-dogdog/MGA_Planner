"""
A universal Lambert solver supporting single- and multi-revolution transfers.

Python translation of the MATLAB routine 'universalLambert.m'.

References
----------
[1] McElreath, J., Down, I. M., Majji, M., "A universal approach for solving
    the multi-revolution Lambert's Problem", Celestial Mechanics and
    Dynamical Astronomy, 2025.
"""

import warnings

import numpy as np


def universal_lambert(r1, r2, dt, mu, eta, N):
    r"""
    Solves the Lambert problem using the universal variable formulation.

    Parameters
    ----------
    r1: numpy.array
        Initial position vector [3,] (units: distance).
    r2: numpy.array
        Final position vector [3,] (units: distance).
    dt: float
        Time of flight (units: time).
    mu: float
        Gravitational parameter of the central body (units: distance^3/time^2).
    eta: int
        Transfer direction: ``0`` for short way, ``1`` for long way.
    N: int
        Number of revolutions (non-negative integer).

    Returns
    -------
    v1: numpy.array
        Initial velocity vector(s), shape ``(3, M)`` where M is the number of
        solutions found (use ``v1[:, 0]`` for the single-solution case).
    v2: numpy.array
        Final velocity vector(s), shape ``(3, M)``.
    xout: numpy.array
        Final converged solution(s) of the solving parameter, shape ``(M,)``.

    Notes
    -----
    Supports elliptic, parabolic and hyperbolic trajectories, as well as
    multi-revolution solutions. Complex arithmetic is used internally to
    handle the hyperbolic initial guesses, exactly as in the original
    MATLAB code.
    """
    # Max iterations and tolerances
    maxit = 25
    cnv1 = 1e-10
    cnv2 = 1e-2

    # Error checks
    if eta not in (0, 1):
        raise ValueError(
            "eta must equal zero or one, eta = 0 for short-way and eta = 1 for long-way."
        )
    if N % 1 != 0 or N < 0:
        raise ValueError("N must be an integer and greater than or equal to zero.")
    if dt < 0:
        raise ValueError("Transfer time must be greater than zero.")
    N = int(N)

    # Constants
    sqmu = np.sqrt(mu)
    r1n = np.linalg.norm(r1)
    r2n = np.linalg.norm(r2)
    r1pr2 = r1n + r2n
    r1tr2 = r1n * r2n
    ir1 = np.asarray(r1, dtype=float) / r1n
    ir2 = np.asarray(r2, dtype=float) / r2n
    c = np.linalg.norm(np.asarray(r2, dtype=float) - np.asarray(r1, dtype=float))
    ic = (np.asarray(r2, dtype=float) - np.asarray(r1, dtype=float)) / c
    TwoPI = 2 * np.pi
    theta = TwoPI * eta + (-1) ** eta * np.arccos(np.dot(r1, r2) / (r1n * r2n))
    cosHalfTh = np.cos(0.5 * theta)
    costheta = 2 * cosHalfTh**2 - 1
    sintheta = np.sin(theta)
    TwoPIN = TwoPI * N
    m = (-1) ** N * np.sqrt(r1tr2) * cosHalfTh
    Ttil = dt * sqmu

    # Initial guess for the solution
    if N == 0:
        # Single revolution solution
        M = 1
        chi0 = np.sqrt(2 * (r1pr2 - 2 * m))  # parabolic transfer time parameter
        T0 = chi0 * m + chi0**3 / 6
        if T0 < Ttil:
            # Elliptic solution: interpolate for initial guess
            Tth = _tlamb(theta, r1pr2, m, 0)[0]
            beta = np.sqrt((Ttil - T0) / (Tth - T0)) * theta / (TwoPI - theta)
            xi = TwoPI * beta / (1 + beta)
            bds = [0.0, TwoPI]
        elif T0 > Ttil:
            # Hyperbolic solution: interpolate for initial guess
            if theta < np.pi:
                tau = r1pr2**2 / (2 * m**2) - 1
                xtau = 1j * np.log(tau + np.sqrt(tau**2 - 1 + 0j))
                x1 = xtau * np.exp(-2.5)
                T1 = _tlamb(x1, r1pr2, m, 0)[0]
                xi = xtau * (1 - (Ttil / T0) ** 1.6) ** (
                    2.5 / (np.log(-(T0**1.6) + 0j) - np.log(T1**1.6 - T0**1.6))
                )
                bds = [0.0, xtau.imag]
            else:  # theta >= pi
                Tth = _tlamb(1j * theta, r1pr2, m, 0)[0]
                xi = 1j * theta * np.sqrt(np.log(Ttil / T0) / np.log(Tth / T0))
                bds = [0.0, np.inf]
        else:
            # Parabolic solution found directly
            v1, v2 = _recon_vel(
                0.0, r1pr2, r1tr2, c, ir1, ir2, ic, m, sqmu, costheta, sintheta
            )
            return v1, v2, np.array([0.0])
    else:
        # Multi-revolution solution
        M = 2
        sqN = np.sqrt(N)
        bds = [TwoPIN, TwoPIN + TwoPI]
        if theta < np.pi:
            xm0 = 2 * (theta + (1 / 6) ** (2.5 * sqN)) ** 0.4 - (1 / 3) ** sqN + TwoPIN
        else:
            xm0 = (
                -2 * (-theta + TwoPI + 0.05 ** np.sqrt(N)) ** 0.4
                + (TwoPI - 0.25**sqN)
                + TwoPIN
            )

        # Solve for the multi-rev minimum transfer time (Halley iterations)
        x = xm0
        for k in range(1, maxit + 1):
            Tm, Tp, Tpp, Tppp = _tlamb(x, r1pr2, m, 1)
            dx = -2 * Tp * Tpp / (2 * Tpp**2 - Tp * Tppp)
            xold = x
            x = x + dx
            if x > bds[1] or x < bds[0]:
                # Use Newton if dx jumps bounds
                x = xold - Tp / Tpp
            err = abs((x - xold) / (x - TwoPIN + TwoPI))
            dT = (Ttil - Tm) / Ttil
            if dT < 1e-4 and cnv2 != 1e-8:
                cnv2 = 1e-8
            if err < cnv2:
                break
        xm = x
        Tm = _tlamb(xm, r1pr2, m, 0)[0]

        # Check if converged
        if k == maxit:
            warnings.warn(
                "Failed to converge within max number of iterations when "
                "solving for the multi-rev minimum transfer time."
            )

        if Tm < Ttil:
            # Two multi-rev solutions exist
            xth = theta + TwoPIN
            N1 = N + 1
            Tth = _tlamb(xth, r1pr2, m, 0)[0]
            d = (xm - xth) ** 2 / (
                abs(Tth - Tm)
                * (TwoPIN - xth) ** 2
                * (TwoPIN + TwoPI - xth) ** 2
            )
            gam = -np.sqrt((Ttil - Tm) * d)
            b_s = -TwoPI * (N + N1) * gam + 1
            c_s = TwoPI**2 * N * N1 * gam - xm
            b_l = TwoPI * (N + N1) * gam + 1
            c_l = -(TwoPI**2) * N * N1 * gam - xm
            xi = [
                np.real((-b_s + np.sqrt(b_s**2 - 4 * gam * c_s)) / (2 * gam)),
                np.real((-b_l + np.sqrt(b_l**2 + 4 * gam * c_l)) / (-2 * gam)),
            ]
        elif Tm == Ttil:
            # One multi-rev solution found
            v1, v2 = _recon_vel(
                xm, r1pr2, r1tr2, c, ir1, ir2, ic, m, sqmu, costheta, sintheta
            )
            return v1, v2, np.array([xm])
        else:
            raise ValueError(
                "Requested Number of Revolutions is not feasible for the given "
                "transfer time. Increase transfer time or decrease number of "
                "revolutions."
            )

    # Solve for the desired transfer times and reconstruct velocities
    xi = np.atleast_1d(xi)
    xout = np.zeros(M, dtype=complex)
    v1 = np.zeros((3, M))
    v2 = np.zeros((3, M))

    for i in range(M):
        x = xi[i]
        for k in range(1, maxit + 1):
            T, Tp, Tpp, Tppp = _tlamb(x, r1pr2, m, 1)
            T = T - Ttil
            # Householder step (quartic)
            dx = -(6 * T * Tp**2 - 3 * T**2 * Tpp) / (
                6 * Tp**3 - 6 * T * Tp * Tpp + T**2 * Tppp
            )
            xold = x
            x = x + dx
            # Use lower order update if jumps bounds (for robustness)
            if x.real + x.imag > bds[1] or x.real + x.imag < bds[0]:
                # Halley step
                x = xold - 2 * T * Tp / (2 * Tp**2 - T * Tpp)
                if x.real + x.imag > bds[1] or x.real + x.imag < bds[0]:
                    # Newton step
                    x = xold - T / Tp
                    if x.real + x.imag > bds[1] or x.real + x.imag < bds[0]:
                        sgn = np.sign((-T / Tp).real + (-T / Tp).imag)
                        x = xold + 0.5 * (bds[1 + int(sgn > 0)] - xold)
            err = abs((x - xold) / (x - TwoPIN + TwoPI))
            if err < cnv1:
                break
        xout[i] = x
        if k == maxit:
            warnings.warn(
                "Failed to converge within max number of iterations when "
                "solving for the desired solution."
            )
        v1[:, i], v2[:, i] = _recon_vel(
            x, r1pr2, r1tr2, c, ir1, ir2, ic, m, sqmu, costheta, sintheta
        )

    return v1, v2, xout


def _tlamb(x, r1pr2, m, flag):
    """
    Computes the scaled transfer time and its 1st, 2nd, 3rd derivatives
    with respect to the universal variable x.
    """
    Tp = Tpp = Tppp = None

    # Stumpff-like auxiliary quantities (trig evaluation)
    U0 = np.cos(x)
    U0s = np.cos(0.5 * x)
    al = (1 - U0) / (r1pr2 - 2 * m * U0s)
    ial = 1 / al
    sqal = np.sqrt(al)
    isqal = 1 / sqal
    U1 = np.sin(x) * isqal
    U1s = np.sin(0.5 * x) * isqal
    U2 = (1 - U0) * ial
    U3 = (x * isqal - U1) * ial

    T = 2 * U1s * m + U3

    if flag:
        # Derivatives
        iU1s = 1 / U1s
        ap = 0.5 * iU1s * sqal * (2 * U0s - m * al)
        app = -al * 0.5 + (0.5 * U0s * isqal - m * sqal) * ap * iU1s
        appp = (
            (
                0.25 * m * iU1s**2 * (m * al - 0.5 * U0s)
                - 0.5 * isqal * iU1s * (U0s * ial * 0.5 + m) * ap
                - 0.75
            )
            * ap
            + (0.5 * isqal * iU1s) * (U0s - 2 * m * al) * app
        )

        Tp = ial * ((m * U0s + U2) * sqal - 0.5 * (T + 2 * U3) * ap)
        Tpp = ial * (
            (-3 * Tp + 2 * m * U0s * isqal) * ap
            - 0.75 * T * ial * ap**2
            + U1
            - m * al * U1s * 0.5
            - (0.5 * T + U3) * app
        )
        Tppp = ial * (
            U0 * isqal
            - m * U0s * sqal * 0.25
            - ap
            * (
                4.5 * Tpp
                + 1.5 * m * U1s
                + (2.25 * Tp * ial - 0.375 * ial**2 * T * ap) * ap
            )
            + (
                2 * m * U0s * isqal
                - U2 * isqal
                - 3.5 * Tp
                + (U3 * ial - 1.75 * T * ial) * ap
            )
            * app
            - (0.5 * T + U3) * appp
        )

    return T, Tp, Tpp, Tppp


def _recon_vel(x, r1pr2, r1tr2, c, ir1, ir2, ic, m, sqmu, costheta, sintheta):
    """Reconstructs the initial and final velocity vectors from x.

    x may be complex (hyperbolic branch uses complex arithmetic,
    as in the original MATLAB code); the physical velocities are real,
    so any residual imaginary part is discarded here explicitly.
    """
    U0s = np.cos(0.5 * x)
    sqrtp = np.sqrt(r1tr2 * (1 - costheta) / (r1pr2 - 2 * m * U0s))
    vc = sqmu * c * sqrtp / (r1tr2 * sintheta)
    vp = sqmu * (1 - costheta) / (sqrtp * sintheta)
    v1 = vc * ic + vp * ir1
    v2 = vc * ic - vp * ir2

    # 双曲线分支下 x 为复数，物理解的虚部应为 0（仅剩数值残余）
    if max(np.abs(v1.imag).max(), np.abs(v2.imag).max()) > 1e-8:
        warnings.warn("Velocity imaginary part is not negligible: "
                      f"|Im v1|={np.abs(v1.imag).max():.3e}, "
                      f"|Im v2|={np.abs(v2.imag).max():.3e}")
    
    return np.real(v1), np.real(v2)   # <-- 显式取实部


if __name__ == "__main__":
    # Example: Curtis "Orbital Mechanics for Engineering Students", Example 5.2
    # Earth transfer, r1 -> r2 in 76 min (elliptic, single rev)
    mu_earth = 398600.0  # km^3/s^2
    r1 = np.array([15945.34, 0.0, 0.0])  # km
    r2 = np.array([12214.83899, 10249.46731, 0.0])  # km
    tof = 76.0 * 60.0  # s

    v1, v2, xout = universal_lambert(r1, r2, tof, mu_earth, eta=0, N=0)
    print("v1 =", v1[:, 0], "km/s")  # expected ~ [2.0589, 2.9159, 0]
    print("v2 =", v2[:, 0], "km/s")  # expected ~ [-3.7456, 0.4911, 0]
    print("x  =", xout)
