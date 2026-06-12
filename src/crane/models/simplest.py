"""Garcia et al. 1998 'The Simplest Walking Model' (m/M → 0 極限、l=g=1 スケール)。

状態 x = [θ, φ, θ̇, φ̇]。θ: stance 脚角度（slope 法線基準）、
φ: stance 脚から測った swing 脚相対角度。
heel-strike 条件は φ − 2θ = 0 (θ < 0 側)。
"""

from dataclasses import dataclass

import numpy as np


@dataclass(frozen=True)
class SimplestParams:
    gamma: float  # slope angle [rad]


def dynamics(t: float, x, p: SimplestParams):
    """連続相（swing phase）の運動方程式。"""
    theta, phi, theta_dot, _phi_dot = x
    theta_dd = np.sin(theta - p.gamma)
    phi_dd = theta_dd + theta_dot**2 * np.sin(phi) - np.cos(theta - p.gamma) * np.sin(phi)
    return [theta_dot, _phi_dot, theta_dd, phi_dd]


def heelstrike_map(x: np.ndarray) -> np.ndarray:
    """heel-strike の衝突写像（角運動量保存）+ 脚交換。"""
    theta, _phi, theta_dot, _phi_dot = x
    c = np.cos(2.0 * theta)
    return np.array([-theta, -2.0 * theta, c * theta_dot, c * (1.0 - c) * theta_dot])


def lift(theta: float, theta_dot: float) -> np.ndarray:
    """Poincaré 断面座標 y=(θ, θ̇) を全状態へ持ち上げる（衝突直後の拘束）。"""
    return np.array([theta, 2.0 * theta, theta_dot, (1.0 - np.cos(2.0 * theta)) * theta_dot])
