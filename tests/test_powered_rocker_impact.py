import numpy as np

from crane.models.rocker_compass import (
    RockerCompassParams,
    heelstrike_map as passive_heelstrike,
    kinetic_energy,
)
from crane.models.powered_rocker_compass import (
    PoweredRockerCompassParams,
    make_powered_rocker_compass,
    powered_heelstrike_map,
)

NOMINAL = dict(m=1.0, m_h=0.0, c=0.37, rho=0.32, L=1.0, R=0.3, gamma=0.030, g=9.81)


def test_pushoff_zero_equals_passive():
    pp = PoweredRockerCompassParams(**NOMINAL, push_off=0.0)
    pa = RockerCompassParams(**NOMINAL)
    for x in [
        np.array([0.30844, -0.30844, -1.26256, -0.87914]),
        np.array([0.25, -0.25, -1.1, -0.8]),
    ]:
        np.testing.assert_allclose(
            powered_heelstrike_map(x, pp), passive_heelstrike(x, pa), atol=1e-12
        )


def test_pushoff_injects_energy():
    pa = RockerCompassParams(**NOMINAL)
    x = np.array([0.30844, -0.30844, -1.26256, -0.87914])
    pp0 = PoweredRockerCompassParams(**NOMINAL, push_off=0.0)
    ppP = PoweredRockerCompassParams(**NOMINAL, push_off=0.2)
    ke0 = kinetic_energy(powered_heelstrike_map(x, pp0), pa)
    keP = kinetic_energy(powered_heelstrike_map(x, ppP), pa)
    assert keP > ke0


def test_factory_lift_project():
    pp = PoweredRockerCompassParams(**NOMINAL, push_off=0.1)
    model = make_powered_rocker_compass(pp)
    y = np.array([0.3, -1.2, -0.85])
    np.testing.assert_allclose(model.project(model.lift(y)), y, atol=1e-12)
