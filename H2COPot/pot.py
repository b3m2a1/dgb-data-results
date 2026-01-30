import os, numpy as np
from McUtils.Extensions import DynamicFFILibrary, CLoader # we'll load this at runtime

__all__ = [
    "BasePotential",
    "InternalsPotential",
    "Potential"
]

cur_dir = os.path.dirname(os.path.abspath(__file__))
lib_file = os.path.join(cur_dir, "libs", "libh2co.so")
lib_dir = os.path.join(cur_dir, "libs", 'EPAPS')
if not os.path.isfile(lib_file):
    CLoader(cur_dir).custom_make(
        True, lib_dir
    )

# lib_file = os.path.join(test_dir, "libmbpol.so")
BasePotential = DynamicFFILibrary(
    lib_file,
    initialize=dict(
        name='ref_initialize_',
        flag=[int],
        defaults={'flag': -1},
        prep_args=lambda kw: [
            kw.__setitem__('flag', np.array([-1])),
            kw
        ][-1],
        return_type=int,
        return_handler=lambda r, kw: kw['flag']
    ),
    get_pot=dict(
        name='ref_pot_',
        coords=[float],
        energy=[float],
        return_type=float,
        prep_args=lambda kw: [
            kw.__setitem__('energy', np.zeros(len(kw['coords']) if kw['coords'].ndim > 1 else 1)),
            kw][-1],
        defaults={'energy': None},
        return_handler=lambda r, kw: kw['energy']
    ),
    abinit_initialize=dict(
        name='abinit_initialize_',
        flag=[int],
        defaults={'flag': -1},
        return_type=int,
        prep_args=lambda kw: [
            kw.__setitem__('flag', np.array([-1])),
            kw
        ][-1],
        return_handler=lambda r, kw: kw['flag']
    ),
    abinit_pot=dict(
        name='abinit_pot_',
        coords=[float],
        energy=[float],
        return_type=float,
        prep_args=lambda kw: [
            kw.__setitem__('energy', np.zeros(len(kw['coords']) if kw['coords'].ndim > 1 else 1)),
            kw][-1],
        defaults={'energy': None},
        return_handler=lambda r, kw: kw['energy']
    ),
    compiler_options=dict(
        threaded=True
    )
)

class InternalsPotential:
    @classmethod
    def get_pot(cls, coords, *, model='ref', debug=False, threading_mode='omp'):
        cd = os.getcwd()
        try:
            os.chdir(lib_dir)
            base_shape = coords.shape[:-1]
            coords = coords.reshape(-1, coords.shape[-1])

            abs_dihed = np.abs(coords[:, 5])
            wrap_pos = abs_dihed > np.pi #- 1e-8
            # print(wrap_pos)
            if wrap_pos.any():
                coords = coords.copy()
                coords[wrap_pos, 5] = -np.sign(coords[wrap_pos, 5]) * (2*np.pi - abs_dihed[wrap_pos])
                abs_dihed = np.abs(coords[:, 5])
                # wrap_pos = abs_dihed > np.pi  # - 1e-8
            #     print("--->", wrap_pos)

            if model == 'ref':
                BasePotential.initialize()
                vals = BasePotential.get_pot(coords,
                                             threading_vars=['energy', 'coords'],
                                             threading_mode=threading_mode,
                                             debug=debug)
            else:
                BasePotential.abinit_initialize()
                vals = BasePotential.abinit_pot(coords,
                                                threading_vars=['energy', 'coords'],
                                                threading_mode=threading_mode,
                                                debug=debug)

            if base_shape == ():
                return vals[0]
            else:
                return vals.reshape(base_shape)
        finally:
            os.chdir(cd)


class Potential:
    @classmethod
    def from_internals(cls, ints):
        import McUtils.Numputils as nput

        nd = ints.shape[:-1]
        C = np.zeros(nd + (3,))
        O = np.zeros(nd + (3,))
        O[..., 0] = ints[..., 0]
        Y = np.zeros(nd + (3,))
        Y[..., 1] = 1
        H1 = nput.cartesian_from_rad(C, O, Y, ints[..., 1], -ints[..., 3], None)[0]
        H2 = nput.cartesian_from_rad(C, O, H1, ints[..., 2], -ints[..., 4], ints[..., 5])[0]

        return np.moveaxis(np.array([H1, H2, C, O]), 0, -2)

    @classmethod
    def get_internals(cls, coords, *, order=None):
        import McUtils.Numputils as nput
        from McUtils.Data import UnitsData

        if order is None:
            order = [0, 1, 2, 3]
        coords = np.asanyarray(coords)
        H1, H2, C, O = [
            coords[..., i, :]
            for i in order
        ]
        rCH1 = nput.pts_norms(H1, C)
        rCH2 = nput.pts_norms(H2, C)
        rOC = nput.pts_norms(C, O)
        aOCH1 = nput.pts_angles(H1, C, O)[0]
        aOCH2 = nput.pts_angles(H2, C, O)[0]
        dOCHH = nput.pts_dihedrals(H1, O, C, H2)

        internals = np.moveaxis(
            np.array([rOC, rCH1, rCH2, aOCH1, aOCH2, dOCHH]),
            0, -1
        )
        return internals

    @classmethod
    def get_pot(cls, coords, *, order=None, model='ref', debug=False):
        from McUtils.Data import UnitsData

        internals = cls.get_internals(coords, order=order)

        return InternalsPotential.get_pot(internals, model=model, debug=debug)

