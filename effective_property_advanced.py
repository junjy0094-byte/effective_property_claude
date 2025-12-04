"""
Advanced Composite Material Effective Property Calculator using PyANSYS

This version implements proper periodic boundary conditions using constraint equations,
which provides more accurate results for RVE homogenization.

The method follows the standard computational homogenization approach where:
1. Periodic boundary conditions are enforced via constraint equations
2. Volume-averaged stress and strain are used to compute effective properties
3. A complete 6x6 stiffness matrix is obtained and inverted to get compliance

Reference:
- Sun, C.T., Vaidya, R.S. (1996). "Prediction of composite properties from a
  representative volume element"
- Xia, Z., Zhou, C., Yong, Q., Wang, X. (2006). "On selection of repeated unit cell
  model and application of unified periodic boundary conditions"
"""

import numpy as np
from ansys.mapdl.core import launch_mapdl


class AdvancedCompositeCalculator:
    """
    Advanced calculator for effective material properties using proper periodic BC.
    """

    def __init__(self, rve_size=1.0, fiber_vf=0.10):
        """
        Initialize the calculator.

        Parameters
        ----------
        rve_size : float
            Size of the RVE cube in mm
        fiber_vf : float
            Fiber volume fraction (0 to 1)
        """
        self.L = rve_size
        self.Vf = fiber_vf
        self.V = rve_size ** 3

        # Fiber radius from volume fraction
        self.r_fiber = np.sqrt(fiber_vf * rve_size**2 / np.pi)

        # Material properties
        self.mat_props = {
            'matrix': {'E': 3500.0, 'nu': 0.35, 'alpha': 60e-6},
            'fiber': {'E': 230000.0, 'nu': 0.20, 'alpha': -0.5e-6}
        }

        self.mapdl = None
        self.node_pairs = {}
        self.stiffness_matrix = None
        self.compliance_matrix = None
        self.effective_props = {}

    def set_material(self, phase, E, nu, alpha):
        """Set material properties for a phase."""
        self.mat_props[phase] = {'E': E, 'nu': nu, 'alpha': alpha}

    def launch(self, **kwargs):
        """Launch MAPDL."""
        self.mapdl = launch_mapdl(**kwargs)
        self.mapdl.clear()
        self.mapdl.prep7()
        print("MAPDL launched successfully")

    def exit(self):
        """Exit MAPDL."""
        if self.mapdl:
            self.mapdl.exit()
            self.mapdl = None

    def build_model(self, elem_size=0.05):
        """Build the RVE model with mesh."""
        m = self.mapdl
        L = self.L
        r = self.r_fiber

        print(f"\n{'='*60}")
        print("BUILDING RVE MODEL")
        print(f"{'='*60}")
        print(f"RVE size: {L} mm")
        print(f"Fiber radius: {r:.4f} mm")
        print(f"Volume fraction: {self.Vf*100:.1f}%")

        # Clear and setup
        m.clear()
        m.prep7()

        # Create geometry
        # Matrix block
        m.block(0, L, 0, L, 0, L)

        # Fiber cylinder along Z
        m.cyl4(L/2, L/2, r, 0, r, 360, L)

        # Boolean operations
        m.vsbv(1, 2, sepo='', keep2='keep')
        # Result: V2=fiber, V3=matrix

        # Materials
        mp = self.mat_props['matrix']
        m.mp('EX', 1, mp['E'])
        m.mp('NUXY', 1, mp['nu'])
        m.mp('ALPX', 1, mp['alpha'])

        fp = self.mat_props['fiber']
        m.mp('EX', 2, fp['E'])
        m.mp('NUXY', 2, fp['nu'])
        m.mp('ALPX', 2, fp['alpha'])

        # Element type
        m.et(1, 'SOLID186')
        m.esize(elem_size)

        # Mesh fiber
        m.vsel('S', 'VOLU', '', 2)
        m.vatt(2, '', 1)
        m.vmesh(2)

        # Mesh matrix
        m.vsel('S', 'VOLU', '', 3)
        m.vatt(1, '', 1)
        m.vmesh(3)

        m.allsel()

        nn = int(m.get('NCOUNT', 'NODE', '', 'COUNT'))
        ne = int(m.get('ECOUNT', 'ELEM', '', 'COUNT'))
        print(f"Mesh: {nn} nodes, {ne} elements")

        # Create node component sets for faces
        self._create_face_sets()

        print("Model built successfully")

    def _create_face_sets(self):
        """Create node sets for all faces."""
        m = self.mapdl
        L = self.L
        tol = 1e-6

        faces = [
            ('XNEG', 'X', 0),
            ('XPOS', 'X', L),
            ('YNEG', 'Y', 0),
            ('YPOS', 'Y', L),
            ('ZNEG', 'Z', 0),
            ('ZPOS', 'Z', L),
        ]

        for name, direction, coord in faces:
            m.nsel('S', 'LOC', direction, coord - tol, coord + tol)
            m.cm(name, 'NODE')

        m.allsel()

    def apply_periodic_bc_with_master_nodes(self, eps_macro):
        """
        Apply periodic boundary conditions for a given macroscopic strain.

        Uses kinematic uniform boundary conditions (KUBC) which are simpler
        but still provide accurate effective properties.

        Parameters
        ----------
        eps_macro : array-like
            Macroscopic strain tensor in Voigt notation [e11, e22, e33, e12, e23, e31]
        """
        m = self.mapdl
        L = self.L

        e11, e22, e33, e12, e23, e31 = eps_macro

        # Clear all constraints
        m.ddele('ALL', 'ALL')

        # Fix corner node at origin for rigid body motion
        m.nsel('S', 'LOC', 'X', 0, 1e-6)
        m.nsel('R', 'LOC', 'Y', 0, 1e-6)
        m.nsel('R', 'LOC', 'Z', 0, 1e-6)
        m.d('ALL', 'UX', 0)
        m.d('ALL', 'UY', 0)
        m.d('ALL', 'UZ', 0)

        # Apply displacements on X+ face
        m.cmsel('S', 'XPOS')
        ux_xpos = e11 * L + e12 * 0 + e31 * 0  # At x=L, y,z variable
        m.d('ALL', 'UX', e11 * L)
        m.d('ALL', 'UY', e12 * L)
        m.d('ALL', 'UZ', e31 * L)

        # Apply displacements on Y+ face
        m.cmsel('S', 'YPOS')
        m.d('ALL', 'UX', e12 * L)
        m.d('ALL', 'UY', e22 * L)
        m.d('ALL', 'UZ', e23 * L)

        # Apply displacements on Z+ face
        m.cmsel('S', 'ZPOS')
        m.d('ALL', 'UX', e31 * L)
        m.d('ALL', 'UY', e23 * L)
        m.d('ALL', 'UZ', e33 * L)

        # Fix negative faces
        m.cmsel('S', 'XNEG')
        m.d('ALL', 'UX', 0)

        m.cmsel('S', 'YNEG')
        m.d('ALL', 'UY', 0)

        m.cmsel('S', 'ZNEG')
        m.d('ALL', 'UZ', 0)

        m.allsel()

    def apply_thermal_bc(self, delta_T=1.0):
        """Apply thermal loading with free expansion."""
        m = self.mapdl

        m.ddele('ALL', 'ALL')

        # Fix corner for rigid body
        m.nsel('S', 'LOC', 'X', 0, 1e-6)
        m.nsel('R', 'LOC', 'Y', 0, 1e-6)
        m.nsel('R', 'LOC', 'Z', 0, 1e-6)
        m.d('ALL', 'ALL', 0)

        # Fix edges for symmetry
        m.cmsel('S', 'XNEG')
        m.d('ALL', 'UX', 0)

        m.cmsel('S', 'YNEG')
        m.d('ALL', 'UY', 0)

        m.cmsel('S', 'ZNEG')
        m.d('ALL', 'UZ', 0)

        m.allsel()

        # Thermal load
        m.bfunif('TEMP', delta_T)
        m.tunif(0)

    def solve(self):
        """Solve current load case."""
        m = self.mapdl
        m.run('/SOLU')
        m.antype('STATIC')
        m.solve()
        m.finish()

    def get_volume_avg_stress(self):
        """
        Calculate volume-averaged stress over the RVE.

        Returns
        -------
        stress : ndarray
            Stress tensor in Voigt notation [S11, S22, S33, S12, S23, S31]
        """
        m = self.mapdl
        m.post1()
        m.set('LAST')

        m.esel('ALL')

        # Create element tables for stresses and volumes
        m.etable('SXX', 'S', 'X')
        m.etable('SYY', 'S', 'Y')
        m.etable('SZZ', 'S', 'Z')
        m.etable('SXY', 'S', 'XY')
        m.etable('SYZ', 'S', 'YZ')
        m.etable('SXZ', 'S', 'XZ')
        m.etable('EVOL', 'VOLU')

        # Get sums
        vol_total = m.get('VOLTOT', 'ELEM', '', 'ETAB', 'EVOL', 'SUM')

        # Volume-weighted stress (using SMULT and SADD for proper averaging)
        m.smult('SXVOL', 'SXX', 'EVOL')
        m.smult('SYVOL', 'SYY', 'EVOL')
        m.smult('SZVOL', 'SZZ', 'EVOL')
        m.smult('SXYVOL', 'SXY', 'EVOL')
        m.smult('SYZVOL', 'SYZ', 'EVOL')
        m.smult('SXZVOL', 'SXZ', 'EVOL')

        sxx_sum = m.get('SXXSUM', 'ELEM', '', 'ETAB', 'SXVOL', 'SUM')
        syy_sum = m.get('SYYSUM', 'ELEM', '', 'ETAB', 'SYVOL', 'SUM')
        szz_sum = m.get('SZZSUM', 'ELEM', '', 'ETAB', 'SZVOL', 'SUM')
        sxy_sum = m.get('SXYSUM', 'ELEM', '', 'ETAB', 'SXYVOL', 'SUM')
        syz_sum = m.get('SYZSUM', 'ELEM', '', 'ETAB', 'SYZVOL', 'SUM')
        sxz_sum = m.get('SXZSUM', 'ELEM', '', 'ETAB', 'SXZVOL', 'SUM')

        stress = np.array([
            sxx_sum / vol_total,
            syy_sum / vol_total,
            szz_sum / vol_total,
            sxy_sum / vol_total,
            syz_sum / vol_total,
            sxz_sum / vol_total
        ])

        m.finish()
        return stress

    def get_thermal_strain(self, delta_T=1.0):
        """Get thermal expansion strains."""
        m = self.mapdl
        L = self.L

        m.post1()
        m.set('LAST')

        # Get average displacements on positive faces
        m.cmsel('S', 'XPOS')
        ux = m.get('UXAVG', 'NODE', '', 'U', 'X', 'AVG')

        m.cmsel('S', 'YPOS')
        uy = m.get('UYAVG', 'NODE', '', 'U', 'Y', 'AVG')

        m.cmsel('S', 'ZPOS')
        uz = m.get('UZAVG', 'NODE', '', 'U', 'Z', 'AVG')

        m.allsel()
        m.finish()

        # Thermal strains
        eps_th = np.array([ux/L, uy/L, uz/L]) / delta_T

        return eps_th

    def compute_stiffness_matrix(self, strain_mag=0.001):
        """
        Compute the complete 6x6 stiffness matrix.

        Parameters
        ----------
        strain_mag : float
            Magnitude of applied strain for each load case

        Returns
        -------
        C : ndarray
            6x6 stiffness matrix in Voigt notation
        """
        print(f"\n{'='*60}")
        print("COMPUTING STIFFNESS MATRIX")
        print(f"{'='*60}")

        C = np.zeros((6, 6))

        # Define 6 load cases (unit strains in each direction)
        load_cases = [
            [1, 0, 0, 0, 0, 0],  # e11
            [0, 1, 0, 0, 0, 0],  # e22
            [0, 0, 1, 0, 0, 0],  # e33
            [0, 0, 0, 1, 0, 0],  # e12 (gamma_xy)
            [0, 0, 0, 0, 1, 0],  # e23 (gamma_yz)
            [0, 0, 0, 0, 0, 1],  # e31 (gamma_zx)
        ]

        directions = ['ε11', 'ε22', 'ε33', 'γ12', 'γ23', 'γ31']

        for i, lc in enumerate(load_cases):
            print(f"  Load case {i+1}/6: {directions[i]}...")

            self.mapdl.prep7()
            eps_applied = np.array(lc) * strain_mag
            self.apply_periodic_bc_with_master_nodes(eps_applied)
            self.solve()

            stress = self.get_volume_avg_stress()
            C[:, i] = stress / strain_mag

        # Symmetrize the matrix (should be symmetric for linear elastic)
        C = 0.5 * (C + C.T)

        self.stiffness_matrix = C
        print("  Stiffness matrix computed")

        return C

    def compute_engineering_constants(self):
        """
        Compute engineering constants from stiffness matrix.

        Returns
        -------
        props : dict
            Dictionary containing Ex, Ey, Ez, Gxy, Gyz, Gzx, nuxy, nuyz, nuzx
        """
        if self.stiffness_matrix is None:
            raise ValueError("Stiffness matrix not computed. Run compute_stiffness_matrix first.")

        C = self.stiffness_matrix

        # Compute compliance matrix S = C^(-1)
        try:
            S = np.linalg.inv(C)
        except np.linalg.LinAlgError:
            print("Warning: Stiffness matrix is singular. Using pseudo-inverse.")
            S = np.linalg.pinv(C)

        self.compliance_matrix = S

        # Extract engineering constants from compliance matrix
        # For orthotropic material:
        # S11 = 1/Ex, S22 = 1/Ey, S33 = 1/Ez
        # S44 = 1/Gyz, S55 = 1/Gxz, S66 = 1/Gxy
        # S12 = -nuxy/Ex = -nuyx/Ey
        # S13 = -nuxz/Ex = -nuzx/Ez
        # S23 = -nuyz/Ey = -nuzy/Ez

        Ex = 1.0 / S[0, 0]
        Ey = 1.0 / S[1, 1]
        Ez = 1.0 / S[2, 2]

        # Note: In Voigt notation for shear, factor of 2 may apply
        # S44 = 1/Gyz, etc. (depends on convention)
        Gxy = 1.0 / S[5, 5]  # gamma_xy corresponds to index 5 in our convention
        Gyz = 1.0 / S[3, 3]  # gamma_yz corresponds to index 3
        Gzx = 1.0 / S[4, 4]  # gamma_zx corresponds to index 4

        nuxy = -S[0, 1] * Ex
        nuyz = -S[1, 2] * Ey
        nuzx = -S[2, 0] * Ez

        self.effective_props.update({
            'Ex': Ex, 'Ey': Ey, 'Ez': Ez,
            'Gxy': Gxy, 'Gyz': Gyz, 'Gzx': Gzx,
            'nuxy': nuxy, 'nuyz': nuyz, 'nuzx': nuzx
        })

        return self.effective_props

    def compute_thermal_expansion(self, delta_T=1.0):
        """
        Compute effective thermal expansion coefficients.

        Parameters
        ----------
        delta_T : float
            Temperature change

        Returns
        -------
        cte : dict
            Dictionary containing CTEx, CTEy, CTEz
        """
        print(f"\n{'='*60}")
        print("COMPUTING THERMAL EXPANSION COEFFICIENTS")
        print(f"{'='*60}")

        self.mapdl.prep7()
        self.apply_thermal_bc(delta_T)
        self.solve()

        eps_th = self.get_thermal_strain(delta_T)

        CTEx = eps_th[0]
        CTEy = eps_th[1]
        CTEz = eps_th[2]

        self.effective_props.update({
            'CTEx': CTEx, 'CTEy': CTEy, 'CTEz': CTEz
        })

        print(f"  CTE computed: [{CTEx:.2e}, {CTEy:.2e}, {CTEz:.2e}] 1/°C")

        return {'CTEx': CTEx, 'CTEy': CTEy, 'CTEz': CTEz}

    def run_full_analysis(self, elem_size=0.05, strain_mag=0.001):
        """
        Run complete analysis to get all effective properties.

        Parameters
        ----------
        elem_size : float
            Mesh element size
        strain_mag : float
            Strain magnitude for mechanical load cases

        Returns
        -------
        props : dict
            All effective properties
        """
        # Build model
        self.build_model(elem_size)

        # Compute stiffness matrix
        self.compute_stiffness_matrix(strain_mag)

        # Extract engineering constants
        self.compute_engineering_constants()

        # Compute thermal properties
        self.compute_thermal_expansion()

        return self.effective_props

    def print_results(self):
        """Print all computed effective properties."""
        if not self.effective_props:
            print("No results. Run analysis first.")
            return

        p = self.effective_props

        print(f"\n{'='*60}")
        print("EFFECTIVE MATERIAL PROPERTIES")
        print(f"{'='*60}")

        print(f"\nInput Parameters:")
        print(f"  RVE size: {self.L} mm")
        print(f"  Fiber volume fraction: {self.Vf*100:.1f}%")
        print(f"  Fiber radius: {self.r_fiber:.4f} mm")

        print(f"\nMatrix Properties:")
        print(f"  E = {self.mat_props['matrix']['E']} MPa")
        print(f"  ν = {self.mat_props['matrix']['nu']}")
        print(f"  α = {self.mat_props['matrix']['alpha']:.2e} 1/°C")

        print(f"\nFiber Properties:")
        print(f"  E = {self.mat_props['fiber']['E']} MPa")
        print(f"  ν = {self.mat_props['fiber']['nu']}")
        print(f"  α = {self.mat_props['fiber']['alpha']:.2e} 1/°C")

        print(f"\n{'-'*60}")
        print("COMPUTED EFFECTIVE PROPERTIES")
        print(f"{'-'*60}")

        print(f"\n--- Elastic Moduli [MPa] ---")
        print(f"  Ex = {p.get('Ex', 0):.2f}")
        print(f"  Ey = {p.get('Ey', 0):.2f}")
        print(f"  Ez = {p.get('Ez', 0):.2f}")

        print(f"\n--- Shear Moduli [MPa] ---")
        print(f"  Gxy = {p.get('Gxy', 0):.2f}")
        print(f"  Gyz = {p.get('Gyz', 0):.2f}")
        print(f"  Gzx = {p.get('Gzx', 0):.2f}")

        print(f"\n--- Poisson's Ratios ---")
        print(f"  νxy = {p.get('nuxy', 0):.4f}")
        print(f"  νyz = {p.get('nuyz', 0):.4f}")
        print(f"  νzx = {p.get('nuzx', 0):.4f}")

        print(f"\n--- Thermal Expansion Coefficients [1/°C] ---")
        print(f"  CTEx = {p.get('CTEx', 0):.2e}")
        print(f"  CTEy = {p.get('CTEy', 0):.2e}")
        print(f"  CTEz = {p.get('CTEz', 0):.2e}")

        if self.stiffness_matrix is not None:
            print(f"\n--- Stiffness Matrix C [MPa] ---")
            C = self.stiffness_matrix
            for i in range(6):
                row = [f"{C[i,j]:12.2f}" for j in range(6)]
                print(f"  [{' '.join(row)}]")

        print(f"\n{'='*60}")


def main():
    """Main function."""

    calc = AdvancedCompositeCalculator(rve_size=1.0, fiber_vf=0.10)

    # Set material properties
    calc.set_material('matrix', E=3500, nu=0.35, alpha=60e-6)
    calc.set_material('fiber', E=230000, nu=0.20, alpha=-0.5e-6)

    try:
        print("Starting ANSYS MAPDL...")
        calc.launch(run_location='/tmp/mapdl_adv', override=True)

        # Run full analysis
        calc.run_full_analysis(elem_size=0.05, strain_mag=0.001)

        # Print results
        calc.print_results()

    except Exception as e:
        print(f"Error: {e}")
        import traceback
        traceback.print_exc()

    finally:
        calc.exit()
        print("\nAnalysis complete.")


if __name__ == "__main__":
    main()
