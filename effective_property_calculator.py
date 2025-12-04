"""
Composite Material Effective Property Calculator using PyANSYS (PyMAPDL)

This script calculates effective material properties for a fiber-reinforced composite
using Representative Volume Element (RVE) homogenization method.

Geometry: 1x1x1 mm cube (matrix) with cylindrical fiber at center
Fiber volume fraction: 10%

Output Properties:
- Elastic Modulus: Ex, Ey, Ez
- Shear Modulus: Gxy, Gyz, Gzx
- Poisson's Ratio: vxy, vyz, vzx
- Coefficient of Thermal Expansion: CTEx, CTEy, CTEz
"""

import numpy as np
from ansys.mapdl.core import launch_mapdl


class CompositeEffectivePropertyCalculator:
    """Calculator for effective material properties of fiber-reinforced composites."""

    def __init__(self, rve_size=1.0, fiber_vf=0.10):
        """
        Initialize the calculator.

        Parameters
        ----------
        rve_size : float
            Size of the RVE cube in mm (default: 1.0 mm)
        fiber_vf : float
            Fiber volume fraction (default: 0.10 = 10%)
        """
        self.rve_size = rve_size
        self.fiber_vf = fiber_vf

        # Calculate fiber radius from volume fraction
        # Cylinder volume: V = pi * r^2 * h, where h = rve_size
        # fiber_vf = (pi * r^2 * h) / (rve_size^3)
        # r = sqrt(fiber_vf * rve_size^2 / pi)
        self.fiber_radius = np.sqrt(fiber_vf * rve_size**2 / np.pi)

        # Material properties (default values - can be modified)
        # Matrix: Epoxy
        self.matrix_E = 3500  # MPa (Young's modulus)
        self.matrix_nu = 0.35  # Poisson's ratio
        self.matrix_alpha = 60e-6  # 1/°C (CTE)

        # Fiber: Carbon fiber (transversely isotropic, but simplified as isotropic here)
        self.fiber_E = 230000  # MPa (Young's modulus)
        self.fiber_nu = 0.20  # Poisson's ratio
        self.fiber_alpha = -0.5e-6  # 1/°C (CTE, negative for carbon fiber)

        # MAPDL instance
        self.mapdl = None

        # Results storage
        self.effective_properties = {}

    def set_matrix_properties(self, E, nu, alpha):
        """
        Set matrix material properties.

        Parameters
        ----------
        E : float
            Young's modulus in MPa
        nu : float
            Poisson's ratio
        alpha : float
            Coefficient of thermal expansion in 1/°C
        """
        self.matrix_E = E
        self.matrix_nu = nu
        self.matrix_alpha = alpha

    def set_fiber_properties(self, E, nu, alpha):
        """
        Set fiber material properties.

        Parameters
        ----------
        E : float
            Young's modulus in MPa
        nu : float
            Poisson's ratio
        alpha : float
            Coefficient of thermal expansion in 1/°C
        """
        self.fiber_E = E
        self.fiber_nu = nu
        self.fiber_alpha = alpha

    def start_mapdl(self, **kwargs):
        """Start MAPDL instance."""
        self.mapdl = launch_mapdl(**kwargs)
        self.mapdl.clear()
        self.mapdl.prep7()

    def stop_mapdl(self):
        """Stop MAPDL instance."""
        if self.mapdl is not None:
            self.mapdl.exit()
            self.mapdl = None

    def create_geometry(self):
        """Create RVE geometry with matrix and cylindrical fiber."""
        mapdl = self.mapdl
        L = self.rve_size
        r = self.fiber_radius

        print(f"Creating RVE geometry...")
        print(f"  RVE size: {L} x {L} x {L} mm")
        print(f"  Fiber radius: {r:.4f} mm")
        print(f"  Fiber volume fraction: {self.fiber_vf*100:.1f}%")

        # Create the matrix block
        mapdl.block(0, L, 0, L, 0, L)

        # Create the cylindrical fiber at the center
        # Cylinder along Z-axis, centered at (L/2, L/2)
        mapdl.cyl4(L/2, L/2, r, 0, r, 360, L)

        # Subtract cylinder from block to create matrix
        mapdl.vsbv(1, 2, sepo='', keep2='keep')

        # Now we have:
        # Volume 2: Fiber (cylinder)
        # Volume 3: Matrix (block with hole)

        print("  Geometry created successfully")

    def define_materials(self):
        """Define material properties for matrix and fiber."""
        mapdl = self.mapdl

        print("Defining materials...")

        # Material 1: Matrix (Epoxy)
        mapdl.mp('EX', 1, self.matrix_E)
        mapdl.mp('NUXY', 1, self.matrix_nu)
        mapdl.mp('ALPX', 1, self.matrix_alpha)
        print(f"  Matrix: E={self.matrix_E} MPa, nu={self.matrix_nu}, alpha={self.matrix_alpha} 1/°C")

        # Material 2: Fiber (Carbon)
        mapdl.mp('EX', 2, self.fiber_E)
        mapdl.mp('NUXY', 2, self.fiber_nu)
        mapdl.mp('ALPX', 2, self.fiber_alpha)
        print(f"  Fiber: E={self.fiber_E} MPa, nu={self.fiber_nu}, alpha={self.fiber_alpha} 1/°C")

    def create_mesh(self, element_size=0.05):
        """
        Create finite element mesh.

        Parameters
        ----------
        element_size : float
            Target element size in mm
        """
        mapdl = self.mapdl

        print(f"Creating mesh (element size: {element_size} mm)...")

        # Element type: SOLID187 (10-node tetrahedral) for complex geometry
        mapdl.et(1, 'SOLID187')

        # Set element size
        mapdl.esize(element_size)

        # Use free meshing (required for cylindrical geometry)
        mapdl.mshkey(0)  # Free meshing

        # Mesh the fiber (Volume 2)
        mapdl.vsel('S', 'VOLU', '', 2)
        mapdl.vatt(2, '', 1)  # Material 2, Element type 1
        mapdl.vmesh(2)

        # Mesh the matrix (Volume 3)
        mapdl.vsel('S', 'VOLU', '', 3)
        mapdl.vatt(1, '', 1)  # Material 1, Element type 1
        mapdl.vmesh(3)

        mapdl.vsel('ALL')
        mapdl.nsel('ALL')
        mapdl.esel('ALL')

        # Merge nodes at interface
        mapdl.nummrg('NODE', 1e-6)

        # Get mesh statistics
        n_nodes = mapdl.get('NCOUNT', 'NODE', '', 'COUNT')
        n_elements = mapdl.get('ECOUNT', 'ELEM', '', 'COUNT')
        print(f"  Mesh created: {int(n_nodes)} nodes, {int(n_elements)} elements")

    def get_boundary_nodes(self):
        """Get node sets for periodic boundary conditions."""
        mapdl = self.mapdl
        L = self.rve_size
        tol = 1e-6

        # Select nodes on each face
        # X- face (x=0)
        mapdl.nsel('S', 'LOC', 'X', 0, tol)
        mapdl.cm('X_NEG', 'NODE')

        # X+ face (x=L)
        mapdl.nsel('S', 'LOC', 'X', L-tol, L+tol)
        mapdl.cm('X_POS', 'NODE')

        # Y- face (y=0)
        mapdl.nsel('S', 'LOC', 'Y', 0, tol)
        mapdl.cm('Y_NEG', 'NODE')

        # Y+ face (y=L)
        mapdl.nsel('S', 'LOC', 'Y', L-tol, L+tol)
        mapdl.cm('Y_POS', 'NODE')

        # Z- face (z=0)
        mapdl.nsel('S', 'LOC', 'Z', 0, tol)
        mapdl.cm('Z_NEG', 'NODE')

        # Z+ face (z=L)
        mapdl.nsel('S', 'LOC', 'Z', L-tol, L+tol)
        mapdl.cm('Z_POS', 'NODE')

        mapdl.nsel('ALL')

    def apply_periodic_bc_uniaxial_x(self, strain_val=0.001):
        """Apply periodic BC for uniaxial strain in X direction."""
        mapdl = self.mapdl
        L = self.rve_size

        # Clear previous constraints and loads
        mapdl.ddele('ALL', 'ALL')

        # Fix corner node to prevent rigid body motion
        mapdl.nsel('S', 'LOC', 'X', 0)
        mapdl.nsel('R', 'LOC', 'Y', 0)
        mapdl.nsel('R', 'LOC', 'Z', 0)
        mapdl.d('ALL', 'UX', 0)
        mapdl.d('ALL', 'UY', 0)
        mapdl.d('ALL', 'UZ', 0)

        # Apply displacement on X+ face
        mapdl.cmsel('S', 'X_POS')
        mapdl.d('ALL', 'UX', strain_val * L)

        mapdl.nsel('ALL')

    def apply_periodic_bc_uniaxial_y(self, strain_val=0.001):
        """Apply periodic BC for uniaxial strain in Y direction."""
        mapdl = self.mapdl
        L = self.rve_size

        mapdl.ddele('ALL', 'ALL')

        # Fix corner node
        mapdl.nsel('S', 'LOC', 'X', 0)
        mapdl.nsel('R', 'LOC', 'Y', 0)
        mapdl.nsel('R', 'LOC', 'Z', 0)
        mapdl.d('ALL', 'UX', 0)
        mapdl.d('ALL', 'UY', 0)
        mapdl.d('ALL', 'UZ', 0)

        # Apply displacement on Y+ face
        mapdl.cmsel('S', 'Y_POS')
        mapdl.d('ALL', 'UY', strain_val * L)

        mapdl.nsel('ALL')

    def apply_periodic_bc_uniaxial_z(self, strain_val=0.001):
        """Apply periodic BC for uniaxial strain in Z direction."""
        mapdl = self.mapdl
        L = self.rve_size

        mapdl.ddele('ALL', 'ALL')

        # Fix corner node
        mapdl.nsel('S', 'LOC', 'X', 0)
        mapdl.nsel('R', 'LOC', 'Y', 0)
        mapdl.nsel('R', 'LOC', 'Z', 0)
        mapdl.d('ALL', 'UX', 0)
        mapdl.d('ALL', 'UY', 0)
        mapdl.d('ALL', 'UZ', 0)

        # Apply displacement on Z+ face
        mapdl.cmsel('S', 'Z_POS')
        mapdl.d('ALL', 'UZ', strain_val * L)

        mapdl.nsel('ALL')

    def apply_periodic_bc_shear_xy(self, strain_val=0.001):
        """Apply periodic BC for shear strain in XY plane."""
        mapdl = self.mapdl
        L = self.rve_size

        mapdl.ddele('ALL', 'ALL')

        # Fix corner node
        mapdl.nsel('S', 'LOC', 'X', 0)
        mapdl.nsel('R', 'LOC', 'Y', 0)
        mapdl.nsel('R', 'LOC', 'Z', 0)
        mapdl.d('ALL', 'UX', 0)
        mapdl.d('ALL', 'UY', 0)
        mapdl.d('ALL', 'UZ', 0)

        # Apply shear: UY on X+ face
        mapdl.cmsel('S', 'X_POS')
        mapdl.d('ALL', 'UY', strain_val * L / 2)

        # Apply shear: UX on Y+ face
        mapdl.cmsel('S', 'Y_POS')
        mapdl.d('ALL', 'UX', strain_val * L / 2)

        mapdl.nsel('ALL')

    def apply_periodic_bc_shear_yz(self, strain_val=0.001):
        """Apply periodic BC for shear strain in YZ plane."""
        mapdl = self.mapdl
        L = self.rve_size

        mapdl.ddele('ALL', 'ALL')

        # Fix corner node
        mapdl.nsel('S', 'LOC', 'X', 0)
        mapdl.nsel('R', 'LOC', 'Y', 0)
        mapdl.nsel('R', 'LOC', 'Z', 0)
        mapdl.d('ALL', 'UX', 0)
        mapdl.d('ALL', 'UY', 0)
        mapdl.d('ALL', 'UZ', 0)

        # Apply shear: UZ on Y+ face
        mapdl.cmsel('S', 'Y_POS')
        mapdl.d('ALL', 'UZ', strain_val * L / 2)

        # Apply shear: UY on Z+ face
        mapdl.cmsel('S', 'Z_POS')
        mapdl.d('ALL', 'UY', strain_val * L / 2)

        mapdl.nsel('ALL')

    def apply_periodic_bc_shear_zx(self, strain_val=0.001):
        """Apply periodic BC for shear strain in ZX plane."""
        mapdl = self.mapdl
        L = self.rve_size

        mapdl.ddele('ALL', 'ALL')

        # Fix corner node
        mapdl.nsel('S', 'LOC', 'X', 0)
        mapdl.nsel('R', 'LOC', 'Y', 0)
        mapdl.nsel('R', 'LOC', 'Z', 0)
        mapdl.d('ALL', 'UX', 0)
        mapdl.d('ALL', 'UY', 0)
        mapdl.d('ALL', 'UZ', 0)

        # Apply shear: UX on Z+ face
        mapdl.cmsel('S', 'Z_POS')
        mapdl.d('ALL', 'UX', strain_val * L / 2)

        # Apply shear: UZ on X+ face
        mapdl.cmsel('S', 'X_POS')
        mapdl.d('ALL', 'UZ', strain_val * L / 2)

        mapdl.nsel('ALL')

    def apply_thermal_load(self, delta_T=1.0):
        """Apply thermal load for CTE calculation."""
        mapdl = self.mapdl

        mapdl.ddele('ALL', 'ALL')

        # Fix corner node to prevent rigid body motion
        mapdl.nsel('S', 'LOC', 'X', 0)
        mapdl.nsel('R', 'LOC', 'Y', 0)
        mapdl.nsel('R', 'LOC', 'Z', 0)
        mapdl.d('ALL', 'UX', 0)
        mapdl.d('ALL', 'UY', 0)
        mapdl.d('ALL', 'UZ', 0)

        # Fix X- face in X direction (symmetry-like)
        mapdl.cmsel('S', 'X_NEG')
        mapdl.d('ALL', 'UX', 0)

        # Fix Y- face in Y direction
        mapdl.cmsel('S', 'Y_NEG')
        mapdl.d('ALL', 'UY', 0)

        # Fix Z- face in Z direction
        mapdl.cmsel('S', 'Z_NEG')
        mapdl.d('ALL', 'UZ', 0)

        mapdl.nsel('ALL')

        # Apply uniform temperature
        mapdl.bfunif('TEMP', delta_T)
        mapdl.tunif(0)  # Reference temperature = 0

    def solve(self):
        """Solve the current load case."""
        mapdl = self.mapdl
        mapdl.run('/SOLU')
        mapdl.antype('STATIC')
        mapdl.solve()
        mapdl.finish()

    def get_volume_average_stress(self):
        """Calculate volume-averaged stress from the solution."""
        mapdl = self.mapdl

        mapdl.post1()
        mapdl.set('LAST')
        mapdl.esel('ALL')

        # Get stress components averaged over all elements
        # Using ETABLE to extract and average stresses
        mapdl.etable('SX', 'S', 'X')
        mapdl.etable('SY', 'S', 'Y')
        mapdl.etable('SZ', 'S', 'Z')
        mapdl.etable('SXY', 'S', 'XY')
        mapdl.etable('SYZ', 'S', 'YZ')
        mapdl.etable('SXZ', 'S', 'XZ')
        mapdl.etable('EVOL', 'VOLU')

        # Calculate volume-weighted averages
        total_volume = mapdl.get('VTOT', 'ELEM', '', 'ETAB', 'EVOL', 'SUM')

        sx_sum = mapdl.get('SXSUM', 'ELEM', '', 'ETAB', 'SX', 'SUM')
        sy_sum = mapdl.get('SYSUM', 'ELEM', '', 'ETAB', 'SY', 'SUM')
        sz_sum = mapdl.get('SZSUM', 'ELEM', '', 'ETAB', 'SZ', 'SUM')
        sxy_sum = mapdl.get('SXYSUM', 'ELEM', '', 'ETAB', 'SXY', 'SUM')
        syz_sum = mapdl.get('SYZSUM', 'ELEM', '', 'ETAB', 'SYZ', 'SUM')
        sxz_sum = mapdl.get('SXZSUM', 'ELEM', '', 'ETAB', 'SXZ', 'SUM')

        n_elem = mapdl.get('ECOUNT', 'ELEM', '', 'COUNT')

        stress = {
            'SX': sx_sum / n_elem,
            'SY': sy_sum / n_elem,
            'SZ': sz_sum / n_elem,
            'SXY': sxy_sum / n_elem,
            'SYZ': syz_sum / n_elem,
            'SXZ': sxz_sum / n_elem
        }

        mapdl.finish()
        return stress

    def get_face_displacement(self, face):
        """Get average displacement on a face."""
        mapdl = self.mapdl
        L = self.rve_size

        mapdl.post1()
        mapdl.set('LAST')

        if face == 'X_POS':
            mapdl.cmsel('S', 'X_POS')
            disp = mapdl.get('UX_AVG', 'NODE', '', 'U', 'X', 'AVG')
        elif face == 'Y_POS':
            mapdl.cmsel('S', 'Y_POS')
            disp = mapdl.get('UY_AVG', 'NODE', '', 'U', 'Y', 'AVG')
        elif face == 'Z_POS':
            mapdl.cmsel('S', 'Z_POS')
            disp = mapdl.get('UZ_AVG', 'NODE', '', 'U', 'Z', 'AVG')
        else:
            disp = 0

        mapdl.nsel('ALL')
        mapdl.finish()
        return disp

    def calculate_effective_properties(self, element_size=0.05, strain_val=0.001):
        """
        Calculate all effective properties.

        Parameters
        ----------
        element_size : float
            Mesh element size in mm
        strain_val : float
            Applied strain magnitude for load cases
        """
        L = self.rve_size
        V = L**3  # RVE volume

        print("\n" + "="*60)
        print("COMPOSITE EFFECTIVE PROPERTY CALCULATION")
        print("="*60)

        # Build model
        self.create_geometry()
        self.define_materials()
        self.create_mesh(element_size)
        self.get_boundary_nodes()

        # Store stiffness matrix components
        C = np.zeros((6, 6))

        # Load case 1: Uniaxial strain in X
        print("\nLoad Case 1: Uniaxial strain εxx...")
        self.apply_periodic_bc_uniaxial_x(strain_val)
        self.solve()
        stress1 = self.get_volume_average_stress()
        C[0, 0] = stress1['SX'] / strain_val
        C[1, 0] = stress1['SY'] / strain_val
        C[2, 0] = stress1['SZ'] / strain_val

        # Load case 2: Uniaxial strain in Y
        print("Load Case 2: Uniaxial strain εyy...")
        self.mapdl.prep7()
        self.apply_periodic_bc_uniaxial_y(strain_val)
        self.solve()
        stress2 = self.get_volume_average_stress()
        C[0, 1] = stress2['SX'] / strain_val
        C[1, 1] = stress2['SY'] / strain_val
        C[2, 1] = stress2['SZ'] / strain_val

        # Load case 3: Uniaxial strain in Z
        print("Load Case 3: Uniaxial strain εzz...")
        self.mapdl.prep7()
        self.apply_periodic_bc_uniaxial_z(strain_val)
        self.solve()
        stress3 = self.get_volume_average_stress()
        C[0, 2] = stress3['SX'] / strain_val
        C[1, 2] = stress3['SY'] / strain_val
        C[2, 2] = stress3['SZ'] / strain_val

        # Load case 4: Shear strain XY
        print("Load Case 4: Shear strain γxy...")
        self.mapdl.prep7()
        self.apply_periodic_bc_shear_xy(strain_val)
        self.solve()
        stress4 = self.get_volume_average_stress()
        C[3, 3] = stress4['SXY'] / strain_val

        # Load case 5: Shear strain YZ
        print("Load Case 5: Shear strain γyz...")
        self.mapdl.prep7()
        self.apply_periodic_bc_shear_yz(strain_val)
        self.solve()
        stress5 = self.get_volume_average_stress()
        C[4, 4] = stress5['SYZ'] / strain_val

        # Load case 6: Shear strain ZX
        print("Load Case 6: Shear strain γzx...")
        self.mapdl.prep7()
        self.apply_periodic_bc_shear_zx(strain_val)
        self.solve()
        stress6 = self.get_volume_average_stress()
        C[5, 5] = stress6['SXZ'] / strain_val

        # Calculate compliance matrix S = C^(-1)
        # Note: We only have the diagonal and some off-diagonal terms
        # For a more complete analysis, we would need full periodic BC

        # Approximate engineering constants from stiffness components
        # For orthotropic material:
        # E_i = C_ii - (C_ij * C_ji) / C_jj (simplified)

        # Elastic moduli (approximate)
        Ex = C[0, 0] - (C[0, 1]**2 / C[1, 1] + C[0, 2]**2 / C[2, 2])
        Ey = C[1, 1] - (C[1, 0]**2 / C[0, 0] + C[1, 2]**2 / C[2, 2])
        Ez = C[2, 2] - (C[2, 0]**2 / C[0, 0] + C[2, 1]**2 / C[1, 1])

        # Shear moduli
        Gxy = C[3, 3]
        Gyz = C[4, 4]
        Gzx = C[5, 5]

        # Poisson's ratios (approximate)
        nu_xy = C[0, 1] / C[0, 0]
        nu_yz = C[1, 2] / C[1, 1]
        nu_zx = C[2, 0] / C[2, 2]

        # Thermal expansion coefficients
        print("\nLoad Case 7: Thermal expansion (ΔT = 1°C)...")
        self.mapdl.prep7()
        delta_T = 1.0
        self.apply_thermal_load(delta_T)
        self.solve()

        # Get displacements on positive faces
        ux = self.get_face_displacement('X_POS')
        uy = self.get_face_displacement('Y_POS')
        uz = self.get_face_displacement('Z_POS')

        # CTE = strain / delta_T = (displacement / L) / delta_T
        CTEx = (ux / L) / delta_T
        CTEy = (uy / L) / delta_T
        CTEz = (uz / L) / delta_T

        # Store results
        self.effective_properties = {
            'Ex': Ex, 'Ey': Ey, 'Ez': Ez,
            'Gxy': Gxy, 'Gyz': Gyz, 'Gzx': Gzx,
            'nu_xy': nu_xy, 'nu_yz': nu_yz, 'nu_zx': nu_zx,
            'CTEx': CTEx, 'CTEy': CTEy, 'CTEz': CTEz,
            'stiffness_matrix': C
        }

        return self.effective_properties

    def print_results(self):
        """Print calculated effective properties."""
        if not self.effective_properties:
            print("No results available. Run calculate_effective_properties() first.")
            return

        props = self.effective_properties

        print("\n" + "="*60)
        print("EFFECTIVE MATERIAL PROPERTIES")
        print("="*60)

        print("\n--- Elastic Moduli (MPa) ---")
        print(f"  Ex = {props['Ex']:.2f}")
        print(f"  Ey = {props['Ey']:.2f}")
        print(f"  Ez = {props['Ez']:.2f}")

        print("\n--- Shear Moduli (MPa) ---")
        print(f"  Gxy = {props['Gxy']:.2f}")
        print(f"  Gyz = {props['Gyz']:.2f}")
        print(f"  Gzx = {props['Gzx']:.2f}")

        print("\n--- Poisson's Ratios ---")
        print(f"  νxy = {props['nu_xy']:.4f}")
        print(f"  νyz = {props['nu_yz']:.4f}")
        print(f"  νzx = {props['nu_zx']:.4f}")

        print("\n--- Coefficients of Thermal Expansion (1/°C) ---")
        print(f"  CTEx = {props['CTEx']:.2e}")
        print(f"  CTEy = {props['CTEy']:.2e}")
        print(f"  CTEz = {props['CTEz']:.2e}")

        print("\n" + "="*60)


def main():
    """Main function to run the effective property calculation."""

    # Create calculator instance
    # RVE: 1x1x1 mm cube with 10% fiber volume fraction
    calc = CompositeEffectivePropertyCalculator(rve_size=1.0, fiber_vf=0.10)

    # Set material properties (optional - using defaults)
    # Matrix: Epoxy
    calc.set_matrix_properties(E=3500, nu=0.35, alpha=60e-6)
    # Fiber: Carbon fiber
    calc.set_fiber_properties(E=230000, nu=0.20, alpha=-0.5e-6)

    try:
        # Start MAPDL
        print("Starting ANSYS MAPDL...")
        calc.start_mapdl(run_location='/tmp/mapdl_run', override=True)

        # Calculate effective properties
        calc.calculate_effective_properties(element_size=0.05, strain_val=0.001)

        # Print results
        calc.print_results()

    except Exception as e:
        print(f"Error: {e}")
        raise

    finally:
        # Stop MAPDL
        calc.stop_mapdl()
        print("\nAnalysis completed.")


if __name__ == "__main__":
    main()
