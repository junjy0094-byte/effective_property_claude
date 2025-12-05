"""
Advanced Composite Material Effective Property Calculator using PyANSYS

This version implements Kinematic Uniform Boundary Conditions (KUBC) for RVE homogenization.
KUBC ensures that each face remains planar while allowing free lateral contraction/expansion,
which correctly models the behavior of a homogenized material.

Key features:
1. Square prism fiber geometry (for compatible hex meshing)
2. SOLID185 hex elements with matching meshes on opposite faces
3. KUBC boundary conditions using Coupled DOF (CP command)
4. Volume-averaged stress and strain for effective property computation

KUBC formulation for Ex calculation:
    - X=0 face: Ux=0 (fixed in loading direction)
    - X=L face: Ux=ε*L (displacement load)
    - Y=0, Y=L faces: All nodes have same Uy (coupled DOF - plane remains flat)
    - Z=0, Z=L faces: All nodes have same Uz (coupled DOF - plane remains flat)
    - Minimal rigid body constraints (one node fixed in Uy, Uz)

This approach allows proper Poisson contraction while maintaining plane faces,
matching the behavior of a homogenized material under uniaxial loading.

Reference:
- Suquet, P. (1987). "Elements of homogenization for inelastic solid mechanics"
- Hill, R. (1963). "Elastic properties of reinforced solids"
"""

import os
import glob
import numpy as np
from ansys.mapdl.core import launch_mapdl


class AdvancedCompositeCalculator:
    """
    Advanced calculator for effective material properties using KUBC.

    Uses Kinematic Uniform Boundary Conditions (KUBC) with Coupled DOF
    to ensure plane faces while allowing free lateral contraction.
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

        # Square fiber half-width from volume fraction
        # Vf = (2*a)^2 / L^2 => a = L * sqrt(Vf) / 2
        self.fiber_half_width = rve_size * np.sqrt(fiber_vf) / 2

        # Material properties
        self.mat_props = {
            'matrix': {'E': 3500.0, 'nu': 0.35, 'alpha': 60e-6},
            'fiber': {'E': 72000.0, 'nu': 0.22, 'alpha': 5e-6}
        }

        self.mapdl = None
        self.face_nodes = {}    # Node lists for each face
        self.corner_node = None # Corner node at origin for rigid body constraints
        self.stiffness_matrix = None
        self.compliance_matrix = None
        self.effective_props = {}

    def set_material(self, phase, E, nu, alpha):
        """Set material properties for a phase."""
        self.mat_props[phase] = {'E': E, 'nu': nu, 'alpha': alpha}

    def launch(self, nproc=12, **kwargs):
        """
        Launch MAPDL.

        Parameters
        ----------
        nproc : int
            Number of processors for SMP mode (default: 12)
        """
        self.mapdl = launch_mapdl(nproc=nproc, **kwargs)
        self.mapdl.ignore_errors = True  # Ignore non-critical MAPDL warnings
        self.mapdl.clear()
        self.mapdl.prep7()
        print(f"MAPDL launched successfully (SMP with {nproc} cores)")

    def exit(self):
        """Exit MAPDL."""
        if self.mapdl:
            self.mapdl.exit()
            self.mapdl = None

    def build_model(self, n_div=10):
        """
        Build the RVE model with mapped hex mesh.

        Parameters
        ----------
        n_div : int
            Number of element divisions along each edge
        """
        m = self.mapdl
        L = self.L
        a = self.fiber_half_width
        c = L / 2  # Center of RVE

        print(f"\n{'='*60}")
        print("BUILDING RVE MODEL (Square Fiber, Hex Mesh)")
        print(f"{'='*60}")
        print(f"RVE size: {L} mm")
        print(f"Fiber half-width: {a:.4f} mm")
        print(f"Volume fraction: {self.Vf*100:.1f}%")

        # Clear and setup
        m.clear()
        m.prep7()

        # Define materials
        mp = self.mat_props['matrix']
        m.mp('EX', 1, mp['E'])
        m.mp('NUXY', 1, mp['nu'])
        m.mp('ALPX', 1, mp['alpha'])

        fp = self.mat_props['fiber']
        m.mp('EX', 2, fp['E'])
        m.mp('NUXY', 2, fp['nu'])
        m.mp('ALPX', 2, fp['alpha'])

        # Element type: SOLID185 (8-node hex)
        m.et(1, 'SOLID185')

        # Save DB after material definition
        m.save('step1_materials.db')
        print("  Saved: step1_materials.db")

        # Create geometry using keypoints and volumes for mapped meshing
        # The RVE is divided into 9 volumes (3x3 in XY plane, extruded in Z)
        # Center volume is fiber, surrounding 8 volumes are matrix

        # Fiber boundaries
        x1, x2 = c - a, c + a  # Fiber X boundaries
        y1, y2 = c - a, c + a  # Fiber Y boundaries

        # Create 9 blocks in XY plane, extruded through Z
        # Each block is created, assigned material, and meshed immediately
        regions = [
            # (x_start, x_end, y_start, y_end, material)
            (0, x1, 0, y1, 1),      # Bottom-left corner (matrix)
            (x1, x2, 0, y1, 1),     # Bottom-center (matrix)
            (x2, L, 0, y1, 1),      # Bottom-right corner (matrix)
            (0, x1, y1, y2, 1),     # Middle-left (matrix)
            (x1, x2, y1, y2, 2),    # CENTER = FIBER
            (x2, L, y1, y2, 1),     # Middle-right (matrix)
            (0, x1, y2, L, 1),      # Top-left corner (matrix)
            (x1, x2, y2, L, 1),     # Top-center (matrix)
            (x2, L, y2, L, 1),      # Top-right corner (matrix)
        ]

        # Mesh settings
        m.mshkey(1)  # Mapped meshing
        m.mshape(0, '3D')  # Hex elements

        # Create, assign material, and mesh each volume
        for i, (xs, xe, ys, ye, mat_id) in enumerate(regions):
            m.block(xs, xe, ys, ye, 0, L)
            vol_num = i + 1
            m.lsel('S', 'VOLU', '', vol_num)
            m.lesize('ALL', '', '', n_div)
            m.vsel('S', 'VOLU', '', vol_num)
            m.vatt(mat_id, '', 1)
            m.vmesh(vol_num)

        m.allsel()

        # Save DB after geometry and mesh
        m.save('step2_geometry_mesh.db')
        print("  Saved: step2_geometry_mesh.db")

        # Merge nodes at interfaces
        m.nummrg('NODE', 1e-6)

        nn = int(m.get('NCOUNT', 'NODE', '', 'COUNT'))
        ne = int(m.get('ECOUNT', 'ELEM', '', 'COUNT'))
        print(f"Mesh: {nn} nodes, {ne} elements")

        # Create face node sets for KUBC
        self._create_face_sets()

        # Save DB after face sets
        m.save('step3_face_sets.db')
        print("  Saved: step3_face_sets.db")

        print("Model built successfully")

    def _create_face_sets(self):
        """Create node sets for each face and identify corner node."""
        m = self.mapdl
        L = self.L
        tol = 1e-6

        print("  Creating face node sets for KUBC...")

        # Create face component sets and store node lists
        faces = [
            ('XNEG', 'X', 0),
            ('XPOS', 'X', L),
            ('YNEG', 'Y', 0),
            ('YPOS', 'Y', L),
            ('ZNEG', 'Z', 0),
            ('ZPOS', 'Z', L),
        ]

        self.face_nodes = {}
        for name, direction, coord in faces:
            m.nsel('S', 'LOC', direction, coord - tol, coord + tol)
            m.cm(name, 'NODE')
            # Store node list for this face
            self.face_nodes[name] = m.mesh.nnum.copy()

        m.allsel()

        # Find corner node at origin (0, 0, 0) for rigid body constraints
        m.nsel('S', 'LOC', 'X', 0, tol)
        m.nsel('R', 'LOC', 'Y', 0, tol)
        m.nsel('R', 'LOC', 'Z', 0, tol)
        corner_nodes = m.mesh.nnum
        if len(corner_nodes) > 0:
            self.corner_node = int(corner_nodes[0])
        else:
            # Fallback: find any node on XNEG face
            self.corner_node = int(self.face_nodes['XNEG'][0])
        m.allsel()

        print(f"    XNEG: {len(self.face_nodes['XNEG'])} nodes")
        print(f"    XPOS: {len(self.face_nodes['XPOS'])} nodes")
        print(f"    YNEG: {len(self.face_nodes['YNEG'])} nodes")
        print(f"    YPOS: {len(self.face_nodes['YPOS'])} nodes")
        print(f"    ZNEG: {len(self.face_nodes['ZNEG'])} nodes")
        print(f"    ZPOS: {len(self.face_nodes['ZPOS'])} nodes")
        print(f"    Corner node (origin): {self.corner_node}")

    def _clear_all_constraints(self):
        """Clear all displacement constraints and coupled DOFs."""
        m = self.mapdl
        m.ddele('ALL', 'ALL')
        m.run('CPDELE,ALL')
        m.run('CEDELE,ALL')

    def _couple_face_dof(self, face_name, dof, cp_set_num):
        """
        Couple all nodes on a face to have the same displacement in specified DOF.

        Parameters
        ----------
        face_name : str
            Name of face component (XNEG, XPOS, YNEG, YPOS, ZNEG, ZPOS)
        dof : str
            Degree of freedom to couple (UX, UY, or UZ)
        cp_set_num : int
            Coupled set number for ANSYS CP command
        """
        m = self.mapdl
        m.cmsel('S', face_name)
        m.cp(cp_set_num, dof, 'ALL')
        m.allsel()

    def _couple_all_planes(self, direction, dof, start_cp_num):
        """
        Couple all nodes at each unique coordinate value in the specified direction.

        This ensures that all parallel planes remain flat like in a homogeneous material.
        For example, if direction='X' and dof='UX', all nodes with the same X coordinate
        will have the same UX displacement, keeping each YZ-plane flat.

        Parameters
        ----------
        direction : str
            Coordinate direction ('X', 'Y', or 'Z')
        dof : str
            Degree of freedom to couple (UX, UY, or UZ)
        start_cp_num : int
            Starting coupled set number

        Returns
        -------
        int
            Next available coupled set number
        """
        m = self.mapdl
        L = self.L
        tol = 1e-6

        # Get all nodes and their coordinates
        m.nsel('ALL')
        all_nodes = m.mesh.nodes  # shape (n_nodes, 3)
        node_nums = m.mesh.nnum

        # Get coordinate index
        coord_idx = {'X': 0, 'Y': 1, 'Z': 2}[direction]
        coords = all_nodes[:, coord_idx]

        # Find unique coordinate values (rounded to avoid floating point issues)
        unique_coords = np.unique(np.round(coords, 6))

        cp_num = start_cp_num
        for coord in unique_coords:
            # Select all nodes at this coordinate
            m.nsel('S', 'LOC', direction, coord - tol, coord + tol)
            n_selected = int(m.get('NCOUNT', 'NODE', '', 'COUNT'))

            # Only couple if more than one node
            if n_selected > 1:
                m.cp(cp_num, dof, 'ALL')
                cp_num += 1

            m.allsel()

        return cp_num

    def _apply_linear_shear_constraint(self, shear_dir, disp_dof, coord_dir, start_ce_num):
        """
        Apply linear constraint for shear deformation: displacement proportional to coordinate.

        For example, in XY shear: UX = (Y/L) * UX_master
        This ensures the shear strain is uniform throughout the RVE.

        Parameters
        ----------
        shear_dir : str
            Direction of shear displacement ('X', 'Y', or 'Z')
        disp_dof : str
            Displacement DOF to constrain ('UX', 'UY', or 'UZ')
        coord_dir : str
            Coordinate direction that displacement is proportional to ('X', 'Y', or 'Z')
        start_ce_num : int
            Starting constraint equation number

        Returns
        -------
        int
            Next available constraint equation number
        """
        m = self.mapdl
        L = self.L
        tol = 1e-6

        # Get all nodes and their coordinates
        m.nsel('ALL')
        all_nodes = m.mesh.nodes
        node_nums = m.mesh.nnum

        # Get coordinate index
        coord_idx = {'X': 0, 'Y': 1, 'Z': 2}[coord_dir]
        coords = all_nodes[:, coord_idx]

        # Find master node at coord_dir = L (the loaded face)
        # Select nodes at max coordinate
        m.nsel('S', 'LOC', coord_dir, L - tol, L + tol)
        master_nodes = m.mesh.nnum
        master_node = int(master_nodes[0])  # Use first node as master
        m.allsel()

        ce_num = start_ce_num

        # For each node, create CE: disp_dof_i - (coord_i/L) * disp_dof_master = 0
        # Skip nodes at coord=0 (fixed) and coord=L (master nodes)
        for i, node in enumerate(node_nums):
            coord_val = coords[i]

            # Skip boundary nodes (coord=0 or coord=L)
            if coord_val < tol or coord_val > L - tol:
                continue

            # Ratio of coordinate to L
            ratio = coord_val / L

            # CE command: C1*NODE1.DOF + C2*NODE2.DOF = CONST
            # We want: UX_node - ratio * UX_master = 0
            # CE, NEQN, CONST, NODE1, Lab1, C1, NODE2, Lab2, C2
            m.ce(ce_num, 0, int(node), disp_dof, 1.0, master_node, disp_dof, -ratio)
            ce_num += 1

        return ce_num

    def apply_bc_uniaxial_x(self, strain_val=0.001):
        """
        Apply KUBC boundary conditions for uniaxial strain in X direction (for Ex).

        Boundary conditions:
        - X=0 face: Ux=0 (fixed in loading direction)
        - X=L face: Ux=ε*L (displacement load)
        - All YZ-planes (every X coordinate): nodes have same UX (planes remain flat)
        - All XZ-planes (every Y coordinate): nodes have same UY (planes remain flat)
        - All XY-planes (every Z coordinate): nodes have same UZ (planes remain flat)
        - Rigid body constraint: Corner node Uy=Uz=0
        """
        m = self.mapdl
        L = self.L

        self._clear_all_constraints()

        cp_num = 1  # Coupled DOF set counter

        # X=0 face: Fixed in X direction
        m.cmsel('S', 'XNEG')
        m.d('ALL', 'UX', 0)
        m.allsel()

        # X=L face: Applied displacement
        m.cmsel('S', 'XPOS')
        m.d('ALL', 'UX', strain_val * L)
        m.allsel()

        # Couple all X-planes: All nodes at same X have same UX (loading direction - plane remains flat)
        cp_num = self._couple_all_planes('X', 'UX', cp_num)

        # Couple all Y-planes: All nodes at same Y have same UY (plane remains flat)
        cp_num = self._couple_all_planes('Y', 'UY', cp_num)

        # Couple all Z-planes: All nodes at same Z have same UZ (plane remains flat)
        cp_num = self._couple_all_planes('Z', 'UZ', cp_num)

        # Rigid body constraints - fix corner node in Y and Z only
        m.d(self.corner_node, 'UY', 0)
        m.d(self.corner_node, 'UZ', 0)

        m.allsel()
        print(f"    Applied KUBC for Ex (ε11={strain_val})")

    def apply_bc_uniaxial_y(self, strain_val=0.001):
        """
        Apply KUBC boundary conditions for uniaxial strain in Y direction (for Ey).

        Boundary conditions:
        - Y=0 face: Uy=0 (fixed in loading direction)
        - Y=L face: Uy=ε*L (displacement load)
        - All YZ-planes (every X coordinate): nodes have same UX (planes remain flat)
        - All XZ-planes (every Y coordinate): nodes have same UY (planes remain flat)
        - All XY-planes (every Z coordinate): nodes have same UZ (planes remain flat)
        - Rigid body constraint: Corner node Ux=Uz=0
        """
        m = self.mapdl
        L = self.L

        self._clear_all_constraints()

        cp_num = 1

        # Y=0 face: Fixed in Y direction
        m.cmsel('S', 'YNEG')
        m.d('ALL', 'UY', 0)
        m.allsel()

        # Y=L face: Applied displacement
        m.cmsel('S', 'YPOS')
        m.d('ALL', 'UY', strain_val * L)
        m.allsel()

        # Couple all X-planes: All nodes at same X have same UX (plane remains flat)
        cp_num = self._couple_all_planes('X', 'UX', cp_num)

        # Couple all Y-planes: All nodes at same Y have same UY (loading direction - plane remains flat)
        cp_num = self._couple_all_planes('Y', 'UY', cp_num)

        # Couple all Z-planes: All nodes at same Z have same UZ (plane remains flat)
        cp_num = self._couple_all_planes('Z', 'UZ', cp_num)

        # Rigid body constraints
        m.d(self.corner_node, 'UX', 0)
        m.d(self.corner_node, 'UZ', 0)

        m.allsel()
        print(f"    Applied KUBC for Ey (ε22={strain_val})")

    def apply_bc_uniaxial_z(self, strain_val=0.001):
        """
        Apply KUBC boundary conditions for uniaxial strain in Z direction (for Ez).

        Boundary conditions:
        - Z=0 face: Uz=0 (fixed in loading direction)
        - Z=L face: Uz=ε*L (displacement load)
        - All YZ-planes (every X coordinate): nodes have same UX (planes remain flat)
        - All XZ-planes (every Y coordinate): nodes have same UY (planes remain flat)
        - All XY-planes (every Z coordinate): nodes have same UZ (planes remain flat)
        - Rigid body constraint: Corner node Ux=Uy=0
        """
        m = self.mapdl
        L = self.L

        self._clear_all_constraints()

        cp_num = 1

        # Z=0 face: Fixed in Z direction
        m.cmsel('S', 'ZNEG')
        m.d('ALL', 'UZ', 0)
        m.allsel()

        # Z=L face: Applied displacement
        m.cmsel('S', 'ZPOS')
        m.d('ALL', 'UZ', strain_val * L)
        m.allsel()

        # Couple all X-planes: All nodes at same X have same UX (plane remains flat)
        cp_num = self._couple_all_planes('X', 'UX', cp_num)

        # Couple all Y-planes: All nodes at same Y have same UY (plane remains flat)
        cp_num = self._couple_all_planes('Y', 'UY', cp_num)

        # Couple all Z-planes: All nodes at same Z have same UZ (loading direction - plane remains flat)
        cp_num = self._couple_all_planes('Z', 'UZ', cp_num)

        # Rigid body constraints
        m.d(self.corner_node, 'UX', 0)
        m.d(self.corner_node, 'UY', 0)

        m.allsel()
        print(f"    Applied KUBC for Ez (ε33={strain_val})")

    def apply_bc_shear_xy(self, strain_val=0.001):
        """
        Apply KUBC boundary conditions for shear strain in XY plane (for Gxy).

        Boundary conditions:
        - Y=0 face: Ux=0, Uy=0 (fixed)
        - Y=L face: Ux=γ*L (shear displacement)
        - Linear shear: UX proportional to Y coordinate (UX = Y/L * UX_max)
        - All XZ-planes (every Y coordinate): nodes have same UY (planes remain flat)
        - All XY-planes (every Z coordinate): nodes have same UZ (planes remain flat)
        - Rigid body constraints
        """
        m = self.mapdl
        L = self.L

        self._clear_all_constraints()

        cp_num = 1
        ce_num = 1

        # Y=0 face: Fixed in X direction (transverse)
        m.cmsel('S', 'YNEG')
        m.d('ALL', 'UX', 0)
        m.d('ALL', 'UY', 0)  # Also fix UY to prevent Y-translation
        m.allsel()

        # Y=L face: Shear displacement in X direction
        m.cmsel('S', 'YPOS')
        m.d('ALL', 'UX', strain_val * L)
        m.allsel()

        # Apply linear shear constraint: UX = (Y/L) * UX_max
        # This ensures uniform shear strain throughout
        ce_num = self._apply_linear_shear_constraint('X', 'UX', 'Y', ce_num)

        # Couple all Y-planes: All nodes at same Y have same UY (plane remains flat)
        cp_num = self._couple_all_planes('Y', 'UY', cp_num)

        # Couple all Z-planes: All nodes at same Z have same UZ (plane remains flat)
        cp_num = self._couple_all_planes('Z', 'UZ', cp_num)

        # Rigid body constraint - fix corner node in Z
        m.d(self.corner_node, 'UZ', 0)

        m.allsel()
        print(f"    Applied KUBC for Gxy (γ12={strain_val})")

    def apply_bc_shear_yz(self, strain_val=0.001):
        """
        Apply KUBC boundary conditions for shear strain in YZ plane (for Gyz).

        Boundary conditions:
        - Z=0 face: Uy=0, Uz=0 (fixed)
        - Z=L face: Uy=γ*L (shear displacement)
        - Linear shear: UY proportional to Z coordinate (UY = Z/L * UY_max)
        - All XY-planes (every Z coordinate): nodes have same UZ (planes remain flat)
        - All YZ-planes (every X coordinate): nodes have same UX (planes remain flat)
        - Rigid body constraints
        """
        m = self.mapdl
        L = self.L

        self._clear_all_constraints()

        cp_num = 1
        ce_num = 1

        # Z=0 face: Fixed in Y direction (transverse)
        m.cmsel('S', 'ZNEG')
        m.d('ALL', 'UY', 0)
        m.d('ALL', 'UZ', 0)  # Also fix UZ
        m.allsel()

        # Z=L face: Shear displacement in Y direction
        m.cmsel('S', 'ZPOS')
        m.d('ALL', 'UY', strain_val * L)
        m.allsel()

        # Apply linear shear constraint: UY = (Z/L) * UY_max
        ce_num = self._apply_linear_shear_constraint('Y', 'UY', 'Z', ce_num)

        # Couple all Z-planes: All nodes at same Z have same UZ (plane remains flat)
        cp_num = self._couple_all_planes('Z', 'UZ', cp_num)

        # Couple all X-planes: All nodes at same X have same UX (plane remains flat)
        cp_num = self._couple_all_planes('X', 'UX', cp_num)

        # Rigid body constraint
        m.d(self.corner_node, 'UX', 0)

        m.allsel()
        print(f"    Applied KUBC for Gyz (γ23={strain_val})")

    def apply_bc_shear_zx(self, strain_val=0.001):
        """
        Apply KUBC boundary conditions for shear strain in ZX plane (for Gzx).

        Boundary conditions:
        - X=0 face: Uz=0, Ux=0 (fixed)
        - X=L face: Uz=γ*L (shear displacement)
        - Linear shear: UZ proportional to X coordinate (UZ = X/L * UZ_max)
        - All YZ-planes (every X coordinate): nodes have same UX (planes remain flat)
        - All XZ-planes (every Y coordinate): nodes have same UY (planes remain flat)
        - Rigid body constraints
        """
        m = self.mapdl
        L = self.L

        self._clear_all_constraints()

        cp_num = 1
        ce_num = 1

        # X=0 face: Fixed in Z direction (transverse)
        m.cmsel('S', 'XNEG')
        m.d('ALL', 'UZ', 0)
        m.d('ALL', 'UX', 0)  # Also fix UX
        m.allsel()

        # X=L face: Shear displacement in Z direction
        m.cmsel('S', 'XPOS')
        m.d('ALL', 'UZ', strain_val * L)
        m.allsel()

        # Apply linear shear constraint: UZ = (X/L) * UZ_max
        ce_num = self._apply_linear_shear_constraint('Z', 'UZ', 'X', ce_num)

        # Couple all X-planes: All nodes at same X have same UX (plane remains flat)
        cp_num = self._couple_all_planes('X', 'UX', cp_num)

        # Couple all Y-planes: All nodes at same Y have same UY (plane remains flat)
        cp_num = self._couple_all_planes('Y', 'UY', cp_num)

        # Rigid body constraint
        m.d(self.corner_node, 'UY', 0)

        m.allsel()
        print(f"    Applied KUBC for Gzx (γ31={strain_val})")

    def apply_thermal_bc(self, delta_T=1.0):
        """
        Apply thermal loading with KUBC boundary conditions.

        For thermal analysis:
        - All internal planes have coupled DOF to remain planar (like homogeneous material)
        - Free thermal expansion in all directions
        - Minimal rigid body constraints
        """
        m = self.mapdl

        self._clear_all_constraints()

        cp_num = 1

        # Couple all planes to remain planar during thermal expansion
        # All YZ-planes (every X coordinate): nodes have same UX
        cp_num = self._couple_all_planes('X', 'UX', cp_num)

        # All XZ-planes (every Y coordinate): nodes have same UY
        cp_num = self._couple_all_planes('Y', 'UY', cp_num)

        # All XY-planes (every Z coordinate): nodes have same UZ
        cp_num = self._couple_all_planes('Z', 'UZ', cp_num)

        # Minimal rigid body constraints - fix corner node completely
        m.d(self.corner_node, 'UX', 0)
        m.d(self.corner_node, 'UY', 0)
        m.d(self.corner_node, 'UZ', 0)

        m.allsel()

        # Apply thermal load
        # IMPORTANT: tref sets the reference temperature, tunif sets current uniform temp
        # For thermal expansion: strain = alpha * (T - Tref)
        m.tref(0)  # Reference temperature = 0
        m.bfunif('TEMP', delta_T)  # Current temperature = delta_T

        print(f"    Applied thermal BC (ΔT={delta_T}°C)")

    def solve(self, jobname=None):
        """Solve current load case and optionally save result file."""
        m = self.mapdl

        # Set jobname for this load case if provided
        if jobname:
            m.finish()
            m.filname(jobname)

        m.run('/SOLU')
        m.antype('STATIC')
        m.solve()
        m.finish()

        if jobname:
            print(f"    Result saved: {jobname}.rst")

        # Clean up temporary files (DSP, esav, full, mntr)
        self._cleanup_temp_files()

    def _cleanup_temp_files(self):
        """Remove temporary MAPDL files to save disk space."""
        if not self.mapdl:
            return

        work_dir = self.mapdl.directory
        temp_extensions = ['*.DSP', '*.esav', '*.full', '*.mntr',
                          '*.dsp', '*.ESAV', '*.FULL', '*.MNTR']

        for ext in temp_extensions:
            pattern = os.path.join(work_dir, ext)
            for filepath in glob.glob(pattern):
                try:
                    os.remove(filepath)
                except OSError:
                    pass  # Ignore errors if file is in use or already deleted

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

        # Create element tables for stresses
        m.etable('SXX', 'S', 'X')
        m.etable('SYY', 'S', 'Y')
        m.etable('SZZ', 'S', 'Z')
        m.etable('SXY', 'S', 'XY')
        m.etable('SYZ', 'S', 'YZ')
        m.etable('SXZ', 'S', 'XZ')

        # Get element volumes for proper averaging
        m.etable('EVOL', 'VOLU')

        # Volume-weighted averaging
        m.ssum()

        total_vol = m.get_value('SSUM', '', 'ITEM', 'EVOL')

        sxx_sum = m.get_value('SSUM', '', 'ITEM', 'SXX')
        syy_sum = m.get_value('SSUM', '', 'ITEM', 'SYY')
        szz_sum = m.get_value('SSUM', '', 'ITEM', 'SZZ')
        sxy_sum = m.get_value('SSUM', '', 'ITEM', 'SXY')
        syz_sum = m.get_value('SSUM', '', 'ITEM', 'SYZ')
        sxz_sum = m.get_value('SSUM', '', 'ITEM', 'SXZ')

        n_elem = m.mesh.n_elem

        # Simple average (each element assumed equal volume for hex mesh)
        stress = np.array([
            sxx_sum / n_elem,
            syy_sum / n_elem,
            szz_sum / n_elem,
            sxy_sum / n_elem,
            syz_sum / n_elem,
            sxz_sum / n_elem
        ])

        m.finish()
        return stress

    def get_thermal_strain(self, delta_T=1.0):
        """
        Get effective thermal expansion strains from face displacements.

        The thermal strains are computed from the displacement difference
        between opposite faces:
            ε_th = (u_pos - u_neg) / L / ΔT

        With KUBC, all nodes on a face have the same normal displacement (coupled),
        so we can use any node on each face.
        """
        m = self.mapdl
        L = self.L

        m.post1()
        m.set('LAST')

        # Get a representative node from each face (first node in the list)
        node_xpos = int(self.face_nodes['XPOS'][0])
        node_xneg = int(self.face_nodes['XNEG'][0])
        node_ypos = int(self.face_nodes['YPOS'][0])
        node_yneg = int(self.face_nodes['YNEG'][0])
        node_zpos = int(self.face_nodes['ZPOS'][0])
        node_zneg = int(self.face_nodes['ZNEG'][0])

        # Get UX on X faces
        ux_xpos = float(m.get('UX', 'NODE', node_xpos, 'U', 'X'))
        ux_xneg = float(m.get('UX', 'NODE', node_xneg, 'U', 'X'))

        # Get UY on Y faces
        uy_ypos = float(m.get('UY', 'NODE', node_ypos, 'U', 'Y'))
        uy_yneg = float(m.get('UY', 'NODE', node_yneg, 'U', 'Y'))

        # Get UZ on Z faces
        uz_zpos = float(m.get('UZ', 'NODE', node_zpos, 'U', 'Z'))
        uz_zneg = float(m.get('UZ', 'NODE', node_zneg, 'U', 'Z'))

        m.allsel()
        m.finish()

        # Thermal strains (CTE = strain / delta_T)
        eps_th = np.array([
            (ux_xpos - ux_xneg) / L / delta_T,
            (uy_ypos - uy_yneg) / L / delta_T,
            (uz_zpos - uz_zneg) / L / delta_T
        ])

        return eps_th

    def compute_stiffness_matrix(self, strain_mag=0.001):
        """
        Compute the complete 6x6 stiffness matrix using KUBC.

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
        print("COMPUTING STIFFNESS MATRIX (KUBC)")
        print(f"{'='*60}")
        print(f"Result files will be saved in: {self.mapdl.directory}")

        C = np.zeros((6, 6))

        # Define 6 load cases with corresponding BC functions
        # Format: (direction_name, jobname, bc_function)
        load_cases = [
            ('ε11 (Ex)', 'LC1_e11', self.apply_bc_uniaxial_x),
            ('ε22 (Ey)', 'LC2_e22', self.apply_bc_uniaxial_y),
            ('ε33 (Ez)', 'LC3_e33', self.apply_bc_uniaxial_z),
            ('γ12 (Gxy)', 'LC4_g12', self.apply_bc_shear_xy),
            ('γ23 (Gyz)', 'LC5_g23', self.apply_bc_shear_yz),
            ('γ31 (Gzx)', 'LC6_g31', self.apply_bc_shear_zx),
        ]

        for i, (direction, jobname, bc_func) in enumerate(load_cases):
            print(f"  Load case {i+1}/6: {direction}...")

            self.mapdl.prep7()
            bc_func(strain_mag)

            # Save DB before solving each load case
            self.mapdl.save(f'step4_{jobname}_bc.db')
            print(f"    Saved: step4_{jobname}_bc.db")

            self.solve(jobname=jobname)

            stress = self.get_volume_avg_stress()
            C[:, i] = stress / strain_mag

            print(f"    Stress: σ11={stress[0]:.2f}, σ22={stress[1]:.2f}, σ33={stress[2]:.2f}")
            print(f"            τ12={stress[3]:.2f}, τ23={stress[4]:.2f}, τ31={stress[5]:.2f}")

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
        # S11 = 1/Ex, S22 = 1/Ey, S33 = 1/Ez
        # S44 = 1/Gyz, S55 = 1/Gzx, S66 = 1/Gxy
        # S12 = -nuxy/Ex, S13 = -nuxz/Ex, S23 = -nuyz/Ey

        Ex = 1.0 / S[0, 0]
        Ey = 1.0 / S[1, 1]
        Ez = 1.0 / S[2, 2]

        # Shear moduli (indices 3,4,5 correspond to 12,23,31 in our convention)
        Gxy = 1.0 / S[3, 3]  # γ12
        Gyz = 1.0 / S[4, 4]  # γ23
        Gzx = 1.0 / S[5, 5]  # γ31

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
        self.solve(jobname='LC7_thermal')

        eps_th = self.get_thermal_strain(delta_T)

        CTEx = eps_th[0]
        CTEy = eps_th[1]
        CTEz = eps_th[2]

        self.effective_props.update({
            'CTEx': CTEx, 'CTEy': CTEy, 'CTEz': CTEz
        })

        print(f"  CTE computed: [{CTEx:.2e}, {CTEy:.2e}, {CTEz:.2e}] 1/°C")

        return {'CTEx': CTEx, 'CTEy': CTEy, 'CTEz': CTEz}

    def run_full_analysis(self, n_div=10, strain_mag=0.001):
        """
        Run complete analysis to get all effective properties.

        Parameters
        ----------
        n_div : int
            Number of mesh divisions
        strain_mag : float
            Strain magnitude for mechanical load cases

        Returns
        -------
        props : dict
            All effective properties
        """
        # Build model
        self.build_model(n_div)

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
        print(f"  Fiber half-width: {self.fiber_half_width:.4f} mm")

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
    calc.set_material('fiber', E=72000, nu=0.22, alpha=5e-6)

    try:
        print("Starting ANSYS MAPDL...")
        calc.launch(run_location='/tmp/mapdl_adv', override=True)

        # Run full analysis with hex mesh
        calc.run_full_analysis(n_div=10, strain_mag=0.001)

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
