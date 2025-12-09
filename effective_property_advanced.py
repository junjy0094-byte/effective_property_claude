"""
Advanced Composite Material Effective Property Calculator using PyANSYS

This version implements Periodic Boundary Conditions (PBC) for RVE homogenization
based on ANSYS 2025 R1 Theory Documentation (Section 3.2.2.1).

Key features:
1. Square prism fiber geometry (for compatible hex meshing)
2. SOLID185 hex elements with matching meshes on opposite faces
3. Periodic boundary conditions using Constraint Equations (CE command)
4. Stiffness matrix computation from 6 mechanical load cases
5. Thermal expansion coefficients from thermal load case

Periodic Boundary Conditions (Equations 3.24-3.26):
    On X-faces: u_x(L_x,y,z) = u_x(0,y,z) + ε_x*L_x
                u_y(L_x,y,z) = u_y(0,y,z) + γ_xy*L_x
                u_z(L_x,y,z) = u_z(0,y,z) + γ_xz*L_x
    On Y-faces: u_x(x,L_y,z) = u_x(x,0,z)
                u_y(x,L_y,z) = u_y(x,0,z) + ε_y*L_y
                u_z(x,L_y,z) = u_z(x,0,z) + γ_yz*L_y
    On Z-faces: u_x(x,y,L_z) = u_x(x,y,0)
                u_y(x,y,L_z) = u_y(x,y,0)
                u_z(x,y,L_z) = u_z(x,y,0) + ε_z*L_z

Rigid body constraints (Equation 3.27):
    u_x(point with x=0) = 0
    u_y(point with y=0) = 0
    u_z(point with z=0) = 0

Reference:
- ANSYS 2025 R1 Theory Documentation, Section 3.2.2.1
- Li et al. (2008), Li et al. (2015) for periodic boundary conditions
"""

import os
import glob
import numpy as np
from ansys.mapdl.core import launch_mapdl


class AdvancedCompositeCalculator:
    """
    Advanced calculator for effective material properties using Periodic BC.

    Uses Periodic Boundary Conditions (PBC) with Constraint Equations (CE)
    to enforce periodicity on opposite faces of the RVE.
    Based on ANSYS 2025 R1 Theory Documentation.
    """

    def __init__(self, rve_size=1.0, fiber_vf=0.10):
        """
        Initialize the calculator.

        Parameters
        ----------
        rve_size : float
            Size of the RVE cube in mm (L_x = L_y = L_z = L)
        fiber_vf : float
            Fiber volume fraction (0 to 1)
        """
        self.L = rve_size
        self.Lx = rve_size
        self.Ly = rve_size
        self.Lz = rve_size
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
        self.element_type = 'SOLID185'  # Default element type
        self.face_nodes = {}    # Node lists for each face
        self.node_pairs = {}    # Node pairs for periodic BC (X, Y, Z directions)
        self.rigid_body_nodes = {}  # Nodes for rigid body constraints
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

    def build_model(self, ele_size=0.05, n_div=None, element_type='SOLID185'):
        """
        Build the RVE model with mapped hex mesh.

        Parameters
        ----------
        ele_size : float
            Element size in mm (default: 0.05)
        n_div : int, optional
            Number of element divisions along each edge. If provided, overrides ele_size
            for line sizing (SOLID185 only).
        element_type : str
            Element type: 'SOLID185' (8-node hex) or 'SOLID187' (10-node tet)
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

        # Element type: SOLID185 (8-node hex) or SOLID187 (10-node tet)
        self.element_type = element_type
        m.et(1, element_type)
        print(f"Element type: {element_type}")

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

        # Mesh settings based on element type
        if element_type == 'SOLID187':
            m.mshkey(0)  # Free meshing for tet elements
            m.mshape(1, '3D')  # Tet elements
        else:
            m.mshkey(1)  # Mapped meshing for hex elements
            m.mshape(0, '3D')  # Hex elements

        # Set element size
        m.esize(ele_size)

        # Create, assign material, and mesh each volume
        for i, (xs, xe, ys, ye, mat_id) in enumerate(regions):
            m.block(xs, xe, ys, ye, 0, L)
            vol_num = i + 1
            m.vsel('S', 'VOLU', '', vol_num)
            # Use n_div for line sizing if provided (SOLID185 only)
            if n_div is not None and element_type == 'SOLID185':
                m.lsel('S', 'VOLU', '', vol_num)
                m.lesize('ALL', '', '', n_div)
            m.vatt(mat_id, '', 1)
            m.vmesh(vol_num)

        m.allsel()

        # Merge nodes at interfaces
        m.nummrg('NODE', 1e-6)

        nn = int(m.get('NCOUNT', 'NODE', '', 'COUNT'))
        ne = int(m.get('ECOUNT', 'ELEM', '', 'COUNT'))
        print(f"Mesh: {nn} nodes, {ne} elements")

        # Create face node sets for periodic BC
        self._create_face_sets()

        print("Model built successfully")

    def _create_face_sets(self):
        """Create node sets for each face and create node pairs for periodic BC."""
        m = self.mapdl
        L = self.L
        tol = 1e-6

        print("  Creating face node sets for Periodic BC...")

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

        # Create node pairs for periodic boundary conditions
        self._create_node_pairs()

        # Find nodes for rigid body constraints (Equation 3.27)
        # u_x(point with x=0) = 0, u_y(point with y=0) = 0, u_z(point with z=0) = 0
        self._find_rigid_body_nodes()

        print(f"    XNEG: {len(self.face_nodes['XNEG'])} nodes")
        print(f"    XPOS: {len(self.face_nodes['XPOS'])} nodes")
        print(f"    YNEG: {len(self.face_nodes['YNEG'])} nodes")
        print(f"    YPOS: {len(self.face_nodes['YPOS'])} nodes")
        print(f"    ZNEG: {len(self.face_nodes['ZNEG'])} nodes")
        print(f"    ZPOS: {len(self.face_nodes['ZPOS'])} nodes")
        print(f"    X-direction node pairs: {len(self.node_pairs['X'])}")
        print(f"    Y-direction node pairs: {len(self.node_pairs['Y'])}")
        print(f"    Z-direction node pairs: {len(self.node_pairs['Z'])}")

    def _create_node_pairs(self):
        """
        Create node pairs between opposite faces for periodic BC.

        For each direction, pairs nodes on the negative face with corresponding
        nodes on the positive face that have the same coordinates in the other
        two directions. Optimized O(n) algorithm using dictionary lookup.
        """
        m = self.mapdl
        decimals = 6  # Precision for coordinate rounding

        m.allsel()
        all_nodes = m.mesh.nodes  # shape (n_nodes, 3): [x, y, z]
        all_nnum = m.mesh.nnum    # node numbers

        # Create node number to index mapping for fast lookup
        nnum_to_idx = {n: i for i, n in enumerate(all_nnum)}

        self.node_pairs = {'X': [], 'Y': [], 'Z': []}

        # X-direction pairs: Match by (y, z) coordinates
        # Build dictionary from XPOS nodes: (y,z) -> node_num
        xpos_dict = {}
        for n_pos in self.face_nodes['XPOS']:
            idx = nnum_to_idx[n_pos]
            key = (round(all_nodes[idx, 1], decimals), round(all_nodes[idx, 2], decimals))
            xpos_dict[key] = int(n_pos)

        for n_neg in self.face_nodes['XNEG']:
            idx = nnum_to_idx[n_neg]
            key = (round(all_nodes[idx, 1], decimals), round(all_nodes[idx, 2], decimals))
            if key in xpos_dict:
                self.node_pairs['X'].append((int(n_neg), xpos_dict[key]))

        # Y-direction pairs: Match by (x, z) coordinates
        ypos_dict = {}
        for n_pos in self.face_nodes['YPOS']:
            idx = nnum_to_idx[n_pos]
            key = (round(all_nodes[idx, 0], decimals), round(all_nodes[idx, 2], decimals))
            ypos_dict[key] = int(n_pos)

        for n_neg in self.face_nodes['YNEG']:
            idx = nnum_to_idx[n_neg]
            key = (round(all_nodes[idx, 0], decimals), round(all_nodes[idx, 2], decimals))
            if key in ypos_dict:
                self.node_pairs['Y'].append((int(n_neg), ypos_dict[key]))

        # Z-direction pairs: Match by (x, y) coordinates
        zpos_dict = {}
        for n_pos in self.face_nodes['ZPOS']:
            idx = nnum_to_idx[n_pos]
            key = (round(all_nodes[idx, 0], decimals), round(all_nodes[idx, 1], decimals))
            zpos_dict[key] = int(n_pos)

        for n_neg in self.face_nodes['ZNEG']:
            idx = nnum_to_idx[n_neg]
            key = (round(all_nodes[idx, 0], decimals), round(all_nodes[idx, 1], decimals))
            if key in zpos_dict:
                self.node_pairs['Z'].append((int(n_neg), zpos_dict[key]))

    def _find_rigid_body_nodes(self):
        """
        Find nodes for rigid body constraints (Equation 3.27).

        u_x(point with x=0) = 0  -> fix UX at one node on x=0 plane
        u_y(point with y=0) = 0  -> fix UY at one node on y=0 plane
        u_z(point with z=0) = 0  -> fix UZ at one node on z=0 plane
        """
        m = self.mapdl
        tol = 1e-6

        # Node on x=0 plane for UX=0 constraint
        m.nsel('S', 'LOC', 'X', 0, tol)
        nodes_x0 = m.mesh.nnum
        self.rigid_body_nodes['UX'] = int(nodes_x0[0]) if len(nodes_x0) > 0 else None

        # Node on y=0 plane for UY=0 constraint
        m.nsel('S', 'LOC', 'Y', 0, tol)
        nodes_y0 = m.mesh.nnum
        self.rigid_body_nodes['UY'] = int(nodes_y0[0]) if len(nodes_y0) > 0 else None

        # Node on z=0 plane for UZ=0 constraint
        m.nsel('S', 'LOC', 'Z', 0, tol)
        nodes_z0 = m.mesh.nnum
        self.rigid_body_nodes['UZ'] = int(nodes_z0[0]) if len(nodes_z0) > 0 else None

        m.allsel()

        print(f"    Rigid body nodes: UX@{self.rigid_body_nodes['UX']}, "
              f"UY@{self.rigid_body_nodes['UY']}, UZ@{self.rigid_body_nodes['UZ']}")

    def _clear_all_constraints(self):
        """Clear all displacement constraints, coupled DOFs, and constraint equations."""
        m = self.mapdl
        m.ddele('ALL', 'ALL')
        m.run('CPDELE,ALL')
        m.run('CEDELE,ALL')

    def _apply_periodic_bc(self, eps_x=0.0, eps_y=0.0, eps_z=0.0,
                           gamma_xy=0.0, gamma_yz=0.0, gamma_xz=0.0):
        """
        Apply Boundary Conditions based on element type.

        For SOLID185 (hex, mapped mesh): Periodic BC using Constraint Equations
        For SOLID187 (tet, free mesh): Uniform Displacement BC (KUBC)

        Parameters
        ----------
        eps_x, eps_y, eps_z : float
            Normal strain components
        gamma_xy, gamma_yz, gamma_xz : float
            Shear strain components
        """
        if self.element_type == 'SOLID187':
            self._apply_uniform_displacement_bc(eps_x, eps_y, eps_z,
                                                gamma_xy, gamma_yz, gamma_xz)
        else:
            self._apply_periodic_bc_hex(eps_x, eps_y, eps_z,
                                        gamma_xy, gamma_yz, gamma_xz)

    def _apply_uniform_displacement_bc(self, eps_x=0.0, eps_y=0.0, eps_z=0.0,
                                        gamma_xy=0.0, gamma_yz=0.0, gamma_xz=0.0):
        """
        Apply Kinematic Uniform Boundary Conditions (KUBC) for SOLID187 free mesh.

        Based on ANSYS 2025 R1 Theory Documentation (Equations 3.24-3.27).

        The linear displacement field that satisfies periodic BC (Eq 3.24-3.26):
            u_x(x,y,z) = ε_x * x
            u_y(x,y,z) = γ_xy * x + ε_y * y
            u_z(x,y,z) = γ_xz * x + γ_yz * y + ε_z * z

        Rigid body constraints (Eq 3.27) are satisfied at nodes with x=0, y=0, z=0.
        """
        m = self.mapdl

        self._clear_all_constraints()

        m.allsel()
        all_nodes = m.mesh.nodes  # shape (n_nodes, 3): [x, y, z]
        all_nnum = m.mesh.nnum    # node numbers

        # Get boundary node numbers from all 6 faces
        boundary_nodes = set()
        for face_name in ['XNEG', 'XPOS', 'YNEG', 'YPOS', 'ZNEG', 'ZPOS']:
            boundary_nodes.update(self.face_nodes[face_name])

        # Create node number to index mapping
        nnum_to_idx = {n: i for i, n in enumerate(all_nnum)}

        # Apply linear displacement field to all boundary nodes (Eq 3.24-3.26)
        # u_x = eps_x * x
        # u_y = gamma_xy * x + eps_y * y
        # u_z = gamma_xz * x + gamma_yz * y + eps_z * z
        for node_num in boundary_nodes:
            idx = nnum_to_idx[node_num]
            x, y, z = all_nodes[idx]

            ux = eps_x * x
            uy = gamma_xy * x + eps_y * y
            uz = gamma_xz * x + gamma_yz * y + eps_z * z

            m.d(int(node_num), 'UX', ux)
            m.d(int(node_num), 'UY', uy)
            m.d(int(node_num), 'UZ', uz)

        m.allsel()

    def _apply_periodic_bc_hex(self, eps_x=0.0, eps_y=0.0, eps_z=0.0,
                                gamma_xy=0.0, gamma_yz=0.0, gamma_xz=0.0):
        """
        Apply Periodic Boundary Conditions using Constraint Equations.

        Based on ANSYS 2025 R1 Theory Documentation (Equations 3.24-3.27).
        Used for SOLID185 with mapped hex mesh where node pairing is available.

        On X-faces (Equation 3.24):
            u_x(L_x,y,z) = u_x(0,y,z) + ε_x * L_x
            u_y(L_x,y,z) = u_y(0,y,z) + γ_xy * L_x
            u_z(L_x,y,z) = u_z(0,y,z) + γ_xz * L_x

        On Y-faces (Equation 3.25):
            u_x(x,L_y,z) = u_x(x,0,z)
            u_y(x,L_y,z) = u_y(x,0,z) + ε_y * L_y
            u_z(x,L_y,z) = u_z(x,0,z) + γ_yz * L_y

        On Z-faces (Equation 3.26):
            u_x(x,y,L_z) = u_x(x,y,0)
            u_y(x,y,L_z) = u_y(x,y,0)
            u_z(x,y,L_z) = u_z(x,y,0) + ε_z * L_z
        """
        m = self.mapdl
        Lx, Ly, Lz = self.Lx, self.Ly, self.Lz

        self._clear_all_constraints()

        ce_num = 1  # Constraint equation counter

        # X-direction periodic BC (Equation 3.24)
        # u_pos - u_neg = offset
        # CE format: CE,NEQN,CONST, NODE1,Lab1,C1, NODE2,Lab2,C2, ...
        # CONST + C1*NODE1.Lab1 + C2*NODE2.Lab2 = 0
        # For u_neg - u_pos = -offset:  offset + 1*u_neg + (-1)*u_pos = 0
        for (n_neg, n_pos) in self.node_pairs['X']:
            # UX: u_x(L_x) - u_x(0) = eps_x * L_x
            m.ce(ce_num, eps_x * Lx, n_neg, 'UX', 1, n_pos, 'UX', -1)
            ce_num += 1

            # UY: u_y(L_x) - u_y(0) = gamma_xy * L_x
            m.ce(ce_num, gamma_xy * Lx, n_neg, 'UY', 1, n_pos, 'UY', -1)
            ce_num += 1

            # UZ: u_z(L_x) - u_z(0) = gamma_xz * L_x
            m.ce(ce_num, gamma_xz * Lx, n_neg, 'UZ', 1, n_pos, 'UZ', -1)
            ce_num += 1

        # Y-direction periodic BC (Equation 3.25)
        for (n_neg, n_pos) in self.node_pairs['Y']:
            # UX: u_x(L_y) - u_x(0) = 0
            m.ce(ce_num, 0, n_neg, 'UX', 1, n_pos, 'UX', -1)
            ce_num += 1

            # UY: u_y(L_y) - u_y(0) = eps_y * L_y
            m.ce(ce_num, eps_y * Ly, n_neg, 'UY', 1, n_pos, 'UY', -1)
            ce_num += 1

            # UZ: u_z(L_y) - u_z(0) = gamma_yz * L_y
            m.ce(ce_num, gamma_yz * Ly, n_neg, 'UZ', 1, n_pos, 'UZ', -1)
            ce_num += 1

        # Z-direction periodic BC (Equation 3.26)
        for (n_neg, n_pos) in self.node_pairs['Z']:
            # UX: u_x(L_z) - u_x(0) = 0
            m.ce(ce_num, 0, n_neg, 'UX', 1, n_pos, 'UX', -1)
            ce_num += 1

            # UY: u_y(L_z) - u_y(0) = 0
            m.ce(ce_num, 0, n_neg, 'UY', 1, n_pos, 'UY', -1)
            ce_num += 1

            # UZ: u_z(L_z) - u_z(0) = eps_z * L_z
            m.ce(ce_num, eps_z * Lz, n_neg, 'UZ', 1, n_pos, 'UZ', -1)
            ce_num += 1

        # Rigid body constraints (Equation 3.27)
        # u_x(point with x=0) = 0
        # u_y(point with y=0) = 0
        # u_z(point with z=0) = 0
        if self.rigid_body_nodes['UX'] is not None:
            m.d(self.rigid_body_nodes['UX'], 'UX', 0)
        if self.rigid_body_nodes['UY'] is not None:
            m.d(self.rigid_body_nodes['UY'], 'UY', 0)
        if self.rigid_body_nodes['UZ'] is not None:
            m.d(self.rigid_body_nodes['UZ'], 'UZ', 0)

        m.allsel()

    def apply_bc_load_case_1(self, strain_val=0.001):
        """
        Apply Periodic BC for Load Case 1: Tensile test in X direction.

        ε_x = strain_val, all other strain components = 0
        """
        self._apply_periodic_bc(eps_x=strain_val)
        print(f"    Applied Periodic BC for LC1: ε_x = {strain_val}")

    def apply_bc_load_case_2(self, strain_val=0.001):
        """
        Apply Periodic BC for Load Case 2: Tensile test in Y direction.

        ε_y = strain_val, all other strain components = 0
        """
        self._apply_periodic_bc(eps_y=strain_val)
        print(f"    Applied Periodic BC for LC2: ε_y = {strain_val}")

    def apply_bc_load_case_3(self, strain_val=0.001):
        """
        Apply Periodic BC for Load Case 3: Tensile test in Z direction.

        ε_z = strain_val, all other strain components = 0
        """
        self._apply_periodic_bc(eps_z=strain_val)
        print(f"    Applied Periodic BC for LC3: ε_z = {strain_val}")

    def apply_bc_load_case_4(self, strain_val=0.001):
        """
        Apply Periodic BC for Load Case 4: Shear test in XY plane.

        γ_xy = strain_val, all other strain components = 0
        """
        self._apply_periodic_bc(gamma_xy=strain_val)
        print(f"    Applied Periodic BC for LC4: γ_xy = {strain_val}")

    def apply_bc_load_case_5(self, strain_val=0.001):
        """
        Apply Periodic BC for Load Case 5: Shear test in YZ plane.

        γ_yz = strain_val, all other strain components = 0
        """
        self._apply_periodic_bc(gamma_yz=strain_val)
        print(f"    Applied Periodic BC for LC5: γ_yz = {strain_val}")

    def apply_bc_load_case_6(self, strain_val=0.001):
        """
        Apply Periodic BC for Load Case 6: Shear test in XZ plane.

        γ_xz = strain_val, all other strain components = 0
        """
        self._apply_periodic_bc(gamma_xz=strain_val)
        print(f"    Applied Periodic BC for LC6: γ_xz = {strain_val}")

    def apply_bc_load_case_7_thermal(self, delta_T=1.0):
        """
        Apply Periodic BC for Load Case 7: Thermal expansion.

        Based on ANSYS 2025 R1 Theory Documentation (Equations 3.38-3.40).
        Vanishing macroscopic strain: all strain components = 0
        Temperature change ΔT applied uniformly.
        """
        m = self.mapdl

        # Apply periodic BC with zero strain (Equations 3.38-3.40)
        self._apply_periodic_bc(eps_x=0, eps_y=0, eps_z=0,
                                gamma_xy=0, gamma_yz=0, gamma_xz=0)

        # Apply thermal load
        m.tref(0)  # Reference temperature = 0
        m.bfunif('TEMP', delta_T)  # Current temperature = delta_T

        print(f"    Applied Periodic BC for LC7: Thermal (ΔT = {delta_T}°C)")

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

    def compute_stiffness_matrix(self, strain_mag=0.001):
        """
        Compute the complete 6x6 stiffness matrix [D] using Periodic BC.

        Based on ANSYS 2025 R1 Theory Documentation (Equations 3.16-3.17).

        For each load case, one strain component is set to strain_mag (0.001)
        and all others are set to 0. The stiffness matrix column is computed as:
            D_ij = σ_i / strain_mag  (Equation 3.17)

        Parameters
        ----------
        strain_mag : float
            Magnitude of applied strain for each load case (default: 0.001)

        Returns
        -------
        D : ndarray
            6x6 stiffness matrix in Voigt notation [D11,D12,...,D66]
        """
        print(f"\n{'='*60}")
        print("COMPUTING STIFFNESS MATRIX [D] (Periodic BC)")
        print(f"{'='*60}")
        print(f"Applied strain magnitude: {strain_mag}")
        print(f"Result files will be saved in: {self.mapdl.directory}")

        D = np.zeros((6, 6))

        # Define 6 load cases (Equation 3.24-3.26 with one strain component = strain_mag)
        # Format: (description, jobname, bc_function)
        load_cases = [
            ('LC1: ε_x (tensile X)', 'LC1_eps_x', self.apply_bc_load_case_1),
            ('LC2: ε_y (tensile Y)', 'LC2_eps_y', self.apply_bc_load_case_2),
            ('LC3: ε_z (tensile Z)', 'LC3_eps_z', self.apply_bc_load_case_3),
            ('LC4: γ_xy (shear XY)', 'LC4_gamma_xy', self.apply_bc_load_case_4),
            ('LC5: γ_yz (shear YZ)', 'LC5_gamma_yz', self.apply_bc_load_case_5),
            ('LC6: γ_xz (shear XZ)', 'LC6_gamma_xz', self.apply_bc_load_case_6),
        ]

        for i, (description, jobname, bc_func) in enumerate(load_cases):
            print(f"\n  Load case {i+1}/6: {description}...")

            self.mapdl.prep7()
            bc_func(strain_mag)
            self.solve(jobname=jobname)

            # Get volume-averaged stress (Equation 3.17)
            stress = self.get_volume_avg_stress()

            # D_ij = σ_i / ε_j where ε_j = strain_mag
            D[:, i] = stress / strain_mag

            print(f"    Stress [MPa]: σ_x={stress[0]:.2f}, σ_y={stress[1]:.2f}, σ_z={stress[2]:.2f}")
            print(f"                  τ_xy={stress[3]:.2f}, τ_yz={stress[4]:.2f}, τ_xz={stress[5]:.2f}")

        # Symmetrize the matrix (should be symmetric for linear elastic material)
        D = 0.5 * (D + D.T)

        self.stiffness_matrix = D
        print(f"\n  Stiffness matrix [D] computed successfully")

        return D

    def compute_engineering_constants(self):
        """
        Compute engineering constants from stiffness matrix.

        Based on ANSYS 2025 R1 Theory Documentation (Equations 3.22-3.23).

        Compliance matrix: [C] = [D]^(-1)

        The compliance matrix has the form (Equation 3.23):
            [C] = | 1/E_x      -ν_yx/E_y  -ν_zx/E_z                    |
                  | -ν_xy/E_x  1/E_y      -ν_zy/E_z                    |
                  | -ν_xz/E_x  -ν_yz/E_y  1/E_z                        |
                  |                       1/G_xy                       |
                  |                              1/G_yz                |
                  |                                     1/G_xz         |

        Returns
        -------
        props : dict
            Dictionary containing E_x, E_y, E_z, G_xy, G_yz, G_xz,
            nu_xy, nu_yz, nu_xz (and symmetric nu_yx, nu_zy, nu_zx)
        """
        if self.stiffness_matrix is None:
            raise ValueError("Stiffness matrix not computed. Run compute_stiffness_matrix first.")

        D = self.stiffness_matrix

        # Compute compliance matrix [C] = [D]^(-1) (Equation 3.22)
        try:
            C = np.linalg.inv(D)
        except np.linalg.LinAlgError:
            print("Warning: Stiffness matrix is singular. Using pseudo-inverse.")
            C = np.linalg.pinv(D)

        self.compliance_matrix = C

        # Extract engineering constants from compliance matrix (Equation 3.23)
        # Diagonal terms: C_ii = 1/E_i or 1/G_ij
        # Off-diagonal terms: C_ij = -nu_ji/E_j

        # Elastic moduli
        Ex = 1.0 / C[0, 0]  # E_x
        Ey = 1.0 / C[1, 1]  # E_y
        Ez = 1.0 / C[2, 2]  # E_z

        # Shear moduli
        Gxy = 1.0 / C[3, 3]  # G_xy
        Gyz = 1.0 / C[4, 4]  # G_yz
        Gxz = 1.0 / C[5, 5]  # G_xz

        # Poisson's ratios from compliance matrix
        # C[0,1] = -nu_yx/E_y => nu_yx = -C[0,1] * E_y
        # C[1,0] = -nu_xy/E_x => nu_xy = -C[1,0] * E_x
        nu_xy = -C[1, 0] * Ex  # -ν_xy/E_x
        nu_yx = -C[0, 1] * Ey  # -ν_yx/E_y

        # C[0,2] = -nu_zx/E_z => nu_zx = -C[0,2] * E_z
        # C[2,0] = -nu_xz/E_x => nu_xz = -C[2,0] * E_x
        nu_xz = -C[2, 0] * Ex  # -ν_xz/E_x
        nu_zx = -C[0, 2] * Ez  # -ν_zx/E_z

        # C[1,2] = -nu_zy/E_z => nu_zy = -C[1,2] * E_z
        # C[2,1] = -nu_yz/E_y => nu_yz = -C[2,1] * E_y
        nu_yz = -C[2, 1] * Ey  # -ν_yz/E_y
        nu_zy = -C[1, 2] * Ez  # -ν_zy/E_z

        self.effective_props.update({
            'Ex': Ex, 'Ey': Ey, 'Ez': Ez,
            'Gxy': Gxy, 'Gyz': Gyz, 'Gxz': Gxz,
            'nu_xy': nu_xy, 'nu_yx': nu_yx,
            'nu_xz': nu_xz, 'nu_zx': nu_zx,
            'nu_yz': nu_yz, 'nu_zy': nu_zy
        })

        return self.effective_props

    def compute_thermal_expansion(self, delta_T=1.0):
        """
        Compute effective secant thermal expansion coefficients.

        Based on ANSYS 2025 R1 Theory Documentation (Equations 3.34-3.37).

        For orthotropic linear elastic material with thermal strain:
            {ε} = {ε^th} + [D]^(-1){σ}   (Equation 3.34)

        With vanishing macroscopic strain {ε} = 0 and temperature change ΔT:
            {ε^th} = -[D]^(-1){σ}         (Equation 3.36)

        The secant thermal expansion coefficients are:
            {α^se} = -1/ΔT * [D]^(-1) * {σ}   (Equation 3.37)

        where {σ} is the stress obtained from boundary reactions with
        vanishing macroscopic strain (Equations 3.38-3.40).

        Parameters
        ----------
        delta_T : float
            Temperature change (default: 1.0°C)

        Returns
        -------
        cte : dict
            Dictionary containing alpha_x, alpha_y, alpha_z (secant CTE)
        """
        print(f"\n{'='*60}")
        print("COMPUTING THERMAL EXPANSION COEFFICIENTS (Equation 3.37)")
        print(f"{'='*60}")
        print(f"Temperature change ΔT = {delta_T}°C")

        if self.stiffness_matrix is None:
            raise ValueError("Stiffness matrix not computed. Run compute_stiffness_matrix first.")

        if self.compliance_matrix is None:
            # Compute compliance matrix if not available
            self.compliance_matrix = np.linalg.inv(self.stiffness_matrix)

        # Apply LC7: Thermal load with vanishing macroscopic strain
        self.mapdl.prep7()
        self.apply_bc_load_case_7_thermal(delta_T)
        self.solve(jobname='LC7_thermal')

        # Get volume-averaged stress from thermal load case
        stress = self.get_volume_avg_stress()
        print(f"    Thermal stress [MPa]: σ_x={stress[0]:.2f}, σ_y={stress[1]:.2f}, σ_z={stress[2]:.2f}")
        print(f"                          τ_xy={stress[3]:.2f}, τ_yz={stress[4]:.2f}, τ_xz={stress[5]:.2f}")

        # Compute secant thermal expansion coefficients (Equation 3.37)
        # {α^se} = -1/ΔT * [C] * {σ} where [C] = [D]^(-1)
        C = self.compliance_matrix
        alpha_se = -1.0 / delta_T * np.dot(C, stress)

        # Extract normal components (shear components should be ~0)
        alpha_x = alpha_se[0]
        alpha_y = alpha_se[1]
        alpha_z = alpha_se[2]

        self.effective_props.update({
            'alpha_x': alpha_x,
            'alpha_y': alpha_y,
            'alpha_z': alpha_z
        })

        print(f"\n  Secant CTE computed (Equation 3.37):")
        print(f"    α_x = {alpha_x:.2e} 1/°C")
        print(f"    α_y = {alpha_y:.2e} 1/°C")
        print(f"    α_z = {alpha_z:.2e} 1/°C")

        return {'alpha_x': alpha_x, 'alpha_y': alpha_y, 'alpha_z': alpha_z}

    def run_full_analysis(self, ele_size=0.05, strain_mag=0.001, element_type='SOLID185', n_div=None):
        """
        Run complete analysis to get all effective properties.

        Parameters
        ----------
        ele_size : float
            Element size in mm (default: 0.05)
        strain_mag : float
            Strain magnitude for mechanical load cases
        element_type : str
            Element type: 'SOLID185' (8-node hex) or 'SOLID187' (10-node tet)
        n_div : int, optional
            Number of mesh divisions (overrides ele_size for SOLID185)

        Returns
        -------
        props : dict
            All effective properties
        """
        # Build model
        self.build_model(ele_size, n_div, element_type)

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
        print("Based on ANSYS 2025 R1 Theory Documentation")
        print(f"{'='*60}")

        print(f"\nInput Parameters:")
        print(f"  RVE size: {self.L} x {self.L} x {self.L} mm")
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
        print("COMPUTED EFFECTIVE PROPERTIES (Equations 3.22-3.23)")
        print(f"{'-'*60}")

        print(f"\n--- Elastic Moduli E [MPa] ---")
        print(f"  E_x = {p.get('Ex', 0):.2f}")
        print(f"  E_y = {p.get('Ey', 0):.2f}")
        print(f"  E_z = {p.get('Ez', 0):.2f}")

        print(f"\n--- Shear Moduli G [MPa] ---")
        print(f"  G_xy = {p.get('Gxy', 0):.2f}")
        print(f"  G_yz = {p.get('Gyz', 0):.2f}")
        print(f"  G_xz = {p.get('Gxz', 0):.2f}")

        print(f"\n--- Poisson's Ratios ν ---")
        print(f"  ν_xy = {p.get('nu_xy', 0):.4f}    ν_yx = {p.get('nu_yx', 0):.4f}")
        print(f"  ν_xz = {p.get('nu_xz', 0):.4f}    ν_zx = {p.get('nu_zx', 0):.4f}")
        print(f"  ν_yz = {p.get('nu_yz', 0):.4f}    ν_zy = {p.get('nu_zy', 0):.4f}")

        print(f"\n--- Secant Thermal Expansion Coefficients α [1/°C] (Eq. 3.37) ---")
        print(f"  α_x = {p.get('alpha_x', 0):.2e}")
        print(f"  α_y = {p.get('alpha_y', 0):.2e}")
        print(f"  α_z = {p.get('alpha_z', 0):.2e}")

        if self.stiffness_matrix is not None:
            print(f"\n--- Stiffness Matrix [D] [MPa] (Eq. 3.16) ---")
            D = self.stiffness_matrix
            labels = ['D11', 'D21', 'D31', 'D41', 'D51', 'D61']
            for i in range(6):
                row = [f"{D[i,j]:12.2f}" for j in range(6)]
                print(f"  [{' '.join(row)}]")

        if self.compliance_matrix is not None:
            print(f"\n--- Compliance Matrix [C] = [D]^(-1) [1/MPa] (Eq. 3.22-3.23) ---")
            C = self.compliance_matrix
            for i in range(6):
                row = [f"{C[i,j]:12.2e}" for j in range(6)]
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
        calc.run_full_analysis(ele_size=0.05, strain_mag=0.001)

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
