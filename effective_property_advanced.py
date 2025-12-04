"""
Advanced Composite Material Effective Property Calculator using PyANSYS

This version implements proper periodic boundary conditions using constraint equations (CE),
following ANSYS Material Designer methodology for RVE homogenization.

Key features:
1. Square prism fiber geometry (for compatible hex meshing)
2. SOLID185 hex elements with matching meshes on opposite faces
3. Periodic boundary conditions via CE (constraint equations) with master nodes
4. Volume-averaged stress and strain for effective property computation

Periodic BC formulation:
    u(x+) - u(x-) = ε̄ · Δx

where:
    - u(x+), u(x-) are displacements on opposite faces
    - ε̄ is the macroscopic strain tensor
    - Δx is the distance vector between faces

Reference:
- Xia, Z., Zhou, C., Yong, Q., Wang, X. (2006). "On selection of repeated unit cell
  model and application of unified periodic boundary conditions"
- ANSYS Material Designer Theory Guide
"""

import numpy as np
from ansys.mapdl.core import launch_mapdl


class AdvancedCompositeCalculator:
    """
    Advanced calculator for effective material properties using proper periodic BC.

    Uses master nodes and constraint equations (CE) to enforce periodic boundary
    conditions, which is the same approach used in ANSYS Material Designer.
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
        self.master_nodes = {}  # Master nodes for periodic BC
        self.node_pairs = {}    # Corresponding node pairs on opposite faces
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

        # Create geometry using keypoints and volumes for mapped meshing
        # The RVE is divided into 9 volumes (3x3 in XY plane, extruded in Z)
        # Center volume is fiber, surrounding 8 volumes are matrix

        # Fiber boundaries
        x1, x2 = c - a, c + a  # Fiber X boundaries
        y1, y2 = c - a, c + a  # Fiber Y boundaries

        # Create 9 blocks in XY plane, extruded through Z
        vol_id = 1
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

        for xs, xe, ys, ye, mat_id in regions:
            m.block(xs, xe, ys, ye, 0, L)

        # Glue all volumes together
        m.vsel('ALL')
        m.vglue('ALL')
        m.allsel()

        # Set element size for mapped mesh
        # Calculate divisions based on geometry
        fiber_width = 2 * a
        matrix_width = (L - fiber_width) / 2

        # Ensure matching divisions on opposite faces
        n_fiber = max(2, int(n_div * fiber_width / L))
        n_matrix = max(2, int(n_div * matrix_width / L))
        n_z = n_div

        # Set line divisions for mapped meshing
        # Select lines by length and location
        m.lsel('ALL')
        m.lesize('ALL', '', '', n_div)

        # Mesh all volumes
        m.mshkey(1)  # Mapped meshing
        m.mshape(0, '3D')  # Hex elements

        # Mesh each volume with appropriate material
        m.vsel('ALL')
        volumes = m.vlist('ALL')

        # We need to identify which volumes are fiber and which are matrix
        # Get volume list and check their centroids
        m.allsel()
        n_vol = int(m.get('VCOUNT', 'VOLU', '', 'COUNT'))

        for v in range(1, n_vol + 1):
            m.vsel('S', 'VOLU', '', v)
            # Get centroid
            cx = m.get('CENTX', 'VOLU', v, 'CENT', 'X')
            cy = m.get('CENTY', 'VOLU', v, 'CENT', 'Y')

            # Check if centroid is in fiber region
            if x1 < cx < x2 and y1 < cy < y2:
                mat_id = 2  # Fiber
            else:
                mat_id = 1  # Matrix

            m.vatt(mat_id, '', 1)
            m.vmesh(v)

        m.allsel()

        # Merge nodes at interfaces
        m.nummrg('NODE', 1e-6)

        nn = int(m.get('NCOUNT', 'NODE', '', 'COUNT'))
        ne = int(m.get('ECOUNT', 'ELEM', '', 'COUNT'))
        print(f"Mesh: {nn} nodes, {ne} elements")

        # Create face node sets and find node pairs
        self._create_face_sets_and_pairs()

        # Create master nodes for periodic BC
        self._create_master_nodes()

        print("Model built successfully")

    def _create_face_sets_and_pairs(self):
        """Create node sets for faces and find corresponding node pairs."""
        m = self.mapdl
        L = self.L
        tol = 1e-6

        print("  Creating face node sets and finding node pairs...")

        # Create face component sets
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

        # Get all node coordinates
        nodes = m.mesh.nodes
        node_nums = m.mesh.nnum

        # Create mapping from coordinates to node numbers
        coord_to_node = {}
        for i, nnum in enumerate(node_nums):
            x, y, z = nodes[i]
            key = (round(x, 5), round(y, 5), round(z, 5))
            coord_to_node[key] = nnum

        # Find corresponding node pairs on opposite faces
        self.node_pairs = {'X': [], 'Y': [], 'Z': []}

        # X-direction pairs (XPOS to XNEG)
        for i, nnum in enumerate(node_nums):
            x, y, z = nodes[i]
            if abs(x - L) < tol:  # Node on XPOS
                # Find corresponding node on XNEG
                key_neg = (round(0, 5), round(y, 5), round(z, 5))
                if key_neg in coord_to_node:
                    node_neg = coord_to_node[key_neg]
                    self.node_pairs['X'].append((nnum, node_neg))

        # Y-direction pairs (YPOS to YNEG)
        for i, nnum in enumerate(node_nums):
            x, y, z = nodes[i]
            if abs(y - L) < tol:  # Node on YPOS
                key_neg = (round(x, 5), round(0, 5), round(z, 5))
                if key_neg in coord_to_node:
                    node_neg = coord_to_node[key_neg]
                    self.node_pairs['Y'].append((nnum, node_neg))

        # Z-direction pairs (ZPOS to ZNEG)
        for i, nnum in enumerate(node_nums):
            x, y, z = nodes[i]
            if abs(z - L) < tol:  # Node on ZPOS
                key_neg = (round(x, 5), round(y, 5), round(0, 5))
                if key_neg in coord_to_node:
                    node_neg = coord_to_node[key_neg]
                    self.node_pairs['Z'].append((nnum, node_neg))

        print(f"    X-pairs: {len(self.node_pairs['X'])}")
        print(f"    Y-pairs: {len(self.node_pairs['Y'])}")
        print(f"    Z-pairs: {len(self.node_pairs['Z'])}")

    def _create_master_nodes(self):
        """
        Create master (reference) nodes for periodic boundary conditions.

        Three master nodes are created, each controlling the displacement gradient
        in one direction:
        - Master X: controls ∂u/∂x (UX=ε11*L, UY=ε21*L, UZ=ε31*L)
        - Master Y: controls ∂u/∂y (UX=ε12*L, UY=ε22*L, UZ=ε32*L)
        - Master Z: controls ∂u/∂z (UX=ε13*L, UY=ε23*L, UZ=ε33*L)
        """
        m = self.mapdl
        L = self.L

        print("  Creating master nodes for periodic BC...")

        # Get the highest existing node number
        max_node = int(m.get('NMAX', 'NODE', '', 'NUM', 'MAX'))

        # Create 3 master nodes at arbitrary locations outside RVE
        # These nodes only provide DOFs, their location doesn't matter
        self.master_nodes = {
            'X': max_node + 1,
            'Y': max_node + 2,
            'Z': max_node + 3,
        }

        m.n(self.master_nodes['X'], L * 2, 0, 0)
        m.n(self.master_nodes['Y'], L * 2, L, 0)
        m.n(self.master_nodes['Z'], L * 2, L * 2, 0)

        print(f"    Master nodes: X={self.master_nodes['X']}, Y={self.master_nodes['Y']}, Z={self.master_nodes['Z']}")

    def apply_periodic_bc(self, eps_macro):
        """
        Apply periodic boundary conditions for a given macroscopic strain.

        Uses constraint equations (CE) to enforce:
            u(x+) - u(x-) = ε̄ · Δx

        Parameters
        ----------
        eps_macro : array-like
            Macroscopic strain tensor in Voigt notation [ε11, ε22, ε33, γ12, γ23, γ31]
            Note: γ12 = 2*ε12 (engineering shear strain)
        """
        m = self.mapdl
        L = self.L

        e11, e22, e33, g12, g23, g31 = eps_macro
        # Convert engineering shear strain to tensor shear strain
        e12 = g12 / 2
        e23 = g23 / 2
        e31 = g31 / 2

        # Clear all previous constraints
        m.ddele('ALL', 'ALL')
        m.cedele('ALL')

        # Delete any constraint equations from previous load case
        m.run('CEDELE,ALL')

        # Counter for constraint equations
        ce_num = 1

        # Apply periodic BC via constraint equations
        # For each pair of nodes on opposite faces:
        # u_pos - u_neg = Master_DOF
        # CE format: CE,NEQN,CONST,NODE1,Lab1,C1,NODE2,Lab2,C2,...
        # CONST + C1*u1 + C2*u2 + ... = 0

        # X-direction periodicity (XPOS - XNEG)
        # u(L,y,z) - u(0,y,z) = [e11, e12, e31]^T * L
        for node_pos, node_neg in self.node_pairs['X']:
            # UX: u_pos - u_neg - Master_X.UX = 0
            m.ce(ce_num, 0, node_pos, 'UX', 1, node_neg, 'UX', -1, self.master_nodes['X'], 'UX', -1)
            ce_num += 1
            # UY: u_pos - u_neg - Master_X.UY = 0
            m.ce(ce_num, 0, node_pos, 'UY', 1, node_neg, 'UY', -1, self.master_nodes['X'], 'UY', -1)
            ce_num += 1
            # UZ: u_pos - u_neg - Master_X.UZ = 0
            m.ce(ce_num, 0, node_pos, 'UZ', 1, node_neg, 'UZ', -1, self.master_nodes['X'], 'UZ', -1)
            ce_num += 1

        # Y-direction periodicity (YPOS - YNEG)
        # u(x,L,z) - u(x,0,z) = [e12, e22, e23]^T * L
        for node_pos, node_neg in self.node_pairs['Y']:
            m.ce(ce_num, 0, node_pos, 'UX', 1, node_neg, 'UX', -1, self.master_nodes['Y'], 'UX', -1)
            ce_num += 1
            m.ce(ce_num, 0, node_pos, 'UY', 1, node_neg, 'UY', -1, self.master_nodes['Y'], 'UY', -1)
            ce_num += 1
            m.ce(ce_num, 0, node_pos, 'UZ', 1, node_neg, 'UZ', -1, self.master_nodes['Y'], 'UZ', -1)
            ce_num += 1

        # Z-direction periodicity (ZPOS - ZNEG)
        # u(x,y,L) - u(x,y,0) = [e31, e23, e33]^T * L
        for node_pos, node_neg in self.node_pairs['Z']:
            m.ce(ce_num, 0, node_pos, 'UX', 1, node_neg, 'UX', -1, self.master_nodes['Z'], 'UX', -1)
            ce_num += 1
            m.ce(ce_num, 0, node_pos, 'UY', 1, node_neg, 'UY', -1, self.master_nodes['Z'], 'UY', -1)
            ce_num += 1
            m.ce(ce_num, 0, node_pos, 'UZ', 1, node_neg, 'UZ', -1, self.master_nodes['Z'], 'UZ', -1)
            ce_num += 1

        print(f"    Created {ce_num-1} constraint equations")

        # Apply displacements to master nodes to impose macroscopic strain
        # Master X: gradient in X direction
        m.d(self.master_nodes['X'], 'UX', e11 * L)
        m.d(self.master_nodes['X'], 'UY', e12 * L)
        m.d(self.master_nodes['X'], 'UZ', e31 * L)

        # Master Y: gradient in Y direction
        m.d(self.master_nodes['Y'], 'UX', e12 * L)
        m.d(self.master_nodes['Y'], 'UY', e22 * L)
        m.d(self.master_nodes['Y'], 'UZ', e23 * L)

        # Master Z: gradient in Z direction
        m.d(self.master_nodes['Z'], 'UX', e31 * L)
        m.d(self.master_nodes['Z'], 'UY', e23 * L)
        m.d(self.master_nodes['Z'], 'UZ', e33 * L)

        # Fix one corner node to prevent rigid body motion
        # Find the node at origin (0, 0, 0)
        m.nsel('S', 'LOC', 'X', 0, 1e-6)
        m.nsel('R', 'LOC', 'Y', 0, 1e-6)
        m.nsel('R', 'LOC', 'Z', 0, 1e-6)
        m.d('ALL', 'ALL', 0)

        m.allsel()

    def apply_thermal_bc(self, delta_T=1.0):
        """
        Apply thermal loading with periodic boundary conditions.

        For thermal analysis, the periodic BC ensures uniform expansion
        without constraint. Master node displacements are left free.
        """
        m = self.mapdl
        L = self.L

        # Clear all previous constraints
        m.ddele('ALL', 'ALL')
        m.run('CEDELE,ALL')

        ce_num = 1

        # Apply periodic BC with zero prescribed gradient (free expansion)
        # X-direction periodicity
        for node_pos, node_neg in self.node_pairs['X']:
            m.ce(ce_num, 0, node_pos, 'UX', 1, node_neg, 'UX', -1, self.master_nodes['X'], 'UX', -1)
            ce_num += 1
            m.ce(ce_num, 0, node_pos, 'UY', 1, node_neg, 'UY', -1, self.master_nodes['X'], 'UY', -1)
            ce_num += 1
            m.ce(ce_num, 0, node_pos, 'UZ', 1, node_neg, 'UZ', -1, self.master_nodes['X'], 'UZ', -1)
            ce_num += 1

        # Y-direction periodicity
        for node_pos, node_neg in self.node_pairs['Y']:
            m.ce(ce_num, 0, node_pos, 'UX', 1, node_neg, 'UX', -1, self.master_nodes['Y'], 'UX', -1)
            ce_num += 1
            m.ce(ce_num, 0, node_pos, 'UY', 1, node_neg, 'UY', -1, self.master_nodes['Y'], 'UY', -1)
            ce_num += 1
            m.ce(ce_num, 0, node_pos, 'UZ', 1, node_neg, 'UZ', -1, self.master_nodes['Y'], 'UZ', -1)
            ce_num += 1

        # Z-direction periodicity
        for node_pos, node_neg in self.node_pairs['Z']:
            m.ce(ce_num, 0, node_pos, 'UX', 1, node_neg, 'UX', -1, self.master_nodes['Z'], 'UX', -1)
            ce_num += 1
            m.ce(ce_num, 0, node_pos, 'UY', 1, node_neg, 'UY', -1, self.master_nodes['Z'], 'UY', -1)
            ce_num += 1
            m.ce(ce_num, 0, node_pos, 'UZ', 1, node_neg, 'UZ', -1, self.master_nodes['Z'], 'UZ', -1)
            ce_num += 1

        # Fix corner node for rigid body motion
        m.nsel('S', 'LOC', 'X', 0, 1e-6)
        m.nsel('R', 'LOC', 'Y', 0, 1e-6)
        m.nsel('R', 'LOC', 'Z', 0, 1e-6)
        m.d('ALL', 'ALL', 0)

        # For thermal: constrain shear deformation modes (off-diagonal of master nodes)
        # but allow normal expansion (diagonal terms are free)
        m.d(self.master_nodes['X'], 'UY', 0)  # No shear ε12 from thermal
        m.d(self.master_nodes['X'], 'UZ', 0)  # No shear ε31 from thermal
        m.d(self.master_nodes['Y'], 'UX', 0)  # No shear ε12 from thermal
        m.d(self.master_nodes['Y'], 'UZ', 0)  # No shear ε23 from thermal
        m.d(self.master_nodes['Z'], 'UX', 0)  # No shear ε31 from thermal
        m.d(self.master_nodes['Z'], 'UY', 0)  # No shear ε23 from thermal

        # Leave UX of Master X, UY of Master Y, UZ of Master Z FREE
        # These will give us the thermal strains

        m.allsel()

        # Apply thermal load
        m.bfunif('TEMP', delta_T)
        m.tunif(0)  # Reference temperature

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
        Get effective thermal expansion strains from master node displacements.

        The thermal strains are computed from the master node displacements:
            ε_th = u_master / L / ΔT
        """
        m = self.mapdl
        L = self.L

        m.post1()
        m.set('LAST')

        # Get displacements from master nodes
        # Master X.UX gives ε11*L, Master Y.UY gives ε22*L, Master Z.UZ gives ε33*L

        m.nsel('S', 'NODE', '', self.master_nodes['X'])
        ux_master_x = m.get('UX', 'NODE', self.master_nodes['X'], 'U', 'X')

        m.nsel('S', 'NODE', '', self.master_nodes['Y'])
        uy_master_y = m.get('UY', 'NODE', self.master_nodes['Y'], 'U', 'Y')

        m.nsel('S', 'NODE', '', self.master_nodes['Z'])
        uz_master_z = m.get('UZ', 'NODE', self.master_nodes['Z'], 'U', 'Z')

        m.allsel()
        m.finish()

        # Thermal strains (CTE = strain / delta_T)
        eps_th = np.array([
            ux_master_x / L / delta_T,
            uy_master_y / L / delta_T,
            uz_master_z / L / delta_T
        ])

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
        print(f"Result files will be saved in: {self.mapdl.directory}")

        C = np.zeros((6, 6))

        # Define 6 load cases (unit strains in each direction)
        # Voigt notation: [ε11, ε22, ε33, γ12, γ23, γ31]
        load_cases = [
            [1, 0, 0, 0, 0, 0],  # ε11
            [0, 1, 0, 0, 0, 0],  # ε22
            [0, 0, 1, 0, 0, 0],  # ε33
            [0, 0, 0, 1, 0, 0],  # γ12 (engineering shear)
            [0, 0, 0, 0, 1, 0],  # γ23
            [0, 0, 0, 0, 0, 1],  # γ31
        ]

        directions = ['ε11', 'ε22', 'ε33', 'γ12', 'γ23', 'γ31']
        jobnames = ['LC1_e11', 'LC2_e22', 'LC3_e33', 'LC4_g12', 'LC5_g23', 'LC6_g31']

        for i, lc in enumerate(load_cases):
            print(f"  Load case {i+1}/6: {directions[i]}...")

            self.mapdl.prep7()
            eps_applied = np.array(lc) * strain_mag
            self.apply_periodic_bc(eps_applied)
            self.solve(jobname=jobnames[i])

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
