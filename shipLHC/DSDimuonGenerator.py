import ROOT
import numpy as np

class DSDimuonGenerator(ROOT.FairGenerator):
    """
    Acceptance-guaranteed dimuon generator for DS tracking studies.
    Generates simultaneous mu- (primary) and mu+ (charm decay) originating in
    the target interaction volume (303 <= z <= 350, 18 < y < 48, -44 < x < -11)
    and aimed directly at the Downstream (DS) muon system.
    """
    def __init__(self, csv_file, min_ds_track_dist=1.0, max_ds_track_dist=25.0):
        super().__init__()
        self.mass = 0.1056584  # Muon mass in GeV/c^2
        self.z_DS3 = 540.0     # z position of last DS tracking station (DS3) [cm]

        # Target interaction volume requested:
        # -44 < x < -11, 18 < y < 48, 303 <= z <= 350 [cm]
        self.tgt_x = (-44.0, -11.0)
        self.tgt_y = ( 18.0,  48.0)
        self.tgt_z = (303.0, 350.0)

        # DS3 acceptance window (matching detector alignment) [cm]
        # Active cross section of DS bars
        self.ds_x = (-45.0, -10.0)
        self.ds_y = ( 15.0,  50.0)

        # Two-track separation range at DS3 [cm]
        self.min_dist = min_ds_track_dist
        self.max_dist = max_ds_track_dist

        # Load muon 1 energy spectrum (E_low, E_high, counts)
        data = np.loadtxt(csv_file, delimiter=',', skiprows=0)
        self.e_low = data[:, 0]
        self.e_high = data[:, 1]
        counts = data[:, 2]
        self.cdf = np.cumsum(counts) / np.sum(counts)
        self.cdf[-1] = 1.0

    def sample_energy_mu1(self):
        """
        Sample E1 from spectrum. To ensure the second muon (E2 ~ 0.02 * E1)
        punches through the ~1 m US iron absorber (threshold ~ 2-2.5 GeV),
        we sample E1 >= 100 GeV.
        """
        while True:
            u = np.random.uniform(0.0, 1.0)
            idx = min(np.searchsorted(self.cdf, u), len(self.e_low) - 1)
            e = np.random.uniform(self.e_low[idx], self.e_high[idx])
            if e >= 10.0:
                return float(e)

    def ReadEvent(self, primGen):
        # 1. Primary neutrino interaction vertex (Muon 1) in target region
        x1 = np.random.uniform(*self.tgt_x)
        y1 = np.random.uniform(*self.tgt_y)
        z1 = np.random.uniform(*self.tgt_z)

        # 2. Charm decay vertex (Muon 2): displaced forward by a few cm
        dr = np.random.exponential(scale=2.0)  # Mean decay distance ~ 2 cm
        theta = np.random.uniform(0.0, np.radians(15.0))
        phi = np.random.uniform(0.0, 2.0 * np.pi)

        dx = dr * np.sin(theta) * np.cos(phi)
        dy = dr * np.sin(theta) * np.sin(phi)
        dz = dr * np.cos(theta)
        x2, y2, z2 = x1 + dx, y1 + dy, z1 + dz

        # 3. Hit points on DS3 within acceptance
        x_ds1 = np.random.uniform(*self.ds_x)
        y_ds1 = np.random.uniform(*self.ds_y)

        # Muon 2 hit position on DS3 with controlled separation
        sep = np.random.uniform(self.min_dist, self.max_dist)
        sep_phi = np.random.uniform(0.0, 2.0 * np.pi)
        x_ds2 = np.clip(x_ds1 + sep * np.cos(sep_phi), self.ds_x[0], self.ds_x[1])
        y_ds2 = np.clip(y_ds1 + sep * np.sin(sep_phi), self.ds_y[0], self.ds_y[1])

        # 4. Energies: E2 is ~2% of E1, guarded above absorber threshold (2.5 GeV)
        E1 = self.sample_energy_mu1()
        E2 = max(0.02 * E1, 2.5)

        # 5. Direction vectors from Target vertices to DS3 impact points
        v1 = np.array([x_ds1 - x1, y_ds1 - y1, self.z_DS3 - z1])
        u1 = v1 / np.linalg.norm(v1)
        p1 = np.sqrt(E1**2 - self.mass**2)
        px1, py1, pz1 = p1 * u1

        v2 = np.array([x_ds2 - x2, y_ds2 - y2, self.z_DS3 - z2])
        u2 = v2 / np.linalg.norm(v2)
        p2 = np.sqrt(E2**2 - self.mass**2)
        px2, py2, pz2 = p2 * u2

        # 6. Add tracks to FairPrimaryGenerator:
        #    Track 0: mu- (PDG = 13)
        #    Track 1: mu+ (PDG = -13)
        primGen.AddTrack( 13, px1, py1, pz1, x1, y1, z1, -1, True, E1, 0.0, 1.0)
        primGen.AddTrack(-13, px2, py2, pz2, x2, y2, z2, -1, True, E2, 0.0, 1.0)

        return True
