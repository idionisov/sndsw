import os
import ROOT
import numpy as np
import pandas as pd

class DSDimuonGenerator(ROOT.FairGenerator):
    """
    Acceptance-guaranteed dimuon generator for DS tracking studies.
    Generates two muons from one common vertex randomly distributed in the
    target region (303 <= z <= 350, 18 < y < 48, -44 < x < -11 cm), both
    directed to traverse the Downstream (DS) system.

    Muon 1 energy is sampled according to direct_mu1_energy.csv.
    Muon 2 energy is set by sampling the relative energy ratio (E2/E1)
    from direct_dimuon_e_ratio.csv.
    """
    def __init__(
        self,
        spectrum_csv="/afs/cern.ch/work/i/idioniso/snd-numu-charm/pgun/direct_mu1_energy.csv",
        ratio_csv="/afs/cern.ch/work/i/idioniso/snd-numu-charm/pgun/direct_dimuon_e_ratio.csv",
        min_ds_track_dist=1.0,
        max_ds_track_dist=25.0,
        tgt_x=(-44.0, -11.0),
        tgt_y=(18.0, 48.0),
        tgt_z=(303.0, 350.0),
        ds_x=(-45.0, -10.0),
        ds_y=(15.0, 50.0),
        z_DS3=540.0,
        pdg1=13,
        pdg2=-13,
        csv_file=None,
        **kwargs,
    ):
        super().__init__()
        if csv_file is not None:
            spectrum_csv = csv_file

        self.mass = 0.1056584  # Muon mass in GeV/c^2
        self.pdg1 = pdg1
        self.pdg2 = pdg2
        self.z_DS3 = z_DS3

        # Target interaction volume:
        self.tgt_x = tgt_x
        self.tgt_y = tgt_y
        self.tgt_z = tgt_z

        # DS3 acceptance window (matching detector alignment) [cm]
        self.ds_x = ds_x
        self.ds_y = ds_y

        # Separation at DS3 [cm]
        self.min_dist = min_ds_track_dist
        self.max_dist = max_ds_track_dist

        # Load muon 1 energy spectrum and ratio CDFs
        self.e1_low, self.e1_high, self.e1_cdf = self._load_cdf(spectrum_csv)
        self.ratio_low, self.ratio_high, self.ratio_cdf = self._load_cdf(ratio_csv)

    @staticmethod
    def _load_cdf(csv_file):
        if not os.path.exists(csv_file):
            raise FileNotFoundError(f"CSV file not found: {csv_file}")

        df = pd.read_csv(csv_file)
        if 'bin_low' in df.columns and 'bin_high' in df.columns:
            low = df['bin_low'].to_numpy(dtype=float)
            high = df['bin_high'].to_numpy(dtype=float)
            if 'cumulative_prob' in df.columns:
                cdf = df['cumulative_prob'].to_numpy(dtype=float)
            elif 'probability' in df.columns:
                prob = df['probability'].to_numpy(dtype=float)
                cdf = np.cumsum(prob) / np.sum(prob)
            elif 'content' in df.columns:
                content = df['content'].to_numpy(dtype=float)
                cdf = np.cumsum(content) / np.sum(content)
            else:
                raise ValueError(f"No probability column in {csv_file}")
        else:
            # Fallback for 3-column format without header: low, high, counts
            data = np.loadtxt(csv_file, delimiter=',')
            low = data[:, 0]
            high = data[:, 1]
            counts = data[:, 2]
            cdf = np.cumsum(counts) / np.sum(counts)

        cdf = cdf / cdf[-1]
        cdf[-1] = 1.0
        return low, high, cdf

    def sample_energies(self):
        """Sample E1 from spectrum and E2 = E1 * ratio from ratio distribution."""
        while True:
            u1 = np.random.uniform(0.0, 1.0)
            idx1 = min(np.searchsorted(self.e1_cdf, u1), len(self.e1_cdf) - 1)
            e1 = float(np.random.uniform(self.e1_low[idx1], self.e1_high[idx1]))
            if e1 <= self.mass:
                continue

            u2 = np.random.uniform(0.0, 1.0)
            idx2 = min(np.searchsorted(self.ratio_cdf, u2), len(self.ratio_cdf) - 1)
            ratio = float(np.random.uniform(self.ratio_low[idx2], self.ratio_high[idx2]))
            e2 = e1 * ratio
            if e2 > self.mass:
                return e1, e2

    def ReadEvent(self, primGen):
        # 1. Single vertex in target region for both muons
        x = np.random.uniform(*self.tgt_x)
        y = np.random.uniform(*self.tgt_y)
        z = np.random.uniform(*self.tgt_z)

        # 2. Energies from input distributions
        E1, E2 = self.sample_energies()

        # 3. Hit points on DS3 within active acceptance
        x_ds1 = np.random.uniform(*self.ds_x)
        y_ds1 = np.random.uniform(*self.ds_y)

        # Muon 2 hit position on DS3
        if self.min_dist is not None and self.max_dist is not None and self.max_dist > 0:
            found = False
            for _ in range(100):
                sep = np.random.uniform(self.min_dist, self.max_dist)
                sep_phi = np.random.uniform(0.0, 2.0 * np.pi)
                x_cand = x_ds1 + sep * np.cos(sep_phi)
                y_cand = y_ds1 + sep * np.sin(sep_phi)
                if self.ds_x[0] <= x_cand <= self.ds_x[1] and self.ds_y[0] <= y_cand <= self.ds_y[1]:
                    x_ds2, y_ds2 = x_cand, y_cand
                    found = True
                    break
            if not found:
                x_ds2 = np.clip(x_cand, self.ds_x[0], self.ds_x[1])
                y_ds2 = np.clip(y_cand, self.ds_y[0], self.ds_y[1])
        else:
            x_ds2 = np.random.uniform(*self.ds_x)
            y_ds2 = np.random.uniform(*self.ds_y)

        # 4. Momentum vectors from common vertex to DS3 points
        v1 = np.array([x_ds1 - x, y_ds1 - y, self.z_DS3 - z], dtype=float)
        u1 = v1 / np.linalg.norm(v1)
        p1 = np.sqrt(E1**2 - self.mass**2)
        px1, py1, pz1 = p1 * u1

        v2 = np.array([x_ds2 - x, y_ds2 - y, self.z_DS3 - z], dtype=float)
        u2 = v2 / np.linalg.norm(v2)
        p2 = np.sqrt(E2**2 - self.mass**2)
        px2, py2, pz2 = p2 * u2

        # 5. Add both tracks originating from the same vertex (x, y, z)
        primGen.AddTrack(self.pdg1, px1, py1, pz1, x, y, z, -1, True, E1, 0.0, 1.0)
        primGen.AddTrack(self.pdg2, px2, py2, pz2, x, y, z, -1, True, E2, 0.0, 1.0)

        return True
