import xml.etree.ElementTree as ET
from array import array

import matplotlib.pyplot as plt
import numpy as np
import ROOT
import scipy.ndimage
import shipunit as unit


def hit_finder(slope, intercept, box_centers, box_ds, tol = 0.) :
    """ Finds hits intersected by Hough line """
    d = np.abs(box_centers[0,:,1] - (box_centers[0,:,0]*slope + intercept))
    center_in_box = d < (box_ds[0,:,1]+tol)/2.
    clips_corner = np.abs(slope) > np.abs((d - (box_ds[0,:,1]+tol)/2.)/(box_ds[0,:,0]+tol)/2.)
    hit_mask = np.logical_or(center_in_box, clips_corner)
    return np.where(hit_mask)[0]

class hough() :
    """ Hough transform implementation """
    def __init__(self, n_yH, yH_range, n_xH, xH_range, z_offset, Hformat, space_scale, det_Zlen, squaretheta = False, smooth = True) :
        self.n_yH, self.n_xH = n_yH, n_xH
        self.yH_range, self.xH_range = yH_range, xH_range
        self.z_offset, self.HoughSpace_format, self.space_scale = z_offset, Hformat, space_scale
        self.det_Zlen, self.smooth = det_Zlen, smooth
        self.yH_bins = np.linspace(self.yH_range[0], self.yH_range[1], n_yH)

        if not squaretheta :
            self.xH_bins = np.linspace(self.xH_range[0], self.xH_range[1], n_xH)
        else :
            self.xH_bins = np.linspace(np.sign(self.xH_range[0])*(self.xH_range[0]**0.5), np.sign(self.xH_range[1])*(self.xH_range[1]**0.5), n_xH)
            self.xH_bins = np.sign(self.xH_bins)*np.square(self.xH_bins)

        self.cos_thetas, self.sin_thetas = np.cos(self.xH_bins), np.sin(self.xH_bins)
        self.xH_i = np.array(list(range(n_xH)))
        self.n_yH_scaled, self.n_xH_scaled = int(n_yH*space_scale), int(n_xH*space_scale)
        self.yH_bins_scaled = np.linspace(self.yH_range[0], self.yH_range[1], self.n_yH_scaled)

        if not squaretheta :
            self.xH_bins_scaled= np.linspace(
                self.xH_range[0],
                self.xH_range[1],
                self.n_xH_scaled
            )
        else :
            self.xH_bins_scaled = np.linspace(
                np.sign(self.xH_range[0])*(self.xH_range[0]**0.5),
                np.sign(self.xH_range[1])*(self.xH_range[1]**0.5),
                self.n_xH_scaled
            )
            self.xH_bins_scaled = np.sign(self.xH_bins_scaled)*np.square(self.xH_bins_scaled)

        self.cos_thetas_scaled, self.sin_thetas_scaled = np.cos(self.xH_bins_scaled), np.sin(self.xH_bins_scaled)
        self.xH_i_scaled = np.array(list(range(self.n_xH_scaled)))

    def fit(self, hit_collection, is_scaled, draw, weights = None) :
        if not is_scaled:
           n_xH, n_yH, cos_thetas, sin_thetas, xH_bins, yH_bins, xH_i, res = self.n_xH, self.n_yH, self.cos_thetas, self.sin_thetas, self.xH_bins, self.yH_bins, self.xH_i, self.res
        else:
           n_xH, n_yH, cos_thetas, sin_thetas, xH_bins, yH_bins, xH_i, res = self.n_xH_scaled, self.n_yH_scaled, self.cos_thetas_scaled, self.sin_thetas_scaled, self.xH_bins_scaled, self.yH_bins_scaled, self.xH_i_scaled, self.res*self.space_scale
        self.accumulator = np.zeros((n_yH, n_xH))
        for i_hit, hit in enumerate(hit_collection) :
            shifted_hitZ = hit[0] - self.z_offset
            if self.HoughSpace_format == 'normal': hit_yH = shifted_hitZ*cos_thetas + hit[1]*sin_thetas
            elif self.HoughSpace_format == 'linearSlopeIntercept': hit_yH = hit[1] - shifted_hitZ*xH_bins
            elif self.HoughSpace_format== 'linearIntercepts': hit_yH = (self.det_Zlen*hit[1] - shifted_hitZ*xH_bins)/(self.det_Zlen - shifted_hitZ)
            out_of_range = np.logical_and(hit_yH > self.yH_range[0], hit_yH < self.yH_range[1])
            hit_yH_i = np.floor((hit_yH[out_of_range] - self.yH_range[0])/(self.yH_range[1] - self.yH_range[0])*n_yH).astype(np.int_)
            if weights is not None : self.accumulator[hit_yH_i, xH_i[out_of_range]] += weights[i_hit]
            else : self.accumulator[hit_yH_i, xH_i[out_of_range]] += 1
        if self.smooth_full : self.accumulator = scipy.ndimage.gaussian_filter(self.accumulator, self.sigma, truncate=self.truncate)
        if draw :
            plt.figure(); plt.imshow(self.accumulator, origin = "lower", extent = [self.xH_range[0], self.xH_range[-1], self.yH_range[0], self.yH_range[-1]], aspect = "auto"); plt.tight_layout(); plt.show()
        if self.smooth_full: i_max = np.unravel_index(self.accumulator.argmax(), self.accumulator.shape)
        else:
          maxima = np.argwhere(self.accumulator == np.amax(self.accumulator))
          if len(maxima) == 1: i_max = maxima[0]
          elif len(maxima)==0 or (len(maxima) > 1 and self.HoughSpace_format == 'linearIntercepts'): return(-999, -999)
          elif len(maxima) > 1 and abs(min([x[1] for x in maxima]) - max([x[1] for x in maxima])) < res: i_max = maxima[0]
          else:
            maxima_slopesAxis_list = np.asarray(list([x[1] for x in maxima]))
            quantile = np.quantile(maxima_slopesAxis_list, self.n_quantile)
            Nwithin = 0
            for item in maxima:
               if abs(item[1]-quantile)< res: Nwithin += 1
            if Nwithin/len(maxima) > self.n_quantile:
               i_x = min([x[1] for x in maxima], key=lambda b: abs(b-quantile))
               for im in maxima:
                 if im[1] == i_x : i_max = im
            else: return(-999, -999)
        found_yH, found_xH = yH_bins[int(i_max[0])], xH_bins[int(i_max[1])]
        if self.HoughSpace_format == 'normal':
           slope = -1./np.tan(found_xH)
           interceptShift = found_yH/np.sin(found_xH)
           intercept = (np.tan(found_xH)*interceptShift + self.z_offset)/np.tan(found_xH)
        elif self.HoughSpace_format == 'linearSlopeIntercept':
           slope, intercept = found_xH, found_yH - found_xH*self.z_offset
        elif self.HoughSpace_format == 'linearIntercepts':
           slope = (found_xH - found_yH)/self.det_Zlen
           intercept = found_yH - slope*self.z_offset
        return (slope, intercept)

    def fit_randomize(self, hit_collection, hit_d, n_random, is_scaled, draw, weights = None) :
        if not len(hit_collection) : return (-1, -1)
        if (n_random > 0) :
            random_hit_collection = []
            for i_random in range(n_random) :
                random_hits = np.random.uniform(size = hit_collection.shape) - 0.5
                random_hits *= hit_d
                random_hits += hit_collection
                random_hit_collection.append(random_hits)
            random_hit_collection = np.concatenate(random_hit_collection)
            if weights is not None : weights = np.tile(weights, n_random)
            return self.fit(random_hit_collection, is_scaled, draw, weights)
        else : return self.fit(hit_collection, is_scaled, draw, weights)

def numPlanesHit(systems, detector_ids) :
    scifi_stations = detector_ids[systems == 0]//1000000
    mufi_ds_planes = (detector_ids[systems == 3]%10000)//1000
    mufi_us_planes = (detector_ids[systems == 2]%10000)//1000
    return len(np.unique(scifi_stations)) + len(np.unique(mufi_ds_planes)) + len(np.unique(mufi_us_planes))

class MuonReco(ROOT.FairTask) :
    " Muon reconstruction "
    def Init(self) :
        self.logger = ROOT.FairLogger.GetLogger()
        self.lsOfGlobals, self.ioman = ROOT.gROOT.GetListOfGlobals(), ROOT.FairRootManager.Instance()
        self.scifiDet, self.mufiDet = self.lsOfGlobals.FindObject('Scifi'), self.lsOfGlobals.FindObject('MuFilter')
        self.isMC = False if self.ioman.GetInTree().GetName() == 'rawConv' else True
        sink = self.ioman.GetSink()
        eventTree = sink.GetOutTree() if sink else None
        if eventTree: self.MuFilterHits, self.ScifiHits, self.EventHeader = eventTree.Digi_MuFilterHits, eventTree.Digi_ScifiHits, eventTree.EventHeader
        else: self.MuFilterHits, self.ScifiHits, self.EventHeader = self.ioman.GetInTree().Digi_MuFilterHits, self.ioman.GetInTree().Digi_ScifiHits, self.ioman.GetInTree().EventHeader
        self.scale, self.events_run = 1, 0
        tree = ET.parse(self.par_file); root = tree.getroot()
        if not hasattr(self, "genfitTrack"): self.genfitTrack = int(root[0].text)
        self.draw = int(root[1].text)
        for case in root.findall('tracking_case'):
            if case.get('name') == self.tracking_case:
               self.Scifi_meas, max_angle = int(case.find('use_Scifi_clust').text), float(case.find('max_angle').text)
               for rep in case.findall('Hough_space_format'):
                   if rep.get('name') == self.Hough_space_format:
                      n_yH, yH_min_xz, yH_max_xz = int(rep.find('N_yH_bins').text), float(rep.find('yH_min_xz').text), float(rep.find('yH_max_xz').text)
                      yH_min_yz, yH_max_yz = float(rep.find('yH_min_yz').text), float(rep.find('yH_max_yz').text)
                      n_xH, xH_min_xz, xH_max_xz = int(rep.find('N_xH_bins').text), float(rep.find('xH_min_xz').text), float(rep.find('xH_max_xz').text)
                      xH_min_yz, xH_max_yz = float(rep.find('xH_min_yz').text), float(rep.find('xH_max_yz').text)
               self.HT_space_scale = 1 if case.find('HT_space_scale')==None else float(case.find('HT_space_scale').text)
               self.n_random, self.muon_weight, self.min_planes_hit = int(case.find('n_random').text), int(case.find('mufi_weight').text), int(case.find('min_planes_hit').text)
               self.max_reco_muons, self.tolerance = int(case.find('max_reco_muons').text), float(case.find('tolerance').text)
               self.hits_to_fit, self.hits_for_triplet = case.find('hits_to_fit').text.strip(), (case.find('hits_for_hough').text.strip() if case.find('hits_to_validate')==None else case.find('hits_to_validate').text.strip())
               self.mask_plane, self.Nhits_per_plane = int(case.find('mask_plane').text), int(case.find('Nhits_per_plane').text)
               self.smooth_full, self.sigma, self.truncate = int(case.find('smooth_full').text), int(case.find('sigma').text), int(case.find('truncate').text)
               self.n_quantile, self.res = float(case.find('n_quantile').text), int(case.find('res').text)
        self.MuFilter_ds_dx, self.MuFilter_ds_dy, self.MuFilter_ds_dz = self.mufiDet.GetConfParF("MuFilter/DownstreamBarY"), self.mufiDet.GetConfParF("MuFilter/DownstreamBarY"), self.mufiDet.GetConfParF("MuFilter/DownstreamBarZ")
        self.MuFilter_us_dx, self.MuFilter_us_dy, self.MuFilter_us_dz = self.mufiDet.GetConfParF("MuFilter/UpstreamBarX"), self.mufiDet.GetConfParF("MuFilter/UpstreamBarY"), self.mufiDet.GetConfParF("MuFilter/UpstreamBarZ")
        self.Scifi_dx, self.Scifi_dy, self.Scifi_dz = self.scifiDet.GetConfParF("Scifi/channel_width"), self.scifiDet.GetConfParF("Scifi/channel_width"), self.scifiDet.GetConfParF("Scifi/epoxymat_z")
        self.MuFilter_us_nSiPMs = self.mufiDet.GetConfParI("MuFilter/UpstreamnSiPMs")*self.mufiDet.GetConfParI("MuFilter/UpstreamnSides")
        self.MuFilter_ds_nSiPMs_hor, self.MuFilter_ds_nSiPMs_vert = self.mufiDet.GetConfParI("MuFilter/DownstreamnSiPMs")*self.mufiDet.GetConfParI("MuFilter/DownstreamnSides"), self.mufiDet.GetConfParI("MuFilter/DownstreamnSiPMs")
        self.Scifi_nPlanes, self.DS_nPlanes = self.scifiDet.GetConfParI("Scifi/nscifi"), self.mufiDet.GetConfParI("MuFilter/NDownstreamPlanes")
        self.max_n_hits_plane, self.max_n_Scifi_hits = 3, 3*2*self.scifiDet.GetConfParI("Scifi/nscifi")
        if self.hits_for_triplet.find('sf') >= 0 and self.hits_for_triplet.find('ds') >= 0:
           det_Zlen = (self.mufiDet.GetConfParF("MuFilter/Muon9Dy") - self.scifiDet.GetConfParF("Scifi/Ypos0"))*unit.cm + 5.0*unit.cm
           z_offset = self.scifiDet.GetConfParF("Scifi/Ypos0")*unit.cm - 2.5*unit.cm
        elif self.hits_for_triplet == 'sf':
           det_Zlen = (self.scifiDet.GetConfParF("Scifi/Ypos4") - self.scifiDet.GetConfParF("Scifi/Ypos0"))*unit.cm + 5.0*unit.cm
           z_offset = self.scifiDet.GetConfParF("Scifi/Ypos0")*unit.cm - 2.5*unit.cm
        elif self.hits_for_triplet == 'ds':
           det_Zlen = (self.mufiDet.GetConfParF("MuFilter/Muon9Dy") - self.mufiDet.GetConfParF("MuFilter/Muon6Dy"))*unit.cm + 5.0*unit.cm
           z_offset = self.mufiDet.GetConfParF("MuFilter/Muon6Dy")*unit.cm - 2.5*unit.cm
        if self.tracking_case.find('nu_') >= 0: z_offset = 0*unit.cm
        if self.Hough_space_format == 'normal':
            self.h_ZX = hough(n_yH, [yH_min_xz, yH_max_xz], n_xH, [-max_angle+np.pi/2., max_angle+np.pi/2.], z_offset, self.Hough_space_format, self.HT_space_scale, det_Zlen)
            self.h_ZY = hough(n_yH, [yH_min_yz, yH_max_yz], n_xH, [-max_angle+np.pi/2., max_angle+np.pi/2.], z_offset, self.Hough_space_format, self.HT_space_scale, det_Zlen)
        else:
            self.h_ZX = hough(n_yH, [yH_min_xz, yH_max_xz], n_xH, [xH_min_xz, xH_max_xz], z_offset, self.Hough_space_format, self.HT_space_scale, det_Zlen)
            self.h_ZY = hough(n_yH, [yH_min_yz, yH_max_yz], n_xH, [xH_min_yz, xH_max_yz], z_offset, self.Hough_space_format, self.HT_space_scale, det_Zlen)
        self.h_ZX.smooth_full, self.h_ZY.smooth_full = self.smooth_full, self.smooth_full
        self.h_ZX.sigma, self.h_ZX.truncate, self.h_ZY.sigma, self.h_ZY.truncate = self.sigma, self.truncate, self.sigma, self.truncate
        self.h_ZX.n_quantile, self.h_ZX.res, self.h_ZY.n_quantile, self.h_ZY.res = self.n_quantile, self.res, self.n_quantile, self.res
        self.track_type = 11 if self.hits_to_fit == "sf" else (13 if self.hits_to_fit == "ds" else 15)
        self.a, self.b = ROOT.TVector3(), ROOT.TVector3()
        if self.ioman.GetObject('Reco_MuonTracks') != None: self.kalman_tracks = self.ioman.GetObject('Reco_MuonTracks')
        else:
           if self.genfitTrack:
              self.kalman_tracks = ROOT.TObjArray(10)
              if hasattr(self, "standalone") and self.standalone: self.ioman.Register("Reco_MuonTracks", self.ioman.GetFolderName(), self.kalman_tracks, ROOT.kTRUE)
           else:
              self.kalman_tracks = ROOT.TClonesArray("sndRecoTrack", 10)
              if hasattr(self, "standalone") and self.standalone: self.ioman.Register("Reco_MuonTracks", "", self.kalman_tracks, ROOT.kTRUE)
        if self.Scifi_meas: self.clusScifi = ROOT.TObjArray(100)
        geoMat, bfield = ROOT.genfit.TGeoMaterialInterface(), ROOT.genfit.ConstField(0, 0, 0)
        ROOT.genfit.FieldManager.getInstance().init(bfield)
        ROOT.genfit.MaterialEffects.getInstance().init(geoMat)
        ROOT.genfit.MaterialEffects.getInstance().setNoEffects()
        self.kalman_fitter = ROOT.genfit.KalmanFitter()
        self.kalman_fitter.setMaxIterations(50)
        self.kalman_sigmaScifi_spatial = self.Scifi_dx / 12**0.5
        self.kalman_sigmaMufiUS_spatial, self.kalman_sigmaMufiDS_spatial = self.MuFilter_us_dy / 12**0.5, self.MuFilter_ds_dy/ 12**0.5
        return 0

    def SetScaleFactor(self, s): self.scale = s
    def SetParFile(self, f): self.par_file = f
    def SetTrackingCase(self, c): self.tracking_case = c
    def SetHoughSpaceFormat(self, h): self.Hough_space_format = h
    def ForceGenfitTrackFormat(self): self.genfitTrack = 1
    def SetStandalone(self): self.standalone = 1

    def Exec(self, opt) :
        self.kalman_tracks.Clear('C'); self.events_run += 1
        hit_collection = {"pos" : [[], [], []], "d" : [[], [], []], "vert" : [], "index" : [], "system" : [], "detectorID" : [], "B" : [[], [], []], "time": [], "mask": []}
        if ("us" in self.hits_to_fit) or ("ds" in self.hits_to_fit) or ("ve" in self.hits_to_fit) :
            for i_h, muHit in enumerate(self.MuFilterHits) :
                if muHit.GetSystem() == 1 and "ve" not in self.hits_to_fit: continue
                if muHit.GetSystem() == 2 and "us" not in self.hits_to_fit: continue
                if muHit.GetSystem() == 3 and "ds" not in self.hits_to_fit: continue
                self.mufiDet.GetPosition(muHit.GetDetectorID(), self.a, self.b)
                hit_collection["pos"][0].append(self.a.X()); hit_collection["pos"][1].append(self.a.Y()); hit_collection["pos"][2].append(self.a.Z())
                hit_collection["B"][0].append(self.b.X()); hit_collection["B"][1].append(self.b.Y()); hit_collection["B"][2].append(self.b.Z())
                hit_collection["vert"].append(muHit.isVertical()); hit_collection["system"].append(muHit.GetSystem())
                hit_collection["d"][0].append(self.MuFilter_ds_dx); hit_collection["d"][2].append(self.MuFilter_ds_dz)
                hit_collection["index"].append(i_h); hit_collection["detectorID"].append(muHit.GetDetectorID()); hit_collection["mask"].append(False)
                times = []
                if muHit.GetSystem() == 3 :
                    hit_collection["d"][1].append(self.MuFilter_ds_dx)
                    for ch in range(self.MuFilter_ds_nSiPMs_hor):
                        if muHit.isVertical() and ch==self.MuFilter_ds_nSiPMs_vert: break
                        times.append(muHit.GetAllTimes()[ch] if self.isMC else muHit.GetAllTimes()[ch]*6.25)
                else :
                    hit_collection["d"][1].append(self.MuFilter_us_dy)
                    for ch in range(self.MuFilter_us_nSiPMs): times.append(muHit.GetAllTimes()[ch] if self.isMC else muHit.GetAllTimes()[ch]*6.25)
                hit_collection["time"].append(times)
        if "sf" in self.hits_to_fit :
            if self.Scifi_meas:
               self.clusScifi.Clear(); self.scifiCluster()
               for i_cl, cl in enumerate(self.clusScifi) :
                   cl.GetPosition(self.a,self.b)
                   hit_collection["pos"][0].append(self.a.X()); hit_collection["pos"][1].append(self.a.Y()); hit_collection["pos"][2].append(self.a.Z())
                   hit_collection["B"][0].append(self.b.X()); hit_collection["B"][1].append(self.b.Y()); hit_collection["B"][2].append(self.b.Z())
                   hit_collection["d"][0].append(cl.GetN()*self.Scifi_dx); hit_collection["d"][1].append(cl.GetN()*self.Scifi_dy); hit_collection["d"][2].append(self.Scifi_dz)
                   hit_collection["vert"].append(True if int(cl.GetFirst()/100000)%10==1 else False)
                   hit_collection["index"].append(i_cl); hit_collection["system"].append(0); hit_collection["detectorID"].append(cl.GetFirst()); hit_collection["mask"].append(False)
                   hit_collection["time"].append([cl.GetTime()/6.25 if self.isMC else cl.GetTime()])
            else:
                 for i_h, sfHit in enumerate(self.ScifiHits) :
                     if not sfHit.isValid(): continue
                     self.scifiDet.GetSiPMPosition(sfHit.GetDetectorID(), self.a, self.b)
                     hit_collection["pos"][0].append(self.a.X()); hit_collection["pos"][1].append(self.a.Y()); hit_collection["pos"][2].append(self.a.Z())
                     hit_collection["B"][0].append(self.b.X()); hit_collection["B"][1].append(self.b.Y()); hit_collection["B"][2].append(self.b.Z())
                     hit_collection["d"][0].append(self.Scifi_dx); hit_collection["d"][1].append(self.Scifi_dy); hit_collection["d"][2].append(self.Scifi_dz)
                     hit_collection["vert"].append(sfHit.isVertical()); hit_collection["index"].append(i_h); hit_collection["system"].append(0); hit_collection["detectorID"].append(sfHit.GetDetectorID()); hit_collection["mask"].append(False)
                     hit_collection["time"].append([sfHit.GetTime() if self.isMC else sfHit.GetTime()*6.25])
        if not hit_collection['pos'][0]: return
        for k, it in hit_collection.items() :
            if k in ['vert', 'mask']:
                dt = np.bool_
            elif k in ["index", "system", "detectorID"]:
                dt = np.int32
            elif k != 'time':
                dt = np.float32
            if k == 'time':
               ln = max(map(len, it))
               hit_collection[k] = np.stack(np.array([xi+[None]*(ln-len(xi)) for xi in it]), axis = 1)
            else: hit_collection[k] = np.array(it, dtype = dt)

        hit_collection["used"] = np.zeros(len(hit_collection["pos"][0]), dtype=np.bool_)
        triplet_sys = [0] if "sf" in self.hits_for_triplet else []
        for s, c in [("ve", 1), ("us", 2), ("ds", 3)]:
            if s in self.hits_for_triplet:
                triplet_sys.append(c)

        prev_ZY, prev_ZX = [], []

        for i_muon in range(self.max_reco_muons) :
            m_h = np.logical_and(~hit_collection["vert"], np.isin(hit_collection["system"], triplet_sys))
            m_v = np.logical_and(hit_collection["vert"], np.isin(hit_collection["system"], triplet_sys))

            def search_hough(is_v):
                m = np.logical_and(
                    np.logical_and(
                        hit_collection["vert"] == is_v, ~hit_collection["mask"]
                    ), ~hit_collection["used"]
                )
                ax = 0 if is_v else 1
                return (self.h_ZX if is_v else self.h_ZY).fit_randomize(
                    np.dstack(
                        [hit_collection["pos"][2][m], hit_collection["pos"][ax][m]]
                    )[0],
                    np.dstack(
                        [hit_collection["d"][2][m], hit_collection["d"][ax][m]]
                    )[0],
                    self.n_random,
                    False,
                    self.draw
                )

            h_ZY_new, h_ZX_new = search_hough(False), search_hough(True)
            tol = self.tolerance

            def eval_q(h, is_v, unused_only=True):
                if h[0] in [-1, -999]:
                    return 0
                m = m_v if is_v else m_h
                ax = 0 if is_v else 1

                # hits returned by hit_finder are relative to the input filtered arrays (hits in subset 'm')
                hits_rel = hit_finder(
                    h[0], h[1],
                    np.dstack([hit_collection["pos"][2][m], hit_collection["pos"][ax][m]]),
                    np.dstack([hit_collection["d"][2][m], hit_collection["d"][ax][m]]),
                    tol
                )
                if not len(hits_rel):
                    return 0

                # Filter hit properties for consistency
                sub_sys, sub_det, sub_used = hit_collection["system"][m], hit_collection["detectorID"][m], hit_collection["used"][m]
                if unused_only:
                    return numPlanesHit(sub_sys[hits_rel][~sub_used[hits_rel]], sub_det[hits_rel][~sub_used[hits_rel]])
                return numPlanesHit(sub_sys[hits_rel], sub_det[hits_rel])

            q_ZY, q_ZX = eval_q(h_ZY_new, False), eval_q(h_ZX_new, True)
            best_ZY, best_ZX = h_ZY_new, h_ZX_new

            # Primary Recycling: If one is strong but other is weak, prefer existing lines
            if q_ZY >= self.min_planes_hit and q_ZX < self.min_planes_hit:
                for p in prev_ZX:
                    if eval_q(p, True, False) >= self.min_planes_hit:
                        best_ZX, q_ZX = p, self.min_planes_hit
                        break
            elif q_ZX >= self.min_planes_hit and q_ZY < self.min_planes_hit:
                for p in prev_ZY:
                    if eval_q(p, False, False) >= self.min_planes_hit:
                        best_ZY, q_ZY = p, self.min_planes_hit
                        break

            if q_ZY < self.min_planes_hit or q_ZX < self.min_planes_hit:
                break

            track_hits_ZY = hit_finder(
                best_ZY[0], best_ZY[1],
                np.dstack([hit_collection["pos"][2][~hit_collection["vert"]], hit_collection["pos"][1][~hit_collection["vert"]]]),
                np.dstack([hit_collection["d"][2][~hit_collection["vert"]], hit_collection["d"][1][~hit_collection["vert"]]]),
                tol
            )
            track_hits_ZX = hit_finder(
                best_ZX[0], best_ZX[1],
                np.dstack([hit_collection["pos"][2][hit_collection["vert"]], hit_collection["pos"][0][hit_collection["vert"]]]),
                np.dstack([hit_collection["d"][2][hit_collection["vert"]], hit_collection["d"][0][hit_collection["vert"]]]),
                tol
            )

            posM, momM, covM = ROOT.TVector3(0, 0, 0.), ROOT.TVector3(0, 0, 100.), ROOT.TMatrixDSym(6)
            rv = self.kalman_sigmaScifi_spatial if self.hits_to_fit.find('sf') >= 0 else self.kalman_sigmaMufiDS_spatial

            for i in range(3):
                covM[i][i] = rv*rv
            for i in range(3,6):
                covM[i][i] = ROOT.TMath.Power(rv / 8. / ROOT.TMath.Sqrt(3), 2)

            rep = ROOT.genfit.RKTrackRep(13)

            state = ROOT.genfit.MeasuredStateOnPlane(rep)
            rep.setPosMomCov(state, posM, momM, covM)

            seedS, seedC = ROOT.TVectorD(6), ROOT.TMatrixDSym(6)
            rep.get6DStateCov(state, seedS, seedC)

            theTrack = ROOT.genfit.Track(rep, seedS, seedC)

            h_z = np.concatenate([hit_collection["pos"][2][hit_collection["vert"]][track_hits_ZX], hit_collection["pos"][2][~hit_collection["vert"]][track_hits_ZY]])
            h_A0 = np.concatenate([hit_collection["pos"][0][hit_collection["vert"]][track_hits_ZX], hit_collection["pos"][0][~hit_collection["vert"]][track_hits_ZY]])
            h_A1 = np.concatenate([hit_collection["pos"][1][hit_collection["vert"]][track_hits_ZX], hit_collection["pos"][1][~hit_collection["vert"]][track_hits_ZY]])
            h_B = [np.concatenate([hit_collection["B"][i][hit_collection["vert"]][track_hits_ZX], hit_collection["B"][i][~hit_collection["vert"]][track_hits_ZY]]) for i in range(3)]
            h_det = np.concatenate([hit_collection["detectorID"][hit_collection["vert"]][track_hits_ZX], hit_collection["detectorID"][~hit_collection["vert"]][track_hits_ZY]])
            k_sig = np.concatenate([hit_collection["d"][0][hit_collection["vert"]][track_hits_ZX] / 12**0.5, hit_collection["d"][1][~hit_collection["vert"]][track_hits_ZY] / 12**0.5])
            k_max = np.concatenate([((hit_collection["d"][0][hit_collection["vert"]][track_hits_ZX]/2.)**2 + (hit_collection["d"][2][hit_collection["vert"]][track_hits_ZX]/2.)**2)**0.5, ((hit_collection["d"][1][~hit_collection["vert"]][track_hits_ZY]/2.)**2 + (hit_collection["d"][2][~hit_collection["vert"]][track_hits_ZY]/2.)**2)**0.5])
            h_time = {ch: np.concatenate([hit_collection["time"][ch][hit_collection["vert"]][track_hits_ZX], hit_collection["time"][ch][~hit_collection["vert"]][track_hits_ZY]]) for ch in range(hit_collection["time"].shape[0])}

            for i in h_z.argsort() :
                tp = ROOT.genfit.TrackPoint()

                hitC = ROOT.TMatrixDSym(7)
                hitC[6][6] = k_sig[i]**2

                m = ROOT.genfit.WireMeasurement(
                    ROOT.TVectorD(
                        7, array('d', [h_A0[i], h_A1[i], h_z[i], h_B[0][i], h_B[1][i], h_B[2][i], 0.])
                    ),
                    hitC, 1, 6, tp
                )

                m.setMaxDistance(k_max[i])
                m.setDetId(int(h_det[i]))
                m.setHitId(int(i))
                tp.addRawMeasurement(m)
                theTrack.insertPoint(tp)

            self.kalman_fitter.processTrack(theTrack)
            if theTrack.getFitStatus().isFitConverged():
               theTrack.SetUniqueID(self.track_type)
               if self.genfitTrack:
                   self.kalman_tracks.Add(theTrack)
               else:
                  reco = ROOT.sndRecoTrack(theTrack)
                  pT = ROOT.std.vector(ROOT.std.vector('float'))()

                  for i_z in h_z.argsort():
                      pT.push_back(
                          [h_time[ch][i_z] for ch in range(len(h_time)) if h_time[ch][i_z] is not None]
                      )

                  reco.setRawMeasTimes(pT)
                  reco.setTrackType(self.track_type)

                  ROOT.std.swap(reco, self.kalman_tracks.ConstructedAt(i_muon))
                  theTrack.Delete()

               if best_ZY == h_ZY_new:
                   prev_ZY.append(best_ZY)
               if best_ZX == h_ZX_new:
                   prev_ZX.append(best_ZX)

            idx_ZX = np.where(np.in1d(hit_collection["detectorID"], hit_collection["detectorID"][hit_collection["vert"]][track_hits_ZX]))[0]
            idx_ZY = np.where(np.in1d(hit_collection["detectorID"], hit_collection["detectorID"][~hit_collection["vert"]][track_hits_ZY]))[0]
            hit_collection["used"][np.concatenate([idx_ZX, idx_ZY])] = True

    def FinishTask(self):
        print("Processed", self.events_run)
        if not self.genfitTrack:
            self.kalman_tracks.Delete()

    def scifiCluster(self):
       clusters, hitDict = [], {
           self.ScifiHits[k].GetDetectorID():
               k for k in range(self.ScifiHits.GetEntries()) if self.ScifiHits[k].isValid()
       }
       hitList = sorted(hitDict.keys())

       if hitList:
              tmp, cprev, last = [hitList[0]], hitList[0], len(hitList)-1
              hitvector = ROOT.std.vector("sndScifiHit*")()

              for i in range(len(hitList)):
                   if i==0 and len(hitList)>1:
                       continue

                   c, neighbour = hitList[i], (hitList[i]-cprev)==1

                   if not neighbour or c==hitList[last]:
                        hitvector.clear()
                        for aH in tmp:
                            hitvector.push_back(self.ScifiHits[hitDict[aH]])

                        clusters.append(ROOT.sndCluster(tmp[0], len(tmp), hitvector, self.scifiDet, False))

                        if c!=hitList[last]:
                            tmp = [c]
                        elif not neighbour:
                            hitvector.clear()
                            hitvector.push_back(self.ScifiHits[hitDict[c]])
                            clusters.append(ROOT.sndCluster(c, 1, hitvector, self.scifiDet, False))
                   cprev = c
       self.clusScifi.Delete()
       for c in clusters:
           self.clusScifi.Add(c)
