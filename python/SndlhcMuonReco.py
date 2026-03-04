import ROOT
import numpy as np
import scipy.ndimage
from array import array
import xml.etree.ElementTree as ET
import matplotlib.pyplot as plt
import shipunit as unit

def hit_finder(slope, intercept, box_centers, box_ds, tol = 0.) :
    """ Finds hits intersected by Hough line """

    # First check if track at center of box is within box limits
    d = np.abs(box_centers[0,:,1] - (box_centers[0,:,0]*slope + intercept))
    center_in_box = d < (box_ds[0,:,1]+tol)/2.

    # Now check if, assuming line is not in box at box center, the slope is large enough for line to clip the box at corner
    clips_corner = np.abs(slope) > np.abs((d - (box_ds[0,:,1]+tol)/2.)/(box_ds[0,:,0]+tol)/2.)
    
    # If either of these is true, line goes through hit:
    hit_mask = np.logical_or(center_in_box, clips_corner)

    # Return indices
    return np.where(hit_mask)[0]

class hough() :
    """ Hough transform implementation """

    def __init__(self, n_yH, yH_range, n_xH, xH_range, z_offset, Hformat, space_scale, det_Zlen, squaretheta = False, smooth = True) :

        self.n_yH = n_yH
        self.n_xH = n_xH

        self.yH_range = yH_range
        self.xH_range = xH_range

        self.z_offset = z_offset
        self.HoughSpace_format = Hformat
        self.space_scale = space_scale
        
        self.det_Zlen = det_Zlen

        self.smooth = smooth

        self.yH_bins = np.linspace(self.yH_range[0], self.yH_range[1], n_yH)
        if not squaretheta :
            self.xH_bins = np.linspace(self.xH_range[0], self.xH_range[1], n_xH)
        else :
            self.xH_bins = np.linspace(np.sign(self.xH_range[0])*(self.xH_range[0]**0.5), np.sign(self.xH_range[1])*(self.xH_range[1]**0.5), n_xH)
            self.xH_bins = np.sign(self.xH_bins)*np.square(self.xH_bins)

        self.cos_thetas = np.cos(self.xH_bins)
        self.sin_thetas = np.sin(self.xH_bins)
        
        self.xH_i = np.array(list(range(n_xH)))

        # A back-up Hough space designed to have more/less bins wrt the default one above.
        # It is useful when fitting some low-E muon tracks, which are curved due to mult. scattering.
        self.n_yH_scaled = int(n_yH*space_scale)
        self.n_xH_scaled = int(n_xH*space_scale)
        self.yH_bins_scaled = np.linspace(self.yH_range[0], self.yH_range[1], self.n_yH_scaled)
        if not squaretheta :
            self.xH_bins_scaled= np.linspace(self.xH_range[0], self.xH_range[1], self.n_xH_scaled)
        else :
            self.xH_bins_scaled = np.linspace(np.sign(self.xH_range[0])*(self.xH_range[0]**0.5), np.sign(self.xH_range[1])*(self.xH_range[1]**0.5), self.n_xH_scaled)
            self.xH_bins_scaled = np.sign(self.xH_bins_scaled)*np.square(self.xH_bins_scaled)

        self.cos_thetas_scaled = np.cos(self.xH_bins_scaled)
        self.sin_thetas_scaled = np.sin(self.xH_bins_scaled)

        self.xH_i_scaled = np.array(list(range(self.n_xH_scaled)))

    def fit(self, hit_collection, is_scaled, draw, weights = None) :

        if not is_scaled:
           n_xH = self.n_xH
           n_yH = self.n_yH
           cos_thetas = self.cos_thetas
           sin_thetas = self.sin_thetas
           xH_bins = self.xH_bins
           yH_bins = self.yH_bins
           xH_i = self.xH_i
           res = self.res
        else:
           n_xH = self.n_xH_scaled
           n_yH = self.n_yH_scaled
           cos_thetas = self.cos_thetas_scaled
           sin_thetas = self.sin_thetas_scaled
           xH_bins = self.xH_bins_scaled
           yH_bins = self.yH_bins_scaled
           xH_i = self.xH_i_scaled
           res = self.res*self.space_scale

        self.accumulator = np.zeros((n_yH, n_xH))
        for i_hit, hit in enumerate(hit_collection) :
            shifted_hitZ = hit[0] - self.z_offset
            if self.HoughSpace_format == 'normal':
                 hit_yH = shifted_hitZ*cos_thetas + hit[1]*sin_thetas
            elif self.HoughSpace_format == 'linearSlopeIntercept':
                 hit_yH = hit[1] - shifted_hitZ*xH_bins
            elif self.HoughSpace_format== 'linearIntercepts':
                 hit_yH = (self.det_Zlen*hit[1] - shifted_hitZ*xH_bins)/(self.det_Zlen - shifted_hitZ)
            out_of_range = np.logical_and(hit_yH > self.yH_range[0], hit_yH < self.yH_range[1]) 
            hit_yH_i = np.floor((hit_yH[out_of_range] - self.yH_range[0])/(self.yH_range[1] - self.yH_range[0])*n_yH).astype(np.int_)

            if weights is not None :
                self.accumulator[hit_yH_i, xH_i[out_of_range]] += weights[i_hit]
            else :
                self.accumulator[hit_yH_i, xH_i[out_of_range]] += 1

        # Smooth accumulator
        if self.smooth_full :
            self.accumulator = scipy.ndimage.gaussian_filter(self.accumulator, self.sigma, truncate=self.truncate)

        # This might be useful for debugging, but leave out for now.
        if draw :
            plt.figure()
            plt.imshow(self.accumulator, origin = "lower", extent = [self.xH_range[0], self.xH_range[-1], self.yH_range[0], self.yH_range[-1]], aspect = "auto")
            plt.tight_layout()
            plt.show()

        if self.smooth_full:
          i_max = np.unravel_index(self.accumulator.argmax(), self.accumulator.shape)
        else:
          maxima = np.argwhere(self.accumulator == np.amax(self.accumulator))
          if len(maxima) == 1:
            i_max = maxima[0]
          elif len(maxima)==0 or (len(maxima) > 1 and self.HoughSpace_format == 'linearIntercepts'):
               return(-999, -999)
          elif len(maxima) > 1 and abs(min([x[1] for x in maxima]) - max([x[1] for x in maxima])) < res:
               i_max = maxima[0]
          else:
            maxima_slopesAxis_list = np.asarray(list([x[1] for x in maxima]))
            quantile = np.quantile(maxima_slopesAxis_list, self.n_quantile)
            Nwithin = 0
            for item in maxima:
               if abs(item[1]-quantile)< res:
                 Nwithin += 1
            if Nwithin/len(maxima) > self.n_quantile: 
               i_x = min([x[1] for x in maxima], key=lambda b: abs(b-quantile))
               for im in maxima:
                 if im[1] == i_x : i_max = im
            else: return(-999, -999)

        found_yH = yH_bins[int(i_max[0])]
        found_xH = xH_bins[int(i_max[1])]
        
        if self.HoughSpace_format == 'normal':
           slope = -1./np.tan(found_xH)
           interceptShift = found_yH/np.sin(found_xH)
           intercept = (np.tan(found_xH)*interceptShift + self.z_offset)/np.tan(found_xH)
        elif self.HoughSpace_format == 'linearSlopeIntercept':
           slope = found_xH
           intercept = found_yH - slope*self.z_offset
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
    scifi_stations = []
    mufi_ds_planes = []
    mufi_us_planes = []
    scifi_stations.append( detector_ids[systems == 0]//1000000 )
    mufi_ds_planes.append( (detector_ids[systems == 3]%10000)//1000 )
    mufi_us_planes.append( (detector_ids[systems == 2]%10000)//1000 )
    return len(np.unique(scifi_stations)) + len(np.unique(mufi_ds_planes)) + len(np.unique(mufi_us_planes))
    
class MuonReco(ROOT.FairTask) :
    " Muon reconstruction "

    def Init(self) :
        self.logger = ROOT.FairLogger.GetLogger()
        self.lsOfGlobals  = ROOT.gROOT.GetListOfGlobals()
        self.scifiDet = self.lsOfGlobals.FindObject('Scifi')
        self.mufiDet = self.lsOfGlobals.FindObject('MuFilter')
        self.ioman = ROOT.FairRootManager.Instance()
        if self.ioman.GetInTree().GetName() == 'rawConv': self.isMC = False
        else: self.isMC = True
        sink = self.ioman.GetSink()
        eventTree = sink.GetOutTree() if sink else None
        if eventTree:
            self.MuFilterHits = eventTree.Digi_MuFilterHits
            self.ScifiHits       = eventTree.Digi_ScifiHits
            self.EventHeader        = eventTree.EventHeader
        else:
            self.MuFilterHits = self.ioman.GetInTree().Digi_MuFilterHits
            self.ScifiHits = self.ioman.GetInTree().Digi_ScifiHits
            self.EventHeader = self.ioman.GetInTree().EventHeader
        self.scale, self.events_run = 1, 0
        tree = ET.parse(self.par_file)
        root = tree.getroot()
        if not hasattr(self, "genfitTrack"): self.genfitTrack = int(root[0].text)
        self.draw = int(root[1].text)
        for case in root.findall('tracking_case'):
            if case.get('name') == self.tracking_case:
               self.Scifi_meas = int(case.find('use_Scifi_clust').text)
               max_angle = float(case.find('max_angle').text)
               for rep in case.findall('Hough_space_format'):
                   if rep.get('name') == self.Hough_space_format:
                      n_accumulator_yH = int(rep.find('N_yH_bins').text)
                      yH_min_xz, yH_max_xz = float(rep.find('yH_min_xz').text), float(rep.find('yH_max_xz').text)
                      yH_min_yz, yH_max_yz = float(rep.find('yH_min_yz').text), float(rep.find('yH_max_yz').text)
                      n_accumulator_xH = int(rep.find('N_xH_bins').text)
                      xH_min_xz, xH_max_xz = float(rep.find('xH_min_xz').text), float(rep.find('xH_max_xz').text)
                      xH_min_yz, xH_max_yz = float(rep.find('xH_min_yz').text), float(rep.find('xH_max_yz').text)
               self.HT_space_scale = 1 if case.find('HT_space_scale')==None else float(case.find('HT_space_scale').text)
               self.n_random = int(case.find('n_random').text)
               self.muon_weight = int(case.find('mufi_weight').text)
               self.min_planes_hit = int(case.find('min_planes_hit').text)
               self.max_reco_muons = int(case.find('max_reco_muons').text)
               self.tolerance = float(case.find('tolerance').text)
               self.hits_to_fit = case.find('hits_to_fit').text.strip()
               self.hits_for_triplet = case.find('hits_for_hough').text.strip() if case.find('hits_to_validate')==None else case.find('hits_to_validate').text.strip()
               self.mask_plane, self.Nhits_per_plane = int(case.find('mask_plane').text), int(case.find('Nhits_per_plane').text)
               self.smooth_full = int(case.find('smooth_full').text)
               self.sigma, self.truncate = int(case.find('sigma').text), int(case.find('truncate').text)
               self.n_quantile, self.res = float(case.find('n_quantile').text), int(case.find('res').text)
        self.MuFilter_ds_dx = self.mufiDet.GetConfParF("MuFilter/DownstreamBarY")
        self.MuFilter_ds_dy = self.mufiDet.GetConfParF("MuFilter/DownstreamBarY")
        self.MuFilter_ds_dz = self.mufiDet.GetConfParF("MuFilter/DownstreamBarZ")
        self.MuFilter_us_dx = self.mufiDet.GetConfParF("MuFilter/UpstreamBarX")
        self.MuFilter_us_dy = self.mufiDet.GetConfParF("MuFilter/UpstreamBarY")
        self.MuFilter_us_dz = self.mufiDet.GetConfParF("MuFilter/UpstreamBarZ")
        self.Scifi_dx = self.scifiDet.GetConfParF("Scifi/channel_width")
        self.Scifi_dy = self.scifiDet.GetConfParF("Scifi/channel_width")
        self.Scifi_dz = self.scifiDet.GetConfParF("Scifi/epoxymat_z")
        self.MuFilter_us_nSiPMs = self.mufiDet.GetConfParI("MuFilter/UpstreamnSiPMs")*self.mufiDet.GetConfParI("MuFilter/UpstreamnSides")
        self.MuFilter_ds_nSiPMs_hor = self.mufiDet.GetConfParI("MuFilter/DownstreamnSiPMs")*self.mufiDet.GetConfParI("MuFilter/DownstreamnSides")
        self.MuFilter_ds_nSiPMs_vert = self.mufiDet.GetConfParI("MuFilter/DownstreamnSiPMs")
        self.Scifi_nPlanes, self.DS_nPlanes = self.scifiDet.GetConfParI("Scifi/nscifi"), self.mufiDet.GetConfParI("MuFilter/NDownstreamPlanes")
        self.max_n_hits_plane = 3
        self.max_n_Scifi_hits = self.max_n_hits_plane*2*self.Scifi_nPlanes
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
            self.h_ZX = hough(n_accumulator_yH, [yH_min_xz, yH_max_xz], n_accumulator_xH, [-max_angle+np.pi/2., max_angle+np.pi/2.], z_offset, self.Hough_space_format, self.HT_space_scale, det_Zlen)
            self.h_ZY = hough(n_accumulator_yH, [yH_min_yz, yH_max_yz], n_accumulator_xH, [-max_angle+np.pi/2., max_angle+np.pi/2.], z_offset, self.Hough_space_format, self.HT_space_scale, det_Zlen)
        else:
            self.h_ZX = hough(n_accumulator_yH, [yH_min_xz, yH_max_xz], n_accumulator_xH, [xH_min_xz, xH_max_xz], z_offset, self.Hough_space_format, self.HT_space_scale, det_Zlen)
            self.h_ZY = hough(n_accumulator_yH, [yH_min_yz, yH_max_yz], n_accumulator_xH, [xH_min_yz, xH_max_yz], z_offset, self.Hough_space_format, self.HT_space_scale, det_Zlen)
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

    def SetScaleFactor(self, scale):
        self.scale = scale

    def SetParFile(self, file_name):
        self.par_file = file_name
    
    def SetTrackingCase(self, case):
        self.tracking_case = case

    def SetHoughSpaceFormat(self, Hspace_format):
        self.Hough_space_format = Hspace_format

    def ForceGenfitTrackFormat(self):
        self.genfitTrack = 1

    def SetStandalone(self):
        self.standalone = 1

    def Exec(self, opt) :
        self.kalman_tracks.Clear('C')
        self.events_run += 1
        hit_collection = {"pos" : [[], [], []], "d" : [[], [], []], "vert" : [], "index" : [], "system" : [], "detectorID" : [], "B" : [[], [], []], "time": [], "mask": []}
        if ("us" in self.hits_to_fit) or ("ds" in self.hits_to_fit) or ("ve" in self.hits_to_fit) :
            for i_hit, muFilterHit in enumerate(self.MuFilterHits) :
                if muFilterHit.GetSystem() == 1 and "ve" not in self.hits_to_fit : continue
                if muFilterHit.GetSystem() == 2 and "us" not in self.hits_to_fit : continue
                if muFilterHit.GetSystem() == 3 and "ds" not in self.hits_to_fit : continue
                self.mufiDet.GetPosition(muFilterHit.GetDetectorID(), self.a, self.b)
                hit_collection["pos"][0].append(self.a.X()); hit_collection["pos"][1].append(self.a.Y()); hit_collection["pos"][2].append(self.a.Z())
                hit_collection["B"][0].append(self.b.X()); hit_collection["B"][1].append(self.b.Y()); hit_collection["B"][2].append(self.b.Z())
                hit_collection["vert"].append(muFilterHit.isVertical()); hit_collection["system"].append(muFilterHit.GetSystem())
                hit_collection["d"][0].append(self.MuFilter_ds_dx); hit_collection["d"][2].append(self.MuFilter_ds_dz)
                hit_collection["index"].append(i_hit); hit_collection["detectorID"].append(muFilterHit.GetDetectorID()); hit_collection["mask"].append(False)
                times = []
                if muFilterHit.GetSystem() == 3 :
                    hit_collection["d"][1].append(self.MuFilter_ds_dx)
                    for ch in range(self.MuFilter_ds_nSiPMs_hor):
                        if muFilterHit.isVertical() and ch==self.MuFilter_ds_nSiPMs_vert: break
                        times.append(muFilterHit.GetAllTimes()[ch] if self.isMC else muFilterHit.GetAllTimes()[ch]*6.25)
                else :
                    hit_collection["d"][1].append(self.MuFilter_us_dy)
                    for ch in range(self.MuFilter_us_nSiPMs): times.append(muFilterHit.GetAllTimes()[ch] if self.isMC else muFilterHit.GetAllTimes()[ch]*6.25)
                hit_collection["time"].append(times)
        if "sf" in self.hits_to_fit :
            if self.Scifi_meas:
               self.clusScifi.Clear(); self.scifiCluster()
               for i_clust, scifiCl in enumerate(self.clusScifi) :
                   scifiCl.GetPosition(self.a,self.b)
                   hit_collection["pos"][0].append(self.a.X()); hit_collection["pos"][1].append(self.a.Y()); hit_collection["pos"][2].append(self.a.Z())
                   hit_collection["B"][0].append(self.b.X()); hit_collection["B"][1].append(self.b.Y()); hit_collection["B"][2].append(self.b.Z())
                   hit_collection["d"][0].append(scifiCl.GetN()*self.Scifi_dx); hit_collection["d"][1].append(scifiCl.GetN()*self.Scifi_dy); hit_collection["d"][2].append(scifiCl.GetN()*self.Scifi_dz)
                   hit_collection["vert"].append(True if int(scifiCl.GetFirst()/100000)%10==1 else False)
                   hit_collection["index"].append(i_clust); hit_collection["system"].append(0); hit_collection["detectorID"].append(scifiCl.GetFirst()); hit_collection["mask"].append(False)
                   hit_collection["time"].append([scifiCl.GetTime()/6.25 if self.isMC else scifiCl.GetTime()])
            else:
                 for i_hit, scifiHit in enumerate(self.ScifiHits) :
                     if not scifiHit.isValid(): continue 
                     self.scifiDet.GetSiPMPosition(scifiHit.GetDetectorID(), self.a, self.b)
                     hit_collection["pos"][0].append(self.a.X()); hit_collection["pos"][1].append(self.a.Y()); hit_collection["pos"][2].append(self.a.Z())
                     hit_collection["B"][0].append(self.b.X()); hit_collection["B"][1].append(self.b.Y()); hit_collection["B"][2].append(self.b.Z())
                     hit_collection["d"][0].append(self.Scifi_dx); hit_collection["d"][1].append(self.Scifi_dy); hit_collection["d"][2].append(self.Scifi_dz)
                     hit_collection["vert"].append(scifiHit.isVertical()); hit_collection["index"].append(i_hit); hit_collection["system"].append(0); hit_collection["detectorID"].append(scifiHit.GetDetectorID()); hit_collection["mask"].append(False)
                     hit_collection["time"].append([scifiHit.GetTime() if self.isMC else scifiHit.GetTime()*6.25])
        if len(hit_collection['pos'][0])==0: return
        for key, item in hit_collection.items() :
            if key == 'vert' or key == 'mask': this_dtype = np.bool_
            elif key in ["index", "system", "detectorID"]: this_dtype = np.int32
            elif key != 'time': this_dtype = np.float32
            if key == 'time':
               length = max(map(len, item))
               hit_collection[key] = np.stack(np.array([xi+[None]*(length-len(xi)) for xi in item]), axis = 1)
            else: hit_collection[key] = np.array(item, dtype = this_dtype)
        triplet_condition_system = [0] if "sf" in self.hits_for_triplet else []
        for sys, code in [("ve", 1), ("us", 2), ("ds", 3)]:
            if sys in self.hits_for_triplet: triplet_condition_system.append(code)
        
        for i_muon in range(self.max_reco_muons) :
            triplet_h = np.logical_and(~hit_collection["vert"], np.isin(hit_collection["system"], triplet_condition_system))
            triplet_v = np.logical_and(hit_collection["vert"], np.isin(hit_collection["system"], triplet_condition_system))
            if numPlanesHit(hit_collection["system"][triplet_h], hit_collection["detectorID"][triplet_h]) < self.min_planes_hit and numPlanesHit(hit_collection["system"][triplet_v], hit_collection["detectorID"][triplet_v]) < self.min_planes_hit: break
            
            mu_h = np.logical_and(np.logical_and(~hit_collection["vert"], ~hit_collection["mask"]), np.isin(hit_collection["system"], [1, 2, 3]))
            mu_v = np.logical_and(np.logical_and(hit_collection["vert"], ~hit_collection["mask"]), np.isin(hit_collection["system"], [1, 2, 3]))
            sf_h = np.logical_and(np.logical_and(~hit_collection["vert"], ~hit_collection["mask"]), np.isin(hit_collection["system"], [0]))
            sf_v = np.logical_and(np.logical_and(hit_collection["vert"], ~hit_collection["mask"]), np.isin(hit_collection["system"], [0]))

            def prepare_hough(m, s, ax):
                p = np.dstack([np.concatenate([np.tile(hit_collection["pos"][2][m], self.muon_weight), hit_collection["pos"][2][s]]), np.concatenate([np.tile(hit_collection["pos"][ax][m], self.muon_weight), hit_collection["pos"][ax][s]])])[0]
                d = np.dstack([np.concatenate([np.tile(hit_collection["d"][2][m], self.muon_weight), hit_collection["d"][2][s]]), np.concatenate([np.tile(hit_collection["d"][ax][m], self.muon_weight), hit_collection["d"][ax][s]])])[0]
                return p, d

            ZY_p, ZY_d = prepare_hough(mu_h, sf_h, 1)
            ZX_p, ZX_d = prepare_hough(mu_v, sf_v, 0)
            ZY_hough = self.h_ZY.fit_randomize(ZY_p, ZY_d, self.n_random, False, self.draw)
            ZX_hough = self.h_ZX.fit_randomize(ZX_p, ZX_d, self.n_random, False, self.draw)
            tol = self.tolerance

            track_hits_tr_ZY = hit_finder(ZY_hough[0], ZY_hough[1], np.dstack([hit_collection["pos"][2][triplet_h], hit_collection["pos"][1][triplet_h]]), np.dstack([hit_collection["d"][2][triplet_h], hit_collection["d"][1][triplet_h]]), tol)
            track_hits_tr_ZX = hit_finder(ZX_hough[0], ZX_hough[1], np.dstack([hit_collection["pos"][2][triplet_v], hit_collection["pos"][0][triplet_v]]), np.dstack([hit_collection["d"][2][triplet_v], hit_collection["d"][0][triplet_v]]), tol)
            n_zy, n_zx = numPlanesHit(hit_collection["system"][triplet_h][track_hits_tr_ZY], hit_collection["detectorID"][triplet_h][track_hits_tr_ZY]), numPlanesHit(hit_collection["system"][triplet_v][track_hits_tr_ZX], hit_collection["detectorID"][triplet_v][track_hits_tr_ZX])
            
            if (n_zy < self.min_planes_hit and n_zx < self.min_planes_hit) or (n_zy < 2 or n_zx < 2): break

            track_hits_ZY = hit_finder(ZY_hough[0], ZY_hough[1], np.dstack([hit_collection["pos"][2][~hit_collection["vert"]], hit_collection["pos"][1][~hit_collection["vert"]]]), np.dstack([hit_collection["d"][2][~hit_collection["vert"]], hit_collection["d"][1][~hit_collection["vert"]]]), tol)
            track_hits_ZX = hit_finder(ZX_hough[0], ZX_hough[1], np.dstack([hit_collection["pos"][2][hit_collection["vert"]], hit_collection["pos"][0][hit_collection["vert"]]]), np.dstack([hit_collection["d"][2][hit_collection["vert"]], hit_collection["d"][0][hit_collection["vert"]]]), tol)
            
            posM, momM, covM = ROOT.TVector3(0, 0, 0.), ROOT.TVector3(0,0,100.), ROOT.TMatrixDSym(6)
            res_val = self.kalman_sigmaScifi_spatial if self.hits_to_fit.find('sf') >= 0 else self.kalman_sigmaMufiDS_spatial
            for i in range(3): covM[i][i] = res_val*res_val
            for i in range(3,6): covM[i][i] = ROOT.TMath.Power(res_val / 8. / ROOT.TMath.Sqrt(3), 2)
            rep = ROOT.genfit.RKTrackRep(13)
            state = ROOT.genfit.MeasuredStateOnPlane(rep)
            rep.setPosMomCov(state, posM, momM, covM)
            seedState, seedCov = ROOT.TVectorD(6), ROOT.TMatrixDSym(6)
            rep.get6DStateCov(state, seedState, seedCov)
            theTrack = ROOT.genfit.Track(rep, seedState, seedCov)

            hit_z = np.concatenate([hit_collection["pos"][2][hit_collection["vert"]][track_hits_ZX], hit_collection["pos"][2][~hit_collection["vert"]][track_hits_ZY]])
            hit_A0 = np.concatenate([hit_collection["pos"][0][hit_collection["vert"]][track_hits_ZX], hit_collection["pos"][0][~hit_collection["vert"]][track_hits_ZY]])
            hit_A1 = np.concatenate([hit_collection["pos"][1][hit_collection["vert"]][track_hits_ZX], hit_collection["pos"][1][~hit_collection["vert"]][track_hits_ZY]])
            hit_B0, hit_B1, hit_B2 = np.concatenate([hit_collection["B"][0][hit_collection["vert"]][track_hits_ZX], hit_collection["B"][0][~hit_collection["vert"]][track_hits_ZY]]), np.concatenate([hit_collection["B"][1][hit_collection["vert"]][track_hits_ZX], hit_collection["B"][1][~hit_collection["vert"]][track_hits_ZY]]), np.concatenate([hit_collection["B"][2][hit_collection["vert"]][track_hits_ZX], hit_collection["B"][2][~hit_collection["vert"]][track_hits_ZY]])
            hit_detid = np.concatenate([hit_collection["detectorID"][hit_collection["vert"]][track_hits_ZX], hit_collection["detectorID"][~hit_collection["vert"]][track_hits_ZY]])
            kalman_sigma = np.concatenate([hit_collection["d"][0][hit_collection["vert"]][track_hits_ZX] / 12**0.5, hit_collection["d"][1][~hit_collection["vert"]][track_hits_ZY] / 12**0.5])
            kalman_max_d = np.concatenate([((hit_collection["d"][0][hit_collection["vert"]][track_hits_ZX]/2.)**2 + (hit_collection["d"][2][hit_collection["vert"]][track_hits_ZX]/2.)**2)**0.5, ((hit_collection["d"][1][~hit_collection["vert"]][track_hits_ZY]/2.)**2 + (hit_collection["d"][2][~hit_collection["vert"]][track_hits_ZY]/2.)**2)**0.5])
            hitID, hit_time = 0, {}
            for ch in range(hit_collection["time"].shape[0]): hit_time[ch] = np.concatenate([hit_collection["time"][ch][hit_collection["vert"]][track_hits_ZX], hit_collection["time"][ch][~hit_collection["vert"]][track_hits_ZY]])

            for i in hit_z.argsort() :
                tp = ROOT.genfit.TrackPoint(); hitCov = ROOT.TMatrixDSym(7); hitCov[6][6] = kalman_sigma[i]**2
                meas = ROOT.genfit.WireMeasurement(ROOT.TVectorD(7, array('d', [hit_A0[i], hit_A1[i], hit_z[i], hit_B0[i], hit_B1[i], hit_B2[i], 0.])), hitCov, 1, 6, tp)
                meas.setMaxDistance(kalman_max_d[i]); meas.setDetId(int(hit_detid[i])); meas.setHitId(int(hitID)); hitID += 1; tp.addRawMeasurement(meas); theTrack.insertPoint(tp)

            self.kalman_fitter.processTrack(theTrack)
            fitStatus = theTrack.getFitStatus()
            if fitStatus.isFitConverged():
               theTrack.SetUniqueID(self.track_type)
               if self.genfitTrack: self.kalman_tracks.Add(theTrack)
               else :
                  this_track = ROOT.sndRecoTrack(theTrack)
                  pointTimes = ROOT.std.vector(ROOT.std.vector('float'))()
                  for i_z in hit_z.argsort():
                      t_per_hit = [hit_time[ch][i_z] for ch in range(len(hit_time)) if hit_time[ch][i_z] != None]
                      pointTimes.push_back(t_per_hit)
                  this_track.setRawMeasTimes(pointTimes); this_track.setTrackType(self.track_type)
                  ROOT.std.swap(this_track, self.kalman_tracks.ConstructedAt(i_muon)); theTrack.Delete()
            
            idx_ZX = np.where(np.in1d(hit_collection["detectorID"], hit_collection["detectorID"][hit_collection["vert"]][track_hits_ZX]))[0]
            idx_ZY = np.where(np.in1d(hit_collection["detectorID"], hit_collection["detectorID"][~hit_collection["vert"]][track_hits_ZY]))[0]
            idx = np.concatenate([idx_ZX, idx_ZY])
            for k in hit_collection.keys():
                if len(hit_collection[k].shape) == 1: hit_collection[k] = np.delete(hit_collection[k], idx)
                elif len(hit_collection[k].shape) == 2: hit_collection[k] = np.delete(hit_collection[k], idx, axis = 1)

    def FinishTask(self) :
        print("Processed" ,self.events_run)
        if not self.genfitTrack : self.kalman_tracks.Delete()

    def scifiCluster(self):
       clusters, hitDict = [], {}
       for k in range(self.ScifiHits.GetEntries()):
            if self.ScifiHits[k].isValid(): hitDict[self.ScifiHits[k].GetDetectorID()] = k
       hitList = sorted(hitDict.keys())
       if len(hitList)>0:
              tmp, cprev, last = [hitList[0]], hitList[0], len(hitList)-1
              hitvector = ROOT.std.vector("sndScifiHit*")()
              for i in range(len(hitList)):
                   if i==0 and len(hitList)>1: continue
                   c, neighbour = hitList[i], (hitList[i]-cprev)==1
                   if not neighbour or c==hitList[last]:
                        hitvector.clear()
                        for aHit in tmp: hitvector.push_back(self.ScifiHits[hitDict[aHit]])
                        clusters.append(ROOT.sndCluster(tmp[0], len(tmp), hitvector, self.scifiDet, False))
                        if c!=hitList[last]: tmp = [c]
                        elif not neighbour:
                            hitvector.clear(); hitvector.push_back(self.ScifiHits[hitDict[c]])
                            clusters.append(ROOT.sndCluster(c, 1, hitvector, self.scifiDet, False))
                   cprev = c
       self.clusScifi.Delete()
       for c in clusters: self.clusScifi.Add(c)
