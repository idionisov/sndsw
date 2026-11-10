import ROOT, os, csv, json
from AnalysisFunctions import Analysis
from pathlib import Path
import numpy as np
from scipy.optimize import curve_fit, minimize_scalar

def tw_func(qdc, t0, alpha, beta, qdc0, gamma):
    denom = beta * qdc - qdc0
    denom = np.where(denom <= 1e-4, 1e-4, denom)
    return t0 + alpha / denom + gamma * qdc

def nll_objective(sigma_sys, x_data, y_data, sigma_stat, params):
    y_pred = tw_func(x_data, *params)
    var_total = sigma_stat**2 + sigma_sys**2
    return 0.5 * np.sum(((y_data - y_pred) ** 2) / var_total + np.log(var_total))

deciles=[i/10 for i in range(11)]
SiPM_distance_cut=5

ROOT.gStyle.SetPalette(ROOT.kRainBow)
ROOT.gStyle.SetNumberContours(255)
ROOT.gStyle.SetTitleXOffset(0.7)
ROOT.gStyle.SetTitleYOffset(0.6)
for i in ('x', 'y', 'z','t'): ROOT.gStyle.SetTitleSize(0.04, i)

class TimeWalk(ROOT.FairTask):

    def __init__(self, options, monitor):

        self.M=monitor
        self.simulation = options.simulation
        self.muAna = Analysis(options)
        self.muAna.SetTask(self)
        if options.numuStudy or options.nueStudy: self.muAna.nuStudy=True
        else: self.muAna.nuStudy=False
        self.options=options

        lsOfGlobals=ROOT.gROOT.GetListOfGlobals()
        self.MuFilter=lsOfGlobals.FindObject('MuFilter')
        self.Scifi=lsOfGlobals.FindObject('Scifi')
        self.nav=ROOT.gGeoManager.GetCurrentNavigator()
        self.runNr = str(options.runNumber).zfill(6)
        self.outputDir = getattr(options, 'outputDir', './Timing').rstrip('/')
        self.outpath = self.outputDir + '/'
        self.paramDir = getattr(options, 'paramDir', None)
        self.EventNumber=-1
        self.subsystemdict={1:'Veto', 2:'US', 3:'DS'}
        self.nchs={1:224, 2:800}

        self.systemAndPlanes = {1:2,2:5,3:7}
        self.systemAndBars={1:7,2:10,3:60}
        self.systemAndChannels={1:[8,0],2:[6,2],3:[1,0]}
        self.sdict={0:'Scifi',1:'Veto',2:'US',3:'DS'}
        self.zPos=self.M.zPos
        self.muAna.zPos=self.zPos
        self.A, self.B = ROOT.TVector3(), ROOT.TVector3()

        freq=160.316E6
        self.TDC2ns=1E9/freq

        self.largeSiPMmap={0:0 ,1:1 ,3:2 ,4:3 ,6:4 ,7:5}
        self.verticalBarDict={0:1, 1:3, 2:5, 3:6}

        self.hists=self.M.h

        run=ROOT.FairRunAna.Instance()

        ioman=ROOT.FairRootManager.Instance()
        self.OT=ioman.GetSink().GetOutTree()

        ref_arg = str(getattr(options, 'referencesystem', '0')).lower()
        if ref_arg in ('1', 'sf', 'scifi', 'sf_only'):
            self.active_refsystems = ['SF']
            self.referencesystem = 1
            self.refsysname = 'SF'
        elif ref_arg in ('3', '2', 'ds', 'ds_only'):
            self.active_refsystems = ['DS']
            self.referencesystem = 3
            self.refsysname = 'DS'
        else:  # 0, 'both', default
            self.active_refsystems = ['DS', 'SF']
            self.referencesystem = 3
            self.refsysname = 'DS'
        self.mode = 'calibration'
        self.state = 'uncorrected'
        self.step = 1
        self.path = 'TI18'
        if options.path.find('testbeam_June2023_H8') > 0:
            self.path = 'H8'

        #### Time walk correction run is independent of time alignment settings.
        self.TWCorrectionRun = str(5408).zfill(6)

        if self.simulation:
            self.timealignment = 'sim'
            self.GetLumis() # loads the simulation luminosities dictionary into tw.simlumis
            self.simMode = options.simMode
            self.cutdists = self.muAna.GetCutDistributions("005408", ('dy', ))

        elif not self.simulation:
            self.timealignment = self.muAna.GetTimeAlignmentType(runNr=self.runNr)

            if options.debug: self.trackevents = []

            ### If no time-walk correction run is provided. Set the default correction run depending on time alignment of the data set
            if options.AlignmentRun is None or self.muAna.GetTimeAlignmentType(runNr=options.AlignmentRun) != self.timealignment:
                if self.timealignment == 'old': self.AlignmentRun = str(5097).zfill(6)
                elif self.timealignment == 'new': self.AlignmentRun = str(5408).zfill(6)
                elif self.timealignment == 'new+LHCsynch': self.AlignmentRun = str(5999).zfill(6)
            self.muAna.AlignmentRun = self.AlignmentRun
            self.cutdists = {}

        self.tds0_min = getattr(options, 'tds0Min', None)
        self.tds0_max = getattr(options, 'tds0Max', None)
        if self.tds0_min is None or self.tds0_max is None:
            time_align = getattr(self, 'timealignment', 'new')
            # Category I ('old'): [20.25, 24.15] ns
            # Category II ('new'): [13.0, 15.5] ns
            # Category III ('new+LHCsynch'): [18.5, 20.5] ns
            if time_align == 'old':
                default_min, default_max = (20.25, 24.15)
            elif time_align == 'new':
                default_min, default_max = (13.0, 15.5)
            else:
                default_min, default_max = (18.5, 20.5)
            if self.tds0_min is None: self.tds0_min = default_min
            if self.tds0_max is None: self.tds0_max = default_max
        self.tds0_mode = None

        if self.mode == 'systemalignment':
            from systemalignment import SystemAlignment
            self.sa = SystemAlignment(options, self)

        if self.mode in ['reconstructmuonposition', 'multimuonreconstruction']:
            from systemalignment import SystemAlignment
            self.sa = SystemAlignment(options, self)

        elif self.mode == 'showerprofiles':
            self.barycentredata={}
            from showerprofiles import ShowerProfiles
            self.sp = ShowerProfiles(options, self)
            if options.fname:
                if options.fname.find('stage2cuts')>-1: options.numuStudy=True

        elif self.mode == 'selectioncriteria':
            from selectioncriteria import MuonSelectionCriteria as SelectionCriteria
            self.sc = SelectionCriteria(options, self)
            # self.sc.cuteffect['no cuts'] = self.M.eventTree.GetEntries()

        elif self.mode.find('extendedreconstruction')>-1:
            from extendedmuonreconstruction import ExtendedMuonReconstruction as ExtendedMuonReconstruction
            self.emr = ExtendedMuonReconstruction(options, self)
        elif self.mode=='struckquark':
            from extendedmuonreconstruction import QuarkVectorExtrapolation as QuarkVectorExtrapolation
            self.qve = QuarkVectorExtrapolation(options, self)

        # Histogram formatting configuration with built-in robust defaults
        histformat_path = getattr(options, 'histFormatFile', None)
        if not histformat_path:
            candidate = os.path.join(self.outputDir, 'TWhistogramformatting.json')
            if os.path.exists(candidate):
                histformat_path = candidate

        if histformat_path and os.path.exists(histformat_path):
            with open(histformat_path, 'r') as x:
                self.histformatting = json.load(x)
        else:
            self.histformatting = {
                "dtvxpred": {
                    "uncorrected": [(100, -80, 20), (200, 0, 50)],
                    "corrected": [(100, -80, 20), (200, -25, 25)],
                    "DS": {
                        "uncorrected": [(100, -80, 20), (200, 0, 50)],
                        "corrected": [(100, -80, 20), (200, -25, 25)]
                    },
                    "SF": {
                        "uncorrected": [(100, -80, 20), (200, -50, 0)],
                        "corrected": [(100, -80, 20), (200, -50, 0)]
                    }
                },
                "dtvqdc": {
                    "uncorrected": [(200, 0, 200), (350, -10, 25)],
                    "corrected": [(200, 0, 200), (350, -10, 25)],
                    "DS": {
                        "uncorrected": [(200, 0, 200), (350, -10, 25)],
                        "corrected": [(200, 0, 200), (350, -10, 25)]
                    },
                    "SF": {
                        "uncorrected": [(200, 0, 200), (500, -45, 5)],
                        "corrected": [(200, 0, 200), (500, -45, 5)]
                    }
                }
            }

        self.notInDS=0

        self.trackRequired=True if self.mode in('zeroth', 'tof', 'full', 'tw', 'res', 'systemalignment', 'reconstructmuonposition', 'selectioncriteria', 'showerprofiles', 'tds0-studies') else False

    def get_histformat(self, hist_type, state, ref_name=None):
        ref = ref_name if ref_name in ('DS', 'SF') else getattr(self, 'refsysname', 'DS')
        if hist_type in self.histformatting:
            entry = self.histformatting[hist_type]
            if isinstance(entry, dict):
                if ref in entry and isinstance(entry[ref], dict) and state in entry[ref]:
                    return entry[ref][state]
                if state in entry and isinstance(entry[state], (list, tuple)):
                    return entry[state]
        if hist_type == "dtvxpred":
            if ref == "SF":
                return [(100, -80, 20), (200, -50, 0)]
            return [(100, -80, 20), (200, 0, 50)]
        elif hist_type == "dtvqdc":
            if ref == "SF":
                return [(200, 0, 200), (500, -45, 5)]
            return [(200, 0, 200), (350, -10, 25)]
        return [(100, -80, 20), (200, 0, 50)]

    def GetEntries(self):
        return self.eventTree.GetEntries()

    def GetEventScaledWeight(self):
        w = self.M.eventTree.GetWeight()
        luminosity = self.simlumis[self.simMode]
        return np.float32(1/(w*luminosity))

    def GetLumis(self):
        simlumi_filename = getattr(self.options, 'simLumiFile', None)
        if not simlumi_filename or not os.path.exists(simlumi_filename):
            candidate = os.path.join(self.outputDir, 'simulationluminosities.json')
            if os.path.exists(candidate):
                simlumi_filename = candidate
        if simlumi_filename and os.path.exists(simlumi_filename):
            with open(simlumi_filename, 'r') as x:
                self.simlumis = json.load(x)
        else:
            self.simlumis = {}

    def ExecuteEvent(self, event):

        # Only accept events associated with colliding IP1 bunches for LHC synchronous runs
        if not self.simulation:
            # M.xing does not exist for simulation runs
            if self.muAna.isLHCsynch(self.runNr):
                if not self.M.xing['IP1']: return

        if not hasattr(self.muAna, 'task'): self.muAna.SetTask(self)

        tracks={1:[], 3:[]}
        inScifi, inDS = False, False

        hits = event.Digi_MuFilterHits
        scifi_hits = event.Digi_ScifiHits

        if self.simulation:
            self.scaled_event_weight = self.GetEventScaledWeight()
            histname = f'{self.simMode}_eventweight'
            if not histname in self.hists:
                title = f'Event weight for: {self.simMode};Event weight;Counts'
                self.hists[histname] = ROOT.TH1D(histname, title, 200, -1E3, 1E3)
            self.hists[histname].Fill(self.M.eventTree.GetWeight())

        for i,track in enumerate(self.M.Reco_MuonTracks):
            if any([not track.getFitStatus().isFitConverged(), track.getFitStatus().getNdf()==0]):
                continue
            track_type = track.getTrackType() if hasattr(track, 'getTrackType') else track.GetUniqueID()
            if track_type == 1:
                inScifi = True
                tracks[1].append(self.M.Reco_MuonTracks[i])
            if track_type == 3:
                inDS = True
                tracks[3].append(self.M.Reco_MuonTracks[i])

        if self.mode in ('zeroth', 'full', 'calibration', 'tof'):
            contexts = {}
            if 'DS' in self.active_refsystems and tracks[3] and not self.simulation:
                ds_track = tracks[3][0] if len(tracks[3]) == 1 else min(tracks[3], key=lambda t: t.getFitStatus().getChi2() / t.getFitStatus().getNdf())
                fs = ds_track.getFitStatus()
                chi2ndf = fs.getChi2() / fs.getNdf() + 1e-10
                fst = ds_track.getFittedState()
                pos = fst.getPos()
                mom = fst.getMom()
                self.track = ds_track
                self.pos, self.mom = pos, mom
                self.hasTrack = True

                Ex_ds = {}
                for plane in range(5):
                    zEx = self.zPos['MuFilter'][20 + plane]
                    lam = (zEx - pos.z()) / mom.z()
                    Ex_ds[plane] = ROOT.TVector3(pos.x() + lam * mom.x(), pos.y() + lam * mom.y(), pos.z() + lam * mom.z())
                self.Ex = Ex_ds

                if chi2ndf < 2 and self.slopecut() and self.fiducialvolumecut():
                    x = self.muAna.GetDSHaverage(hits)
                    if x != -999.:
                        reft_ds, fired_sipms = x

                        if 'reft-DS' not in self.hists:
                            self.hists['reft-DS'] = ROOT.TH1F('reft-DS', 'Average time of DS horizontal bars;DS horizontal average time [ns];Counts', 200, 0, 50)
                        self.hists['reft-DS'].Fill(reft_ds)

                        if 'NDSHSiPMs' not in self.hists:
                            self.hists['NDSHSiPMs'] = ROOT.TH1F('NDSHSiPMs', '#splitline{Number of SiPMs reading out horizontally}{aligned scintillators in the muon system};Number of SiPMs;Counts', 10, 0, 10)
                        self.hists['NDSHSiPMs'].Fill(fired_sipms)

                        if 'reducedchi2-DS' not in self.hists:
                            self.hists['reducedchi2-DS'] = ROOT.TH1F('reducedchi2-DS', '#chi^{2}_{#nu} (DS tracks);#chi^{2}_{#nu};Counts', 200, 0, 50)
                        self.hists['reducedchi2-DS'].Fill(chi2ndf)

                        if 'tds0vxpred' not in self.hists:
                            self.hists['tds0vxpred'] = ROOT.TH2F('tds0vxpred', 't_{0}^{DS} vs x-pos at US plane 5;x-position at US plane 5 [cm];t_{DS}^{0} [ns]', 100, -80, 20, 1000, 0, 50)
                        self.hists['tds0vxpred'].Fill(Ex_ds[4].x(), reft_ds)

                        if self.tds0cut(reft_ds):
                            contexts['DS'] = {'track': ds_track, 'pos': pos, 'mom': mom, 'Ex': Ex_ds, 'reft': reft_ds}

            if 'SF' in self.active_refsystems and tracks[1] and not self.simulation:
                sf_track = tracks[1][0] if len(tracks[1]) == 1 else min(tracks[1], key=lambda t: t.getFitStatus().getChi2() / t.getFitStatus().getNdf())
                fs = sf_track.getFitStatus()
                chi2ndf = fs.getChi2() / fs.getNdf() + 1e-10
                fst = sf_track.getFittedState()
                pos = fst.getPos()
                mom = fst.getMom()
                self.track = sf_track
                self.pos, self.mom = pos, mom
                self.hasTrack = True

                Ex_sf = {}
                for plane in range(5):
                    zEx = self.zPos['MuFilter'][20 + plane]
                    lam = (zEx - pos.z()) / mom.z()
                    Ex_sf[plane] = ROOT.TVector3(pos.x() + lam * mom.x(), pos.y() + lam * mom.y(), pos.z() + lam * mom.z())
                self.Ex = Ex_sf

                if chi2ndf < 2 and self.slopecut() and self.fiducialvolumecut():
                    reft_sf = self.muAna.GetScifiAverageTime(self.Scifi, scifi_hits)
                    if reft_sf is not None and reft_sf != -999.:
                        if 'reft-SF' not in self.hists:
                            self.hists['reft-SF'] = ROOT.TH1F('reft-SF', 'Average time of SiPMs in Scifi;Time [ns];Counts', 200, -25, 25)
                        self.hists['reft-SF'].Fill(reft_sf)

                        if 'reducedchi2-SF' not in self.hists:
                            self.hists['reducedchi2-SF'] = ROOT.TH1F('reducedchi2-SF', '#chi^{2}_{#nu} (SciFi tracks);#chi^{2}_{#nu};Counts', 200, 0, 50)
                        self.hists['reducedchi2-SF'].Fill(chi2ndf)

                        contexts['SF'] = {'track': sf_track, 'pos': pos, 'mom': mom, 'Ex': Ex_sf, 'reft': reft_sf}

            if not contexts:
                return

            for ref_name, ctx in contexts.items():
                self.track = ctx['track']
                self.pos = ctx['pos']
                self.mom = ctx['mom']
                self.hasTrack = True
                ctx['expected_bar'] = {plane: self.muAna.GetExtrapolatedBarDetID(plane) for plane in range(5)}

            for hit in event.Digi_MuFilterHits:
                if not hit.isValid(): continue
                detID = hit.GetDetectorID()
                s, p, b = self.muAna.parseDetID(detID)
                if s != 2: continue

                channels_t = hit.GetAllTimes()
                signals_dict = dict(hit.GetAllSignals())

                for ref_name, ctx in contexts.items():
                    if ctx['expected_bar'].get(p) != detID:
                        continue
                    self.track = ctx['track']
                    self.pos = ctx['pos']
                    self.mom = ctx['mom']
                    self.Ex = ctx['Ex']
                    self.reft = ctx['reft']
                    self.hasTrack = True

                    for channel in channels_t:
                        SiPM, clock = channel
                        if self.muAna.IsSmallSiPMchannel(SiPM): continue
                        if f'{detID}_{SiPM}' == '21005_15': continue
                        qdc = signals_dict.get(SiPM, -999.)
                        if qdc == -999.: continue
                        fixed_ch = f'{detID}_{SiPM}'

                        if self.step == 1:
                            self.zeroth(fixed_ch, clock, qdc, ref_name=ref_name)
                        elif self.step == 2:
                            self.tof(fixed_ch, clock, qdc, ref_name=ref_name)

            if self.step == 1:
                self.MakeExplantoryHist(event.Digi_MuFilterHits)
            return

        ### If there are more than 1 DS track, take the track with the lowest chi2/Ndf
        if len(tracks[self.referencesystem])==1: self.track = tracks[self.referencesystem][0]
        elif len(tracks[self.referencesystem])>1:
            tmp={i:i.getFitStatus().getChi2()/i.getFitStatus().getNdf() for i in tracks[self.referencesystem]}
            self.track=tmp[ min(tmp) ]

        self.hasTrack=True
        if (self.referencesystem==3 and not inDS) or (self.referencesystem==1 and not inScifi):
            self.hasTrack=False

        # Only throw event out if we need a track but don't have one
        if self.trackRequired and not self.hasTrack: return

        elif self.hasTrack:

            # cut flow hist
            if self.mode=='selectioncriteria': self.sc.cuteffect['has track']+=1

            fitStatus=self.track.getFitStatus()
            fstate=self.track.getFittedState()

            self.pos=fstate.getPos()
            self.mom=fstate.getMom()
            self.trackchi2NDF=fitStatus.getChi2()/fitStatus.getNdf()+1E-10
            self.chi2pndf_hist(fitStatus)

            self.Ex = {}
            for plane in range(5):
                zEx=self.zPos['MuFilter'][20+plane]
                lam=(zEx-self.pos.z())/self.mom.z()
                self.Ex[plane]=ROOT.TVector3(self.pos.x()+lam*self.mom.x(), self.pos.y()+lam*self.mom.y(), self.pos.z()+lam*self.mom.z())

            self.passfiducialvolumecut = self.fiducialvolumecut()

            # if self.mode=='selectioncriteria' and self.passfiducialvolumecut and self.hasTrack:
            #     self.sc.cuteffect['fiducial']+=1

            # cut flow hist
            self.passredchi2cut = False
            if self.trackchi2NDF < 2: self.passredchi2cut = True
            # if self.mode=='selectioncriteria' and self.passredchi2cut and self.hasTrack and self.passfiducialvolumecut: self.sc.cuteffect['#chi^{2}_{#nu}']+=1
            if self.mode=='selectioncriteria' and self.passredchi2cut and self.hasTrack: self.sc.cuteffect['#chi^{2}_{#nu}']+=1

            # Get DS event t0
            if self.referencesystem==3 and self.hasTrack and not self.simulation:
                x = self.muAna.GetDSHaverage(hits) # Now returns in ns
                if x==-999.: return
                self.reft, firedDSHSiPMs = x
                if not 'reft-DS' in self.hists:
                    self.hists['reft-DS']=ROOT.TH1F('reft','Average time of DS horizontal bars;DS horizontal average time [ns];Counts', 200, 0, 50)
                self.hists['reft-DS'].Fill(self.reft)

                if not 'NDSHSiPMs' in self.hists:
                    self.hists['NDSHSiPMs']=ROOT.TH1F('NDSHSiPMs','#splitline{Number of SiPMs reading out horizontally}{aligned scintillators in the muon system}; Number of SiPMs;Counts', 10, 0, 10)
                self.hists['NDSHSiPMs'].Fill(firedDSHSiPMs)

                if not self.tds0cut(self.reft): return

            elif self.referencesystem==1 and self.hasTrack and not self.simulation:
                self.reft = self.muAna.GetScifiAverageTime(self.Scifi, scifi_hits)
                if not self.reft: return
                if not 'reft-SF' in self.hists:
                    self.hists['reft-SF']=ROOT.TH1F('reft','Average time of SiPMs in Scifi;Time [ns];Counts', 200, -25, 25)
                self.hists['reft-SF'].Fill(self.reft)

            self.passslopecut = self.slopecut()
            # if self.mode=='selectioncriteria' and self.passslopecut and self.passredchi2cut and self.hasTrack and self.passfiducialvolumecut:
            if self.mode=='selectioncriteria' and self.passslopecut and self.passredchi2cut and self.hasTrack:
                self.sc.cuteffect['slope']+=1

        # Get DS event t0
        if self.referencesystem==3 and self.hasTrack and not self.simulation:

            t0vxpred=f'tds0vxpred'
            if not t0vxpred in self.hists:
                title='t_{0}^{DS} v x-position expected at HCAL plane 5;x-position at US plane 5 [cm];t_{DS}^{0} [ns]'
                self.hists[t0vxpred]=ROOT.TH2F(t0vxpred, title, 100, -80, 20, 1000, 0, 20)
            self.hists[t0vxpred].Fill(self.Ex[4].x(), self.reft)

        if self.mode == 'tds0-studies':
            self.tds0_studies(hits)
            return

        if self.mode == 'selectioncriteria':
            self.sc.Filltds0vcuts()
            self.sc.FillHists(hits, scifi_hits)
            return

        # Shower profiles also requires a DS track.
        elif self.mode=='showerprofiles':
            if self.simulation:
                mufilterpoints = event.MuFilterPoint
                hit2mc = event.Digi_MuFilterHits2MCPoints.At(0)
                # self.sp.Compare2MCPoints(hits, mufilterpoints, hit2mc)
                # self.hit_times = self.sp.CheckMCHitTimes(hits, mufilterpoints, hit2mc)
                self.sp.ShowerDirection(hits, scifi_hits, mufilterpoints=mufilterpoints, hit2mc=hit2mc)
                self.sp.FillMuonTrackEnergy()
            else:
                self.sp.ShowerDirection(hits, scifi_hits)

            self.sp.ExtractScifiData(scifi_hits)
            self.sp.EvaluateShowerSlopes() # Quite slow, comment if not needed
            self.sp.DoesBarFire(hits)
            return

        elif self.mode.find('extendedreconstruction')>-1:
            mufilterpoints = event.MuFilterPoint
            hit2mc = event.Digi_MuFilterHits2MCPoints.At(0)
            self.emr.ExtendReconstruction(hits, scifi_hits, mufilterpoints=mufilterpoints, hit2mc=hit2mc)
            return

        elif self.mode == 'struckquark':
            self.qve.StruckQuarkExtrapolation(hits)
            return

        elif self.mode == 'struckquark':
            self.qve.StruckQuarkExtrapolation(hits)
            return

        elif self.mode == 'reconstructmuonposition':
            self.sa.ReconstructMuonPosition(hits)
            return

        for hit in event.Digi_MuFilterHits:
            nLeft, nRight=self.muAna.GetnFiredSiPMs(hit)
            if not hit.isValid(): continue

            detID=hit.GetDetectorID()
            s,p,b=self.muAna.parseDetID(detID)
            if s!=2: continue

            # If the track does not go through this bar, skip it
            if self.muAna.GetExtrapolatedBarDetID(p) != detID: continue

            self.MuFilter.GetPosition(detID,self.A,self.B)

            # if self.yresidual3(detID): self.trackrelated = True
            # else: self.trackrelated=False

            channels_t=hit.GetAllTimes()
            channels_qdc=hit.GetAllSignals()

            if not all([self.passslopecut, self.fiducialvolumecut, self.passredchi2cut]): return

            if self.mode=='systemalignment' and s==2:
                self.sa.FillSiPMHists(hit)
                if self.options.XT: self.sa.XTHists(hit)
                self.sa.FillBarHists(hit)
                # if self.simulation:
                #     mufilterpoints = event.MuFilterPoint
                #     hit2mc = event.Digi_MuFilterHits2MCPoints.At(0)
                    # self.sp.Compare2MCPoints(hits, mufilterpoints, hit2mc)
                continue

            if self.mode in ('zeroth', 'full'): self.MakeExplantoryHist(event.Digi_MuFilterHits)

            for channel in channels_t:
                SiPM,clock=channel
                if self.muAna.IsSmallSiPMchannel(SiPM): continue
                if f'{detID}_{SiPM}' == '21005_15': continue # Dead SiPM
                qdc=self.muAna.GetChannelVal(SiPM, channels_qdc)
                if qdc==-999.: continue
                fixed_ch=f'{detID}_{SiPM}'
                if self.mode in ('zeroth', 'full'): self.zeroth(fixed_ch, clock, qdc)
                elif self.mode=='tof': self.tof(fixed_ch, clock, qdc)
                elif self.mode=='tw': self.tw(fixed_ch, clock, qdc, meantimecorrection=False)
                elif self.mode=='res': self.res(fixed_ch, clock, qdc)

    def chi2pndf_hist(self, fitStatus):
        histname = 'reducedchi2'
        if not histname in self.hists:
            if self.referencesystem==3: tmp='muon system'
            else: tmp='scifi'
            title = '#chi^{2} per degree of freedom for tracks in '+tmp+';#chi^{2}_{#nu};Counts'
            self.hists[histname] = ROOT.TH1F(histname, title, 200, 0, 50)
        self.hists[histname].Fill(self.trackchi2NDF)

        histname = 'chi2'
        if not histname in self.hists:
            if self.referencesystem==3: tmp='muon system'
            else: tmp='scifi'
            title = '#chi^{2} for tracks in '+tmp+';#chi^{2} [unitless];Counts'
            self.hists[histname] = ROOT.TH1F(histname, title, 100, 0, 100)
        self.hists[histname].Fill(fitStatus.getChi2())

        histname = 'dof'
        if not histname in self.hists:
            if self.referencesystem==3: tmp='muon system'
            else: tmp='scifi'
            title = 'degrees of freedom for tracks in '+tmp+';Degrees of freedom [unitless];Counts'
            self.hists[histname] = ROOT.TH1F(histname, title, 5, 0, 5)
        self.hists[histname].Fill(fitStatus.getNdf())

    def FillChannelRateHists(self):
        hits=self.M.eventTree.Digi_MuFilterHits
        OneHitPerUS=self.muAna.OneHitPerSystem(hits, (2,))
        OneHitPerDS=self.muAna.OneHitPerSystem(hits, (3,))
        for hit in hits:
            detID=hit.GetDetectorID()
            s,p,b=self.muAna.parseDetID(detID)

            if s==3: continue
            for x in hit.GetAllTimes():
                SiPM, clock = x
                channel=self.muAna.GetSiPMNumberInSystem_PCBbyPCB(detID, SiPM)
                if OneHitPerUS: self.hists[f'channelhitrate-{self.subsystemdict[s]}-1USb/p'].Fill(channel)
                if OneHitPerDS: self.hists[f'channelhitrate-{self.subsystemdict[s]}-1DSb/p'].Fill(channel)
                self.hists[f'channelhitrate-{self.subsystemdict[s]}'].Fill(channel)
                if self.yresidual3(detID): self.hists[f'channelhitrate-{self.subsystemdict[s]}-dycut'].Fill(channel)

    def fiducialvolumecut(self):
        # Just computing the distance of the track from the sides in the first plane
        sides_residuals={'left':True, 'right':True}

        for plane in range(5):

            self.MuFilter.GetPosition(int(f'2{plane}004'), self.A, self.B)

            for side in ('left', 'right'):
                # Do not accept hits where the muon extrapolates to within SiPM_distance_cut cm of the bar end
                if side == 'left' and (self.A.x() - self.Ex[plane].x()) < SiPM_distance_cut:
                    return False
                elif side == 'right' and (self.Ex[plane].x() - self.B.x()) < SiPM_distance_cut:
                    return False
            return True

    def zeroth(self, fixed_ch, clock, qdc, ref_name='DS'):

        ReadableFixedCh = self.muAna.MakeHumanReadableFixedCh(fixed_ch)

        detID, SiPM = (int(fixed_ch.split('_')[i]) for i in range(2))
        s, p, b = self.muAna.parseDetID(detID)

        correctedtime = clock * self.TDC2ns
        t_rel = self.reft - correctedtime

        dtvpred = f'dtvxpred_{ref_name}_{fixed_ch}_{self.state}'

        if not dtvpred in self.hists:
            title = f'{ReadableFixedCh} No TW corr t_{{0}}^{{{ref_name}}} - t_{{SiPM}} vs x-position ({ref_name});x_{{predicted}} [cm];t_{{0}}^{{{ref_name}}} - t_{{SiPM}} [ns]'
            histformat = self.get_histformat("dtvxpred", self.state, ref_name)
            self.hists[dtvpred] = ROOT.TH2F(dtvpred, title, *histformat[0], *histformat[1])

        self.hists[dtvpred].Fill(self.Ex[p].x(), t_rel)

    def MakeExplantoryHist(self, hits):

        # print(f'Running makeExplanatoryHist')
        # Make a QDC distribution for hits in each third of the bar

        for hit in hits:
            detID=hit.GetDetectorID()
            if detID != 24004: continue

            s,p,b=self.muAna.parseDetID(detID)
            if not s==2: continue

            for x in hit.GetAllSignals():
                SiPM, qdc = x
                if SiPM!=4: continue # Only interested in 24004_4 for this histogram
                self.MuFilter.GetPosition(detID, self.A, self.B)
                bar_L = (self.A.x() - self.B.x())

                relative_bar_pos = (self.Ex[p].x() - self.B.x()) / bar_L # fraction of bar length
                if relative_bar_pos > 0 and relative_bar_pos <= 1/3: bar_third=1
                elif relative_bar_pos > 1/3 and relative_bar_pos <= 2/3: bar_third=2
                elif relative_bar_pos > 2/3 and relative_bar_pos <= 1: bar_third=3
                else:
                    print(f'Warning: relative_bar_pos={relative_bar_pos:.2f} outside expected range [0,1]')
                    continue

                histname = f'qdcthird{bar_third}_{detID}_{SiPM}'
                if not histname in self.hists:
                    title=f'QDC distribution for hits {bar_third}/3 from the right side;QDC [a.u];Counts'
                    self.hists[histname] = ROOT.TH1F(histname, title, 200, 0, 200)
                self.hists[histname].Fill(qdc)

    def tof(self, fixed_ch, clock, qdc, ref_name='DS'):

        ReadableFixedCh = self.muAna.MakeHumanReadableFixedCh(fixed_ch)

        detID, SiPM = fixed_ch.split('_')
        s, p, b = self.muAna.parseDetID(int(detID))

        cscint_val = None
        if hasattr(self, 'cscintvalues_by_ref') and ref_name in self.cscintvalues_by_ref:
            cscint_val = self.cscintvalues_by_ref[ref_name].get(fixed_ch)
        if not cscint_val and hasattr(self, 'cscintvalues_by_ref') and 'DS' in self.cscintvalues_by_ref:
            cscint_val = self.cscintvalues_by_ref['DS'].get(fixed_ch)
        if not cscint_val and hasattr(self.muAna, 'cscintvalues') and self.muAna.cscintvalues:
            cscint_val = self.muAna.cscintvalues.get(fixed_ch)

        if not cscint_val:
            return

        _, ToFcorrectedtime = self.muAna.correct_ToF(self.MuFilter, fixed_ch, clock, self.Ex[p].x(), cscint_val=cscint_val)
        if ToFcorrectedtime is None: return

        t_rel = self.reft - ToFcorrectedtime
        dtvqdc = f'dtvqdc_{ref_name}_{fixed_ch}_{self.state}'

        if not dtvqdc in self.hists:
            title = f'{ReadableFixedCh} ToF-corrected time vs QDC ({ref_name});QDC_{{SiPM}} [a.u.];dt^{{ToF}}_{{SiPM}} [ns]'
            histformat = self.get_histformat("dtvqdc", self.state, ref_name)
            self.hists[dtvqdc] = ROOT.TH2F(dtvqdc, title, *histformat[0], *histformat[1])

        self.hists[dtvqdc].Fill(qdc, t_rel)

    def tw(self, fixed_ch, clock, qdc, meantimecorrection=False):

        hists=self.hists
        ReadableFixedCh=self.muAna.MakeHumanReadableFixedCh(fixed_ch)

        detID, SiPM=int(fixed_ch.split('_')[0]), int(fixed_ch.split('_')[1])
        s,p,b = self.muAna.parseDetID(detID)

        # ### Times wrt to DS horizontal average
        TWcorrectedtime = self.muAna.MuFilterCorrectedTime(self.MuFilter, fixed_ch, clock, qdc, x=0)
        if not TWcorrectedtime: return
        TWt_rel = self.reft - TWcorrectedtime

        ### Make histograms
        dtvxpred=f'dtvxpred_{fixed_ch}_{self.state}'
        if not dtvxpred in hists and self.referencesystem==3:
            subtitle='{Time-walk corrected t_{0}^{DS}-t^{tw corr}_{SiPM} v x-position w/ run'+str(self.TWCorrectionRun)+'};x_{predicted} [cm];t_{0}^{DS}-t^{tw corr}_{SiPM} [ns]'
            title='#splitline{'+ReadableFixedCh+'}'+subtitle
            histformat = self.get_histformat('dtvxpred', self.state, 'DS')
            hists[dtvxpred]=ROOT.TH2F(dtvxpred,title,*histformat[0], *histformat[1])

        hists[dtvxpred].Fill(self.Ex[p].x(), TWt_rel)

    def res(self, fixed_ch, clock, qdc):

        hists=self.hists
        ReadableFixedCh=self.muAna.MakeHumanReadableFixedCh(fixed_ch)
        detID, SiPM = int(fixed_ch.split('_')[0]), int(fixed_ch.split('_')[1])
        s,p,b=self.muAna.parseDetID(detID)

        ToFTWtime = self.muAna.MuFilterCorrectedTime(self.MuFilter, fixed_ch, clock, qdc, x=self.Ex[p].x())
        if not ToFTWtime: return
        ToFTWt_rel = self.reft - ToFTWtime

        ### Make histograms
        if self.referencesystem==3: dtvqdc=f'dtvqdc_{fixed_ch}_corrected'
        elif self.referencesystem==1: dtvqdc=f'dtSFvqdc_{fixed_ch}_corrected'

        if not dtvqdc in hists:
            title=ReadableFixedCh+'. dt_{SiPM}^{TW ToF} v QDC_{SiPM};QDC_{SiPM} [a.u];dt_{SiPM}^{TW ToF} [ns]'
            # title='#splitline{'+ReadableFixedCh+'}'+subtitle
            ref_res = 'SF' if self.referencesystem == 1 else 'DS'
            histformat = self.get_histformat('dtvqdc', self.state, ref_res)
            hists[dtvqdc]=ROOT.TH2F(dtvqdc, title, *histformat[0], *histformat[1])

        ### Fill histograms
        hists[dtvqdc].Fill(qdc, ToFTWt_rel)

    def tds0_studies(self, hits):

        if self.passslopecut:
            name=f'tds0-slopecut'
            if not name in self.hists:
                title = 't_{0}^{ds} with slope cut;t_{0}^{ds} [ns];Counts'
                self.hists[name]=ROOT.TH1F(name, title, 200, 0, 50)
            self.hists[name].Fill(self.reft)

        if getattr(self, 'passtdcut', False):
            name=f'tds0-tdcut'
            if not name in self.hists:
                title = 't_{0}^{ds} with timing discriminant cut;t_{0}^{ds} [ns];Counts'
                self.hists[name]=ROOT.TH1F(name, title, 200, 0, 50)
            self.hists[name].Fill(self.reft)

        if getattr(self, 'passtdcut', False) and self.passslopecut:
            name=f'tds0-tdcut+slopecut'
            if not name in self.hists:
                title = 't_{0}^{ds} with timing discriminant cut and slope cut;t_{0}^{ds} [ns];Counts'
                self.hists[name]=ROOT.TH1F(name, title, 200, 0, 50)
            self.hists[name].Fill(self.reft)

        ### Apply slope cut to all these to minimise non-IP1 contribution
        if not self.passslopecut: return

        for mode in ('testing-tds0','testing-deltastations'):
            testing_value = self.muAna.GetDSHaverage(hits, mode)

            if mode.find('tds0')!=-1:

                if testing_value==-999.: continue
                tds0, firedDSHbars = testing_value

                if not mode in self.hists:
                    title = f'{mode.replace("-", " ")}; [ns];Counts'
                    bins = (200, 0, 50)
                    self.hists[mode] = ROOT.TH1F(mode, title, *bins)
                self.hists[mode].Fill(tds0)

            else:
                if isinstance(testing_value, bool): continue
                bins = (100, -12.5, 12.5)
                for d in ('32', '21'):
                    name = f'{mode}{d}'
                    if not name in self.hists:
                        title = f'{mode.replace("-", " ")} {d}; [ns];Counts'
                        self.hists[name] = ROOT.TH1F(name, title, *bins)
                    self.hists[name].Fill(testing_value[f'delta{d}'])

    def Scifi_t0_studies(self, scifi_hits):
        name = f'scifi_averagetime'
        if not name in self.hists:
            title = f'{mode.replace("-", " ")} {d}; [ns];Counts'
            self.hists[name] = ROOT.TH1F(name, title, *bins)
        self.hists[name].Fill(testing_value[f'delta{d}'])

    def GetSimEngine(self):
        simEngine = self.options.geoFile.split('.')[1].split('-')[0]
        return simEngine

    def WriteOutHistograms(self):

        if self.mode == 'calibration' or self.mode in ('zeroth', 'tof', 'full', 'tw', 'res'):
            if getattr(self.options, 'outFile', None):
                outfilename = self.options.outFile
            else:
                run_tag = f'run{str(self.runNr).zfill(6)}'
                base_dir = self.outputDir
                if not (base_dir.endswith(run_tag) or base_dir.endswith(f'run{self.runNr}')):
                    base_dir = os.path.join(base_dir, run_tag)
                outfilename = os.path.join(base_dir, 'timewalk_calibration.root')

            os.makedirs(os.path.dirname(os.path.abspath(outfilename)), exist_ok=True)
            print(f"[INFO] Opening output ROOT file: {outfilename}", flush=True)
            f = ROOT.TFile.Open(outfilename, 'recreate')

            # Global run/reference histograms written at the root level of the file
            global_written = set()
            global_hist_keys = [
                'reft-DS', 'reft-SF', 'NDSHSiPMs', 'reducedchi2-DS', 'reducedchi2-SF', 'reducedchi2', 'tds0vxpred'
            ]
            for gkey in global_hist_keys:
                if gkey in self.hists:
                    f.WriteObject(self.hists[gkey], gkey, 'kOverwrite')
                    global_written.add(gkey)

            if 'reft-DS' in self.hists:
                h_ds = self.hists['reft-DS']
                c_reft_ds = ROOT.TCanvas('c_reft_DS', f'DS Reference Time Distribution with Selection Window (Run {self.runNr})', 1000, 750)
                c_reft_ds.SetGrid()
                c_reft_ds.SetLogy(1)
                h_ds.SetLineWidth(2)
                h_ds.SetLineColor(ROOT.kBlack)
                h_ds.Draw('HIST')

                max_bin = h_ds.GetMaximumBin()
                mode_t = h_ds.GetBinCenter(max_bin)
                if self.tds0_min is None or self.tds0_max is None:
                    self.tds0_mode = mode_t
                    self.tds0_min = mode_t - 1.75
                    self.tds0_max = mode_t + 2.15
                elif self.tds0_mode is None:
                    self.tds0_mode = mode_t

                y_max = max(10.0, h_ds.GetMaximum() * 2.5)
                y_min = 0.5

                box_ds = ROOT.TBox(self.tds0_min, y_min, self.tds0_max, y_max)
                box_ds.SetFillColorAlpha(ROOT.kGreen + 2, 0.20)
                box_ds.SetLineColor(ROOT.kGreen + 3)
                box_ds.SetLineWidth(1)
                box_ds.Draw('SAME')

                l_min = ROOT.TLine(self.tds0_min, y_min, self.tds0_min, y_max)
                l_min.SetLineColor(ROOT.kRed + 1)
                l_min.SetLineWidth(2)
                l_min.SetLineStyle(2)
                l_min.Draw('SAME')

                l_max = ROOT.TLine(self.tds0_max, y_min, self.tds0_max, y_max)
                l_max.SetLineColor(ROOT.kRed + 1)
                l_max.SetLineWidth(2)
                l_max.SetLineStyle(2)
                l_max.Draw('SAME')

                h_ds.Draw('HIST SAME')
                l_min.Draw('SAME')
                l_max.Draw('SAME')

                h_ds.GetListOfFunctions().Add(l_min)
                h_ds.GetListOfFunctions().Add(l_max)

                n_total = h_ds.GetEntries()
                bin_lo = h_ds.FindBin(self.tds0_min)
                bin_hi = h_ds.FindBin(self.tds0_max)
                n_pass = h_ds.Integral(bin_lo, bin_hi)
                pct_pass = (n_pass / n_total * 100.0) if n_total > 0 else 0.0

                pave = ROOT.TPaveText(0.55, 0.70, 0.88, 0.88, 'NDC')
                pave.SetFillColorAlpha(ROOT.kWhite, 0.85)
                pave.SetBorderSize(1)
                pave.SetTextAlign(12)
                pave.SetTextSize(0.028)
                pave.AddText(f'Run: {self.runNr} ({getattr(self, "timealignment", "unknown")} era)')
                if self.tds0_mode is not None:
                    pave.AddText(f'Mode: {self.tds0_mode:.2f} ns')
                pave.AddText(f't_{{0}}^{{DS}} Window: [{self.tds0_min:.2f}, {self.tds0_max:.2f}] ns')
                pave.AddText(f'Accepted: {int(n_pass):,} / {int(n_total):,} ({pct_pass:.1f}%)')
                pave.Draw('SAME')
                c_reft_ds.Update()

                f.WriteObject(c_reft_ds, 'c_reft_DS', 'kOverwrite')
                global_written.add('c_reft_DS')

                sum_dir = f.GetDirectory('Summary')
                if not sum_dir:
                    sum_dir = f.mkdir('Summary')
                ds_sum_dir = sum_dir.GetDirectory('DS')
                if not ds_sum_dir:
                    ds_sum_dir = sum_dir.mkdir('DS')
                ds_sum_dir.WriteObject(c_reft_ds, 'c_reft_DS', 'kOverwrite')

            if 'reft-SF' in self.hists:
                h_sf = self.hists['reft-SF']
                c_reft_sf = ROOT.TCanvas('c_reft_SF', f'SF Reference Time Distribution (Run {self.runNr})', 1000, 750)
                c_reft_sf.SetGrid()
                c_reft_sf.SetLogy(1)
                h_sf.SetLineWidth(2)
                h_sf.SetLineColor(ROOT.kBlack)
                h_sf.Draw('HIST')
                c_reft_sf.Update()

                f.WriteObject(c_reft_sf, 'c_reft_SF', 'kOverwrite')
                global_written.add('c_reft_SF')

                sum_dir = f.GetDirectory('Summary')
                if not sum_dir:
                    sum_dir = f.mkdir('Summary')
                sf_sum_dir = sum_dir.GetDirectory('SF')
                if not sf_sum_dir:
                    sf_sum_dir = sum_dir.mkdir('SF')
                sf_sum_dir.WriteObject(c_reft_sf, 'c_reft_SF', 'kOverwrite')

            for bar_third in (1, 2, 3):
                histname = f'qdcthird{bar_third}_24004_4'
                if histname in self.hists:
                    d_exp = f.GetDirectory('24004_4')
                    if not d_exp:
                        d_exp = f.mkdir('24004_4')
                    d_exp.WriteObject(self.hists[histname], histname, 'kOverwrite')
                    global_written.add(histname)

            print(f"[INFO] Writing {len(self.hists)} 2D histograms into channel directories...", flush=True)
            # Write channel-specific histograms into their dedicated TDirectory (<detID>_<SiPM>)
            for h in self.hists:
                if h in global_written:
                    continue

                parts = h.split('_')
                if len(parts) >= 5 and parts[1] in ('DS', 'SF'):
                    histkey, ref, detID, SiPM = parts[0], parts[1], parts[2], parts[3]
                    fixed_ch = f'{detID}_{SiPM}'
                elif len(parts) == 4:
                    histkey, detID, SiPM, state = parts
                    fixed_ch = f'{detID}_{SiPM}'
                elif len(parts) >= 2:
                    fixed_ch = '_'.join(parts[1:])
                else:
                    # Write unclassified histograms to root level
                    f.WriteObject(self.hists[h], h, 'kOverwrite')
                    continue

                ch_dir = f.GetDirectory(fixed_ch)
                if not ch_dir:
                    ch_dir = f.mkdir(fixed_ch)
                ch_dir.WriteObject(self.hists[h], self.hists[h].GetName(), 'kOverwrite')

            if getattr(self.options, 'speedOnly', False):
                print("[INFO] Processing Phase 1/1: c_scint speed fits & overview canvases (--speedOnly)...", flush=True)
                self.FitAndSaveCscint(f)
                print("[INFO] Finalizing and closing ROOT file...", flush=True)
                f.Close()
                print(f"[SUCCESS] c_scint calibration complete. Results written to {outfilename}", flush=True)
                return

            print("[INFO] Processing Phase 1/2: c_scint speed fits & overview canvases...", flush=True)
            self.FitAndSaveCscint(f)

            print("[INFO] Processing Phase 2/2: 5-parameter time-walk rational model fits & canvases...", flush=True)
            self.FitAndSaveTimeWalk(f)

            print("[INFO] Finalizing and closing ROOT file...", flush=True)
            f.Close()
            print(f"[SUCCESS] Calibration complete. Results written to {outfilename}", flush=True)

        elif self.mode == 'systemalignment':
            self.sa.WriteOutHistograms()

        elif self.mode == 'showerprofiles':
            self.sp.WriteOutHistograms()

        elif self.mode == 'selectioncriteria':
            self.sc.WriteOutHistograms()

        elif self.mode == 'reconstructmuonposition':
            self.sa.WriteOutHistograms()

        elif self.mode.find('extendedreconstruction')>-1:
            self.emr.WriteOutHistograms()

        elif self.mode.find('struckquark')>-1:
            self.qve.WriteOutHistograms()

        elif self.mode == 'tds0-studies':
            if getattr(self.options, 'outFile', None):
                outfilename = self.options.outFile
            else:
                outfilename = f'{self.outpath}/splitfiles/run{self.runNr}/tds0-studies.root'
            os.makedirs(os.path.dirname(os.path.abspath(outfilename)), exist_ok=True)
            f = ROOT.TFile.Open(outfilename, 'recreate')
            for h in self.hists:
                hist = self.hists[h]
                f.WriteObject(hist, hist.GetName(), 'kOverwrite')
            print(f'Hists written to {outfilename}')
            f.Close()

    def SaveRawHistograms(self, outfilename=None):
        if not outfilename:
            out_dir = getattr(self.options, 'outputDir', './Timing')
            outfilename = getattr(self.options, 'outFile', None) or os.path.join(out_dir, f'raw_histograms_{self.options.nStart}.root')
        os.makedirs(os.path.dirname(os.path.abspath(outfilename)), exist_ok=True)
        print(f"[INFO] Saving {len(self.hists)} raw histograms to {outfilename}...", flush=True)
        f = ROOT.TFile.Open(outfilename, 'recreate')
        for hname, hist in self.hists.items():
            if hist:
                f.WriteObject(hist, hname, 'kOverwrite')
        f.Close()
        print(f"[SUCCESS] Raw histograms saved to {outfilename}", flush=True)

    def LoadHistogramsFromRoot(self, root_filename):
        if not os.path.exists(root_filename):
            print(f"[ERROR] Histogram file {root_filename} does not exist.")
            return False
        f = ROOT.TFile.Open(root_filename)
        if not f or f.IsZombie():
            print(f"[ERROR] Could not open {root_filename}.")
            return False

        count = 0
        def read_directory(tdir):
            nonlocal count
            for key in tdir.GetListOfKeys():
                kname = key.GetName()
                obj = key.ReadObj()
                if isinstance(obj, ROOT.TH1):
                    obj.SetDirectory(0)
                    self.hists[kname] = obj
                    count += 1
                elif isinstance(obj, ROOT.TDirectory):
                    read_directory(obj)

        read_directory(f)
        f.Close()
        print(f"[INFO] Successfully loaded {count} histograms from {root_filename}.", flush=True)
        return True

    def LoadCscintFile(self, cscint_path):
        if not cscint_path or not os.path.exists(cscint_path):
            print(f"Warning: cscint file {cscint_path} does not exist.")
            return False
        with open(cscint_path, 'r') as f:
            data = json.load(f)

        if 'tds0_window' in data and isinstance(data['tds0_window'], dict):
            win = data['tds0_window']
            if getattr(self.options, 'tds0Min', None) is None:
                self.tds0_min = win.get('min', self.tds0_min)
            if getattr(self.options, 'tds0Max', None) is None:
                self.tds0_max = win.get('max', self.tds0_max)
            self.tds0_mode = win.get('mode', None)
            print(f"[INFO] Loaded dynamic tds0 window from {cscint_path}: mode={self.tds0_mode}, [{self.tds0_min:.2f}, {self.tds0_max:.2f}] ns")

        self.cscintvalues_by_ref = {'DS': {}, 'SF': {}}
        if 'DS' in data or 'SF' in data:
            if 'DS' in data and isinstance(data['DS'], dict):
                self.cscintvalues_by_ref['DS'] = data['DS']
            if 'SF' in data and isinstance(data['SF'], dict):
                self.cscintvalues_by_ref['SF'] = data['SF']
            if not self.cscintvalues_by_ref['SF'] and self.cscintvalues_by_ref['DS']:
                print(f"[INFO] c_scint constants for SF empty; falling back to DS bar speeds ({len(self.cscintvalues_by_ref['DS'])} channels).")
                self.cscintvalues_by_ref['SF'] = dict(self.cscintvalues_by_ref['DS'])
        elif 'cscint' in data and isinstance(data['cscint'], dict):
            self.cscintvalues_by_ref['DS'] = data['cscint']
            self.cscintvalues_by_ref['SF'] = data['cscint']
        else:
            self.cscintvalues_by_ref['DS'] = data
            self.cscintvalues_by_ref['SF'] = data

        print(f"Loaded c_scint constants from {cscint_path}: DS={len(self.cscintvalues_by_ref['DS'])} channels, SF={len(self.cscintvalues_by_ref['SF'])} channels.")
        return True

    def FitAndSaveCscint(self, root_file):
        prev_ignore = ROOT.gErrorIgnoreLevel
        ROOT.gErrorIgnoreLevel = ROOT.kError
        min_entries_ds = getattr(self.options, 'minEntries', 1000)
        min_entries_sf = getattr(self.options, 'minEntriesSF', min(150, min_entries_ds))
        fit_min_ds = getattr(self.options, 'fitMinDS', getattr(self.options, 'fitMin', -55.0))
        fit_max_ds = getattr(self.options, 'fitMaxDS', getattr(self.options, 'fitMax', -5.0))
        fit_min_sf = getattr(self.options, 'fitMinSF', -44.0)
        fit_max_sf = getattr(self.options, 'fitMaxSF', -8.0)
        min_slice_entries = getattr(self.options, 'minSliceEntries', 15)
        n_segments = getattr(self.options, 'nSegments', 5)
        segment_overlap = getattr(self.options, 'segmentOverlap', 0.85)

        def meets_cscint_thresh(hname, h2):
            if not hname.startswith('dtvxpred_'):
                return False
            is_sf = '_SF_' in hname
            th = min_entries_sf if is_sf else min_entries_ds
            return h2.GetEntries() >= th

        summary_cscint = {'DS': {}, 'SF': {}}
        qa_cscint = {
            'DS': {'slice_chi2': [], 'fit_chi2': [], 'speeds': [], 'speed_errs': [], 'ch_indices': []},
            'SF': {'slice_chi2': [], 'fit_chi2': [], 'speeds': [], 'speed_errs': [], 'ch_indices': []}
        }
        save_slices = getattr(self.options, 'saveSlices', False)
        target_cscint_hists = [h for h in sorted(self.hists.items()) if meets_cscint_thresh(h[0], h[1])]
        total_cscint = len(target_cscint_hists)

        for idx, (hname, h2) in enumerate(target_cscint_hists):
            if (idx + 1) % 50 == 0 or (idx + 1) == total_cscint:
                print(f"  [Progress c_scint] Processed {idx + 1}/{total_cscint} channels ({(idx + 1) / total_cscint * 100:.0f}%)...", flush=True)

            parts = hname.split('_')
            # format: dtvxpred_<ref>_<detID>_<SiPM>_<state> or dtvxpred_<detID>_<SiPM>_<state>
            if len(parts) >= 5 and parts[1] in ('DS', 'SF'):
                ref = parts[1]
                detID = int(parts[2])
                sipm = int(parts[3])
            elif len(parts) >= 4:
                ref = 'DS'
                detID = int(parts[1])
                sipm = int(parts[2])
            else:
                continue
            fixed_ch = f'{detID}_{sipm}'
            fit_min = fit_min_sf if ref == 'SF' else fit_min_ds
            fit_max = fit_max_sf if ref == 'SF' else fit_max_ds

            if self.muAna.IsSmallSiPMchannel(sipm) and not getattr(self.options, 'fitSmallSiPMs', False):
                continue
            if fixed_ch == '21005_15':
                continue

            gr = ROOT.TGraphErrors()
            gr.SetName(f'gr_mpv_{ref}_{fixed_ch}')
            gr.SetTitle(f'MPV vs x_{{pred}} {fixed_ch} ({ref} Reference);x_{{pred}} [cm];t_{{0}}^{{{ref}}} - t_{{SiPM}} [ns]')

            slice_fits = []
            ip = 0
            for bx in range(1, h2.GetNbinsX() + 1):
                x = h2.GetXaxis().GetBinCenter(bx)
                if x < fit_min or x > fit_max:
                    continue
                py = h2.ProjectionY(f'h_slice_{ref}_{fixed_ch}_bx{bx}', bx, bx)
                if py.GetEntries() < min_slice_entries:
                    continue

                max_bin = py.GetMaximumBin()
                mode = py.GetBinCenter(max_bin)
                rms = py.GetRMS()
                fit_range = max(0.8, min(1.8, 1.5 * rms))

                fg = ROOT.TF1(f'fg_{ref}_{fixed_ch}_bx{bx}', 'gaus', mode - fit_range, mode + fit_range)
                fg.SetParameters(py.GetBinContent(max_bin), mode, max(0.3, rms))
                res_g = py.Fit(fg, 'QNS0')
                mpv = None
                mpv_err = None
                fit_func_to_save = None
                fit_type = 'mode'
                if int(res_g) == 0 and 0 < fg.GetParError(1) < 2.5:
                    p1 = fg.GetParameter(1)
                    if mode - 2.5 < p1 < mode + 2.5:
                        mpv = p1
                        mpv_err = fg.GetParError(1)
                        fit_func_to_save = fg
                        fit_type = 'gaus'

                if mpv is None:
                    mpv = mode
                    mpv_err = py.GetBinWidth(max_bin) / (12.0 ** 0.5)

                if fit_func_to_save is not None:
                    ndf_slice = fit_func_to_save.GetNDF()
                    if ndf_slice > 0:
                        qa_cscint[ref]['slice_chi2'].append(fit_func_to_save.GetChisquare() / ndf_slice)

                gr.SetPoint(ip, x, mpv)
                gr.SetPointError(ip, 0.0, mpv_err)
                ip += 1
                slice_fits.append((bx, x, py, fit_func_to_save, mpv, mpv_err, fit_type))

            if gr.GetN() < 5:
                continue

            func = ROOT.TF1(f'fit_cscint_{ref}_{fixed_ch}', 'pol1', fit_min, fit_max)
            func.SetLineColor(ROOT.kRed)
            func.SetLineWidth(2)
            gr.Fit(func, 'QNS0R')

            slope = func.GetParameter(1)
            slope_err = func.GetParError(1)
            const = func.GetParameter(0)
            const_err = func.GetParError(0)
            chi2 = func.GetChisquare()
            ndf = max(1, func.GetNDF())

            if abs(slope) < 1e-7 or slope_err <= 0:
                continue

            side_sign = 1.0 if sipm < 8 else -1.0
            c_scint = side_sign / slope
            delta_c_fit = abs(slope_err / (slope ** 2))

            if c_scint <= 0 or c_scint > 50.0:
                continue

            # Systematic speed uncertainty via overlapping sub-segments (default 5 segments with 85% overlap)
            total_span = fit_max - fit_min
            step_frac = 1.0 - segment_overlap
            denom = 1.0 + (n_segments - 1) * step_frac
            seg_width = total_span / denom if denom > 0 else total_span
            seg_step = seg_width * step_frac

            subrange_speeds = []
            for s_idx in range(n_segments):
                x_low = fit_min + s_idx * seg_step
                x_high = min(fit_max, x_low + seg_width)
                n_pts_in_seg = sum(1 for i in range(gr.GetN()) if x_low <= gr.GetPointX(i) <= x_high)
                if n_pts_in_seg < 3:
                    continue
                seg_func = ROOT.TF1(f'seg_{ref}_{fixed_ch}_{s_idx}', 'pol1', x_low, x_high)
                gr.Fit(seg_func, 'QNS0R')
                m_s = seg_func.GetParameter(1)
                if abs(m_s) > 1e-7:
                    c_s = side_sign / m_s
                    if 5.0 < c_s < 40.0:
                        subrange_speeds.append(float(c_s))

            if len(subrange_speeds) >= 2:
                delta_c_syst = float(np.std(subrange_speeds, ddof=1))
            else:
                delta_c_syst = float(delta_c_fit)

            cname = f'c_cscint_{ref}_{fixed_ch}'
            canv = ROOT.TCanvas(cname, f'Effective Speed Fit {fixed_ch} ({ref} Reference)', 900, 700)
            canv.SetGrid()
            h2.Draw('COLZ')

            gr.SetMarkerStyle(20)
            gr.SetMarkerSize(0.7)
            gr.SetMarkerColor(ROOT.kBlack)
            gr.SetLineColor(ROOT.kBlack)
            gr.Draw('P SAME')

            func.Draw('SAME')

            pave = ROOT.TPaveText(0.15, 0.72, 0.52, 0.88, 'NDC')
            pave.SetFillColorAlpha(ROOT.kWhite, 0.85)
            pave.SetBorderSize(1)
            pave.SetTextAlign(12)
            pave.SetTextSize(0.03)
            pave.AddText(f'Channel: {fixed_ch} ({ref} Reference)')
            pave.AddText(f'c_{{scint}} = {c_scint:.2f} #pm {delta_c_fit:.2f} (fit) #pm {delta_c_syst:.2f} (syst) cm/ns')
            pave.AddText(f'#chi^{{2}}/ndf = {chi2:.1f}/{ndf} ({chi2/ndf:.2f})')
            pave.Draw('SAME')
            canv.Update()

            if root_file is not None:
                ch_dir = root_file.GetDirectory(fixed_ch)
                if not ch_dir:
                    ch_dir = root_file.mkdir(fixed_ch)
                ch_dir.WriteObject(canv, cname, 'kOverwrite')
                ch_dir.WriteObject(gr, gr.GetName(), 'kOverwrite')
                ch_dir.WriteObject(func, func.GetName(), 'kOverwrite')

                prof_x = h2.ProfileX(f'prof_{ref}_{fixed_ch}')
                prof_x.SetTitle(f'ProfileX {fixed_ch} ({ref} Reference);x_{{pred}} [cm];<t_{{0}}^{{{ref}}} - t_{{SiPM}}> [ns]')
                ch_dir.WriteObject(prof_x, prof_x.GetName(), 'kOverwrite')

                canv.Close()
                del canv

                if save_slices:
                    sub_dir = ch_dir.GetDirectory(f'slice_fits_{ref}')
                    if not sub_dir:
                        sub_dir = ch_dir.mkdir(f'slice_fits_{ref}')

                    for bx_idx, x_val, py_hist, f_fit, m_val, m_err, fit_type in slice_fits:
                        if f_fit:
                            f_fit.SetLineColor(ROOT.kRed + 1)
                            f_fit.SetLineWidth(2)
                            py_hist.GetListOfFunctions().Add(f_fit)
                        sub_dir.WriteObject(py_hist, py_hist.GetName(), 'kOverwrite')

            if ref not in summary_cscint:
                summary_cscint[ref] = {}
            summary_cscint[ref][fixed_ch] = [
                float(c_scint),
                float(delta_c_fit),
                float(delta_c_syst)
            ]

            s, p, b = self.muAna.parseDetID(detID)
            ch_idx = (p * 10 + b) * 16 + sipm
            qa_cscint[ref]['fit_chi2'].append(chi2 / ndf)
            qa_cscint[ref]['speeds'].append(float(c_scint))
            qa_cscint[ref]['speed_errs'].append(float(np.sqrt(delta_c_fit**2 + delta_c_syst**2)))
            qa_cscint[ref]['ch_indices'].append(ch_idx)

        if root_file is not None:
            sum_root_dir = root_file.GetDirectory('Summary')
            if not sum_root_dir:
                sum_root_dir = root_file.mkdir('Summary')

            for ref in ('DS', 'SF'):
                data = qa_cscint.get(ref, {})
                if not data or not data.get('speeds'):
                    continue
                ref_sum_dir = sum_root_dir.GetDirectory(ref)
                if not ref_sum_dir:
                    ref_sum_dir = sum_root_dir.mkdir(ref)

                if data['slice_chi2']:
                    h_slice_chi2 = ROOT.TH1F(f'h_chi2_ndf_cscint_slices_{ref}', f'Gaussian Slice Fit #chi^{{2}}_{{#nu}} (c_{{scint}} extraction - {ref});#chi^{{2}} / ndf;Slice Counts', 100, 0.0, 10.0)
                    h_slice_chi2.SetDirectory(0)
                    for sc in data['slice_chi2']:
                        h_slice_chi2.Fill(sc)
                    ref_sum_dir.WriteObject(h_slice_chi2, h_slice_chi2.GetName(), 'kOverwrite')
                else:
                    h_slice_chi2 = None

                h_fit_chi2 = ROOT.TH1F(f'h_chi2_ndf_cscint_fits_{ref}', f'Effective Speed Linear Fit #chi^{{2}}_{{#nu}} ({ref});#chi^{{2}} / ndf;Channels', 100, 0.0, 10.0)
                h_fit_chi2.SetDirectory(0)
                for fc in data['fit_chi2']:
                    h_fit_chi2.Fill(fc)
                ref_sum_dir.WriteObject(h_fit_chi2, h_fit_chi2.GetName(), 'kOverwrite')

                h_speed = ROOT.TH1F(f'h_cscint_{ref}', f'Effective Speed c_{{scint}} ({ref});c_{{scint}} [cm/ns];Channels', 100, 10.0, 20.0)
                h_speed.SetDirectory(0)
                for sp in data['speeds']:
                    h_speed.Fill(sp)
                ref_sum_dir.WriteObject(h_speed, h_speed.GetName(), 'kOverwrite')

                n_pts = len(data['ch_indices'])
                g_speed = ROOT.TGraphErrors(n_pts)
                g_speed.SetName(f'g_cscint_vs_channel_{ref}')
                g_speed.SetTitle(f'Effective Speed c_{{scint}} vs Channel ({ref});Channel Index ((plane#times10+bar)#times16 + sipm);c_{{scint}} [cm/ns]')
                g_speed.SetMarkerStyle(20)
                g_speed.SetMarkerSize(0.6)
                g_speed.SetMarkerColor(ROOT.kAzure + 2)
                for i_pt, (ch_i, sp_val, sp_err) in enumerate(zip(data['ch_indices'], data['speeds'], data['speed_errs'])):
                    g_speed.SetPoint(i_pt, ch_i, sp_val)
                    g_speed.SetPointError(i_pt, 0.0, sp_err)
                ref_sum_dir.WriteObject(g_speed, g_speed.GetName(), 'kOverwrite')

                c_cscint_qa = ROOT.TCanvas(f'c_cscint_qa_{ref}', f'c_scint Quality Assurance ({ref})', 1200, 900)
                c_cscint_qa.Divide(2, 2)
                c_cscint_qa.cd(1)
                h_speed.Draw()
                c_cscint_qa.cd(2)
                if h_slice_chi2:
                    h_slice_chi2.Draw()
                c_cscint_qa.cd(3)
                h_fit_chi2.Draw()
                c_cscint_qa.cd(4)
                g_speed.Draw('AP')
                c_cscint_qa.Update()
                ref_sum_dir.WriteObject(c_cscint_qa, c_cscint_qa.GetName(), 'kOverwrite')
                c_cscint_qa.Close()
                del c_cscint_qa

        missing_sf = 0
        for ch, vals in summary_cscint['DS'].items():
            if ch not in summary_cscint['SF']:
                summary_cscint['SF'][ch] = vals
                missing_sf += 1
        if missing_sf > 0:
            print(f"[INFO] Populated {missing_sf} SF channels with DS c_scint physical fallback values.")

        self.cscintvalues_by_ref = summary_cscint
        self.cscintvalues_summary = summary_cscint
        if self.tds0_min is not None and self.tds0_max is not None:
            self.cscintvalues_summary['tds0_window'] = {
                'mode': self.tds0_mode,
                'min': self.tds0_min,
                'max': self.tds0_max
            }
        self.muAna.cscintvalues = summary_cscint.get('DS', {})
        total_ds = len(summary_cscint.get('DS', {}))
        total_sf = len(summary_cscint.get('SF', {}))
        ROOT.gErrorIgnoreLevel = prev_ignore
        print(f'Fitted c_scint: DS={total_ds} channels, SF={total_sf} channels (minEntries: DS={min_entries_ds}, SF={min_entries_sf}).')
        return summary_cscint

    def FitAndSaveTimeWalk(self, root_file):
        """Fit 5-parameter time-walk model with NLL systematic error floor and export 2-pad canvases for both reference systems."""
        ROOT.gStyle.SetPalette(ROOT.kRainBow)
        ROOT.gStyle.SetNumberContours(255)
        qdc_min = getattr(self.options, 'qdcMin', 1.0)
        min_bin_entries = getattr(self.options, 'minBinEntries', 10)
        min_total_entries_ds = getattr(self.options, 'minTotalEntries', 500)
        min_total_entries_sf = getattr(self.options, 'minTotalEntriesSF', min(150, min_total_entries_ds))

        def meets_tw_threshold(hname, h2):
            if not hname.startswith('dtvqdc_') or hname.endswith('Cut'):
                return False
            is_sf = '_SF_' in hname
            th = min_total_entries_sf if is_sf else min_total_entries_ds
            return h2.GetEntries() >= th

        summary_tw = {'DS': {}, 'SF': {}}
        summary_alignment = {'DS': {}, 'SF': {}}
        qa_tw = {
            'DS': {'chi2_stat': [], 'chi2_tot': [], 'sigma_sys': [], 'd_val': [], 'd_unc': [], 'rms_res': [], 'ch_indices': []},
            'SF': {'chi2_stat': [], 'chi2_tot': [], 'sigma_sys': [], 'd_val': [], 'd_unc': [], 'rms_res': [], 'ch_indices': []}
        }
        target_tw_hists = [h for h in sorted(self.hists.items()) if meets_tw_threshold(h[0], h[1])]
        total_tw = len(target_tw_hists)

        for idx, (hname, h2) in enumerate(target_tw_hists):
            if (idx + 1) % 50 == 0 or (idx + 1) == total_tw:
                print(f"  [Progress Time-Walk] Fitted {idx + 1}/{total_tw} channels ({(idx + 1) / total_tw * 100:.0f}%)...", flush=True)

            parts = hname.split('_')
            # format: dtvqdc_<ref>_<detID>_<SiPM>_<state> or dtvqdc_<detID>_<SiPM>_<state>
            if len(parts) >= 5 and parts[1] in ('DS', 'SF'):
                ref = parts[1]
                detID = int(parts[2])
                sipm = int(parts[3])
            elif len(parts) >= 4:
                ref = 'DS'
                detID = int(parts[1])
                sipm = int(parts[2])
            else:
                continue
            fixed_ch = f'{detID}_{sipm}'

            if self.muAna.IsSmallSiPMchannel(sipm) and not getattr(self.options, 'fitSmallSiPMs', False):
                continue
            if fixed_ch == '21005_15':
                continue

            # 1. Timing Projections per QDC bin: Extract Modal Value (MPV) and SEM
            nbins = h2.GetNbinsX()
            x_pts, y_pts, ey_pts = [], [], []

            for b in range(1, nbins + 1):
                x_val = h2.GetXaxis().GetBinCenter(b)
                if x_val < qdc_min:
                    continue
                py = h2.ProjectionY(f'py_tw_{ref}_{fixed_ch}_{b}', b, b)
                entries = py.GetEntries()
                if entries < min_bin_entries:
                    continue

                max_bin = py.GetMaximumBin()
                y_val = py.GetBinCenter(max_bin)

                rms_val = py.GetRMS()
                ey_val = rms_val / np.sqrt(entries) if entries > 0 else 0.0
                if ey_val <= 0:
                    continue

                x_pts.append(x_val)
                y_pts.append(y_val)
                ey_pts.append(ey_val)

            min_pts_required = 6 if is_sf else 10
            if len(x_pts) < min_pts_required:
                continue

            x_arr = np.array(x_pts, dtype=float)
            y_arr = np.array(y_pts, dtype=float)
            ey_arr = np.array(ey_pts, dtype=float)

            qdc_min_fit = float(x_arr[0])
            qdc_max_fit = float(x_arr[-1])

            # 2. Fit 5-parameter rational function (Thesis Table 4.2 bounds)
            high_qdc_mask = x_arr > (0.6 * qdc_max_fit)
            t0_est = float(np.median(y_arr[high_qdc_mask])) if np.sum(high_qdc_mask) > 3 else float(y_arr[-1])

            p0 = [t0_est, -50.0, 5.0, -2.0, 0.005]
            bounds_lower = [0.5 * t0_est if t0_est > 0 else 1.5 * t0_est, -500.0, 0.01, -10.0, 0.0]
            bounds_upper = [1.5 * t0_est if t0_est > 0 else 0.5 * t0_est, 0.0, 30.0, 0.0, 1.0]

            try:
                popt, _ = curve_fit(tw_func, x_arr, y_arr, p0=p0, bounds=(bounds_lower, bounds_upper), maxfev=10000)
            except Exception:
                continue

            res_nll = minimize_scalar(nll_objective, bounds=(0.001, 1.0), method='bounded', args=(x_arr, y_arr, ey_arr, popt))
            sigma_sys = float(res_nll.x) if res_nll.success else 0.350

            t0, alpha, beta, qdc0, gamma = [float(v) for v in popt]
            y_fit = tw_func(x_arr, *popt)
            ey_tot = np.sqrt(ey_arr**2 + sigma_sys**2)
            ndf = max(1, len(x_arr) - 5)

            chi2_stat = float(np.sum(((y_arr - y_fit) ** 2) / (ey_arr**2)))
            chi2_stat_ndf = chi2_stat / ndf

            chi2_tot = float(np.sum(((y_arr - y_fit) ** 2) / (ey_tot**2)))
            chi2_tot_ndf = chi2_tot / ndf
            rms_res = float(np.sqrt(np.mean((y_arr / y_fit - 1.0) ** 2)))

            ratio = y_arr / y_fit
            eratio_stat = ey_arr / np.abs(y_fit)
            eratio_tot = ey_tot / np.abs(y_fit)

            dx_band = np.zeros_like(x_arr)
            if len(x_arr) > 1:
                dx_band.fill(0.5)

            ROOT.TH1.AddDirectory(False)
            cname = f'c_timewalk_fit_{ref}_{fixed_ch}'
            c = ROOT.TCanvas(cname, f'US Time-Walk Fit - Channel {fixed_ch} ({ref} Reference - Run {self.runNr})', 1100, 900)

            pad1 = ROOT.TPad(f'pad1_{ref}_{fixed_ch}', 'Upper Pad', 0.0, 0.28, 1.0, 1.0)
            pad2 = ROOT.TPad(f'pad2_{ref}_{fixed_ch}', 'Lower Pad', 0.0, 0.0, 1.0, 0.28)

            pad1.SetBottomMargin(0.02)
            pad1.SetLeftMargin(0.12)
            pad1.SetRightMargin(0.12)
            pad1.SetGrid()

            pad2.SetTopMargin(0.03)
            pad2.SetBottomMargin(0.32)
            pad2.SetLeftMargin(0.12)
            pad2.SetRightMargin(0.12)
            pad2.SetGrid()

            pad1.Draw()
            pad2.Draw()

            pad1.cd()
            ROOT.gStyle.SetPalette(ROOT.kRainBow)
            h2_clone = h2.Clone(f'h2_dtvqdc_{ref}_{fixed_ch}')
            h2_clone.SetTitle('')
            x_max_plot = max(180.0, qdc_max_fit * 1.05)
            h2_clone.GetXaxis().SetRangeUser(0, x_max_plot)
            h2_clone.GetXaxis().SetLabelSize(0)
            h2_clone.GetYaxis().SetTitle('dt_{SiPM}^{ToF} [ns]')
            h2_clone.GetYaxis().SetTitleSize(0.055)
            h2_clone.GetYaxis().SetTitleOffset(0.95)
            if ref == 'SF':
                h2_clone.GetYaxis().SetRangeUser(-45.0, 5.0)
            else:
                h2_clone.GetYaxis().SetRangeUser(-10.0, 25.0)
            h2_clone.Draw('COLZ')

            # Hatched error band in ROOT.kBlue + 2 (FillStyle 3004) for optimal visibility over jet (kRainBow) colormap
            g_data_tot = ROOT.TGraphErrors(len(x_arr), x_arr, y_arr, dx_band, ey_tot)
            g_data_tot.SetName(f'g_data_total_error_{ref}_{fixed_ch}')
            g_data_tot.SetFillColor(ROOT.kBlue + 2)
            g_data_tot.SetFillStyle(3004)
            g_data_tot.SetLineColor(ROOT.kBlue + 2)
            g_data_tot.Draw('E2 SAME')

            g_data_stat = ROOT.TGraphErrors(len(x_arr), x_arr, y_arr, np.zeros_like(x_arr), ey_arr)
            g_data_stat.SetName(f'g_data_stat_error_{ref}_{fixed_ch}')
            g_data_stat.SetMarkerStyle(20)
            g_data_stat.SetMarkerSize(0.6)
            g_data_stat.SetMarkerColor(ROOT.kBlack)
            g_data_stat.SetLineColor(ROOT.kBlack)
            g_data_stat.Draw('PE SAME')

            fit_tf1 = ROOT.TF1(f'fit_tw_{ref}_{fixed_ch}', '[0] + [1] / ([2]*x - [3]) + [4]*x', qdc_min_fit, qdc_max_fit)
            fit_tf1.SetParameters(t0, alpha, beta, qdc0, gamma)
            fit_tf1.SetLineColor(ROOT.kRed + 1)
            fit_tf1.SetLineWidth(3)
            fit_tf1.Draw('SAME')

            # Legend in top-left quadrant matching reference plot
            leg = ROOT.TLegend(0.18, 0.58, 0.72, 0.92)
            leg.SetBorderSize(1)
            leg.SetFillColorAlpha(ROOT.kWhite, 0.90)
            leg.SetTextFont(42)
            leg.SetTextSize(0.032)
            leg.AddEntry(g_data_stat, 'Data (Statistical SEM)', 'lep')
            leg.AddEntry(g_data_tot, 'Total Uncertainty (SEM #oplus #sigma_{sys})', 'f')
            leg.AddEntry(fit_tf1, f'Fit: #chi^{{2}}_{{stat}}/#nu = {chi2_stat_ndf:.1f} #rightarrow #chi^{{2}}_{{tot}}/#nu = {chi2_tot_ndf:.2f}', 'l')
            leg.AddEntry('', f't_{{0}} = {t0:.2f} ns, #alpha = {alpha:.2f}, #beta = {beta:.2f}', '')
            leg.AddEntry('', f'QDC_{{0}} = {qdc0:.2f}, #gamma = {gamma:.5f}', '')
            leg.AddEntry('', f'#sigma_{{sys}} = {sigma_sys*1000:.1f} ps, RMS res = {rms_res*100:.2f}%', '')
            leg.Draw()

            # Top-left title box matching reference plot
            pt_title = ROOT.TPaveText(0.0, 0.96, 0.45, 1.0, 'NDC')
            pt_title.SetBorderSize(1)
            pt_title.SetFillColor(ROOT.kWhite)
            pt_title.SetTextAlign(12)
            pt_title.SetTextFont(62)
            pt_title.SetTextSize(0.038)
            pt_title.AddText(f'US Channel {fixed_ch} (Run {str(self.runNr).zfill(6)})')
            pt_title.Draw()

            pad2.cd()
            g_ratio_tot = ROOT.TGraphErrors(len(x_arr), x_arr, ratio, dx_band, eratio_tot)
            g_ratio_tot.SetName(f'g_ratio_total_error_{ref}_{fixed_ch}')
            g_ratio_tot.SetTitle(';QDC_{SiPM} [a.u.];Data / Fit')
            g_ratio_tot.SetFillColor(ROOT.kBlue + 2)
            g_ratio_tot.SetFillStyle(3004)
            g_ratio_tot.SetLineColor(ROOT.kBlue + 2)
            g_ratio_tot.GetXaxis().SetRangeUser(0, x_max_plot)
            g_ratio_tot.GetYaxis().SetRangeUser(0.92, 1.08)
            g_ratio_tot.GetXaxis().SetTitleSize(0.12)
            g_ratio_tot.GetXaxis().SetLabelSize(0.09)
            g_ratio_tot.GetXaxis().SetTitleOffset(0.95)
            g_ratio_tot.GetYaxis().SetTitleSize(0.11)
            g_ratio_tot.GetYaxis().SetLabelSize(0.08)
            g_ratio_tot.GetYaxis().SetTitleOffset(0.50)
            g_ratio_tot.GetYaxis().SetNdivisions(505)
            g_ratio_tot.Draw('AE2')

            g_ratio_stat = ROOT.TGraphErrors(len(x_arr), x_arr, ratio, np.zeros_like(x_arr), eratio_stat)
            g_ratio_stat.SetName(f'g_ratio_stat_error_{ref}_{fixed_ch}')
            g_ratio_stat.SetMarkerStyle(20)
            g_ratio_stat.SetMarkerSize(0.6)
            g_ratio_stat.SetMarkerColor(ROOT.kBlack)
            g_ratio_stat.SetLineColor(ROOT.kBlack)
            g_ratio_stat.Draw('PE SAME')

            line1 = ROOT.TLine(0, 1.0, x_max_plot, 1.0)
            line1.SetLineColor(ROOT.kRed + 1)
            line1.SetLineStyle(2)
            line1.SetLineWidth(2)
            line1.Draw('SAME')

            c.Update()

            # 3. Fine alignment parameter extraction (d_SiPM) via truncated mean of TW-corrected residual
            if ref == 'SF':
                h_res = ROOT.TH1F(f'dt_residual_{ref}_{fixed_ch}', f'{fixed_ch} TW-corrected residual ({ref});dt_{{SiPM}}^{{TW, ToF}} [ns];Counts', 500, -45.0, 5.0)
            else:
                h_res = ROOT.TH1F(f'dt_residual_{ref}_{fixed_ch}', f'{fixed_ch} TW-corrected residual ({ref});dt_{{SiPM}}^{{TW, ToF}} [ns];Counts', 350, -10.0, 25.0)
            h_res.SetDirectory(0)

            for bx in range(1, h2.GetNbinsX() + 1):
                qdc_c = h2.GetXaxis().GetBinCenter(bx)
                if qdc_c < qdc_min_fit or qdc_c > qdc_max_fit:
                    continue
                denom_c = beta * qdc_c - qdc0
                if denom_c <= 1e-4:
                    denom_c = 1e-4
                tw_val = (alpha / denom_c) + gamma * qdc_c

                for by in range(1, h2.GetNbinsY() + 1):
                    cnt = h2.GetBinContent(bx, by)
                    if cnt <= 0:
                        continue
                    dt_c = h2.GetYaxis().GetBinCenter(by)
                    dt_res = dt_c - tw_val
                    h_res.Fill(dt_res, cnt)

            if h_res.GetEntries() >= 10:
                max_bin = h_res.GetMaximumBin()
                t_mode = h_res.GetBinCenter(max_bin)
                std_dev = h_res.GetStdDev()
                low_bin = h_res.FindBin(t_mode - std_dev)
                high_bin = h_res.FindBin(t_mode + std_dev)
                sum_w = 0.0
                sum_wx = 0.0
                sum_wxx = 0.0
                for b in range(low_bin, high_bin + 1):
                    w = h_res.GetBinContent(b)
                    if w <= 0:
                        continue
                    x = h_res.GetBinCenter(b)
                    sum_w += w
                    sum_wx += w * x
                    sum_wxx += w * x * x
                if sum_w > 0:
                    d_val = float(sum_wx / sum_w)
                    var = (sum_wxx / sum_w) - d_val**2
                    sigma_trunc = float(np.sqrt(max(0.0, var)))
                    d_unc = float(sigma_trunc / np.sqrt(sum_w))
                else:
                    d_val = float(t0)
                    d_unc = float(sigma_sys)
            else:
                d_val = float(t0)
                d_unc = float(sigma_sys)

            if ref not in summary_alignment:
                summary_alignment[ref] = {}
            summary_alignment[ref][fixed_ch] = [d_val, d_unc]

            # Save into channel's TDirectory
            if root_file is not None:
                ch_dir = root_file.GetDirectory(fixed_ch)
                if not ch_dir:
                    ch_dir = root_file.mkdir(fixed_ch)
                ch_dir.WriteObject(c, cname, 'kOverwrite')
                ch_dir.WriteObject(fit_tf1, fit_tf1.GetName(), 'kOverwrite')
                ch_dir.WriteObject(g_data_stat, g_data_stat.GetName(), 'kOverwrite')
                ch_dir.WriteObject(g_data_tot, g_data_tot.GetName(), 'kOverwrite')
                ch_dir.WriteObject(g_ratio_stat, g_ratio_stat.GetName(), 'kOverwrite')
                ch_dir.WriteObject(g_ratio_tot, g_ratio_tot.GetName(), 'kOverwrite')
                ch_dir.WriteObject(h_res, h_res.GetName(), 'kOverwrite')
                c.Close()
                del c
            del h_res

            if ref not in summary_tw:
                summary_tw[ref] = {}
            summary_tw[ref][fixed_ch] = {
                't0': t0,
                'alpha': alpha,
                'beta': beta,
                'qdc0': qdc0,
                'gamma': gamma,
                'sigma_sys': sigma_sys,
                'chi2_stat_ndf': chi2_stat_ndf,
                'chi2_tot_ndf': chi2_tot_ndf,
                'rms_residual': rms_res,
                'qdc_limits': [qdc_min_fit, qdc_max_fit]
            }

            s, p, b = self.muAna.parseDetID(detID)
            ch_idx = (p * 10 + b) * 16 + sipm
            qa_tw[ref]['chi2_stat'].append(chi2_stat_ndf)
            qa_tw[ref]['chi2_tot'].append(chi2_tot_ndf)
            qa_tw[ref]['sigma_sys'].append(sigma_sys * 1000.0)
            qa_tw[ref]['d_val'].append(d_val)
            qa_tw[ref]['d_unc'].append(d_unc)
            qa_tw[ref]['rms_res'].append(rms_res * 100.0)
            qa_tw[ref]['ch_indices'].append(ch_idx)

        if root_file is not None:
            sum_root_dir = root_file.GetDirectory('Summary')
            if not sum_root_dir:
                sum_root_dir = root_file.mkdir('Summary')

            for ref in ('DS', 'SF'):
                data = qa_tw.get(ref, {})
                if not data or not data.get('ch_indices'):
                    continue
                ref_sum_dir = sum_root_dir.GetDirectory(ref)
                if not ref_sum_dir:
                    ref_sum_dir = sum_root_dir.mkdir(ref)

                # 1. chi2/ndf timewalk fits (stat and tot)
                h_tw_chi2_stat = ROOT.TH1F(f'h_chi2_ndf_tw_stat_{ref}', f'Time-Walk Fit #chi^{{2}}_{{stat}} / #nu ({ref});#chi^{{2}}_{{stat}} / #nu;Channels', 100, 0.0, 50.0)
                h_tw_chi2_stat.SetDirectory(0)
                for c_stat in data['chi2_stat']:
                    h_tw_chi2_stat.Fill(c_stat)
                ref_sum_dir.WriteObject(h_tw_chi2_stat, h_tw_chi2_stat.GetName(), 'kOverwrite')

                h_tw_chi2_tot = ROOT.TH1F(f'h_chi2_ndf_tw_tot_{ref}', f'Time-Walk Fit #chi^{{2}}_{{tot}} / #nu with #sigma_{{sys}} ({ref});#chi^{{2}}_{{tot}} / #nu;Channels', 100, 0.0, 5.0)
                h_tw_chi2_tot.SetDirectory(0)
                for c_tot in data['chi2_tot']:
                    h_tw_chi2_tot.Fill(c_tot)
                ref_sum_dir.WriteObject(h_tw_chi2_tot, h_tw_chi2_tot.GetName(), 'kOverwrite')

                # 2. Systematic error floor distribution
                h_sigma_sys = ROOT.TH1F(f'h_sigma_sys_{ref}', f'Time-Walk Systematic Error Floor #sigma_{{sys}} ({ref});#sigma_{{sys}} [ps];Channels', 100, 0.0, 500.0)
                h_sigma_sys.SetDirectory(0)
                for s_sys in data['sigma_sys']:
                    h_sigma_sys.Fill(s_sys)
                ref_sum_dir.WriteObject(h_sigma_sys, h_sigma_sys.GetName(), 'kOverwrite')

                d_vals = data['d_val']
                d_min = min(d_vals) - 1.0 if d_vals else -5.0
                d_max = max(d_vals) + 1.0 if d_vals else 20.0
                h_alignment = ROOT.TH1F(f'h_alignment_d_{ref}', f'Fine Alignment Parameter d_{{SiPM}} ({ref});d_{{SiPM}} [ns];Channels', 100, d_min, d_max)
                h_alignment.SetDirectory(0)
                for d_v in d_vals:
                    h_alignment.Fill(d_v)
                ref_sum_dir.WriteObject(h_alignment, h_alignment.GetName(), 'kOverwrite')

                h_rms = ROOT.TH1F(f'h_rms_residual_{ref}', f'Fractional RMS Residuals ({ref});RMS Residual [%];Channels', 100, 0.0, 5.0)
                h_rms.SetDirectory(0)
                for r_v in data['rms_res']:
                    h_rms.Fill(r_v)
                ref_sum_dir.WriteObject(h_rms, h_rms.GetName(), 'kOverwrite')

                n_tw_pts = len(data['ch_indices'])
                g_sigma_ch = ROOT.TGraph(n_tw_pts)
                g_sigma_ch.SetName(f'g_sigma_sys_vs_channel_{ref}')
                g_sigma_ch.SetTitle(f'Time-Walk #sigma_{{sys}} vs Channel ({ref});Channel Index ((plane#times10+bar)#times16 + sipm);#sigma_{{sys}} [ps]')
                g_sigma_ch.SetMarkerStyle(20)
                g_sigma_ch.SetMarkerSize(0.6)
                g_sigma_ch.SetMarkerColor(ROOT.kAzure + 2)
                for i_pt, (ch_i, s_sys) in enumerate(zip(data['ch_indices'], data['sigma_sys'])):
                    g_sigma_ch.SetPoint(i_pt, ch_i, s_sys)
                ref_sum_dir.WriteObject(g_sigma_ch, g_sigma_ch.GetName(), 'kOverwrite')

                g_align_ch = ROOT.TGraphErrors(n_tw_pts)
                g_align_ch.SetName(f'g_alignment_vs_channel_{ref}')
                g_align_ch.SetTitle(f'Fine Alignment d_{{SiPM}} vs Channel ({ref});Channel Index ((plane#times10+bar)#times16 + sipm);d_{{SiPM}} [ns]')
                g_align_ch.SetMarkerStyle(20)
                g_align_ch.SetMarkerSize(0.6)
                g_align_ch.SetMarkerColor(ROOT.kTeal + 2)
                for i_pt, (ch_i, d_v, d_u) in enumerate(zip(data['ch_indices'], data['d_val'], data['d_unc'])):
                    g_align_ch.SetPoint(i_pt, ch_i, d_v)
                    g_align_ch.SetPointError(i_pt, 0.0, d_u)
                ref_sum_dir.WriteObject(g_align_ch, g_align_ch.GetName(), 'kOverwrite')

                c_tw_qa = ROOT.TCanvas(f'c_tw_qa_{ref}', f'Time-Walk Quality Assurance ({ref})', 1200, 900)
                c_tw_qa.Divide(2, 2)
                c_tw_qa.cd(1)
                h_tw_chi2_tot.Draw()
                c_tw_qa.cd(2)
                h_sigma_sys.Draw()
                c_tw_qa.cd(3)
                h_rms.Draw()
                c_tw_qa.cd(4)
                g_sigma_ch.Draw('AP')
                c_tw_qa.Update()
                ref_sum_dir.WriteObject(c_tw_qa, c_tw_qa.GetName(), 'kOverwrite')
                c_tw_qa.Close()
                del c_tw_qa

                c_align_qa = ROOT.TCanvas(f'c_alignment_qa_{ref}', f'Fine Alignment Quality Assurance ({ref})', 1200, 600)
                c_align_qa.Divide(2, 1)
                c_align_qa.cd(1)
                h_alignment.Draw()
                c_align_qa.cd(2)
                g_align_ch.Draw('AP')
                c_align_qa.Update()
                ref_sum_dir.WriteObject(c_align_qa, c_align_qa.GetName(), 'kOverwrite')
                c_align_qa.Close()
                del c_align_qa

        out_dir = self.options.outputDir
        os.makedirs(out_dir, exist_ok=True)
        # Fallback for SF channels lacking sufficient statistics using DS timing reference
        # Δt0 = t0_SF - t0_DS offset calculated from successfully fitted SF channels
        delta_t0_vals = [
            summary_tw['SF'][ch]['t0'] - summary_tw['DS'][ch]['t0']
            for ch in summary_tw.get('SF', {})
            if ch in summary_tw.get('DS', {})
        ]
        delta_t0_median = float(np.median(delta_t0_vals)) if delta_t0_vals else -12.35

        missing_sf_tw = 0
        for ch, tw_ds in summary_tw.get('DS', {}).items():
            if ch not in summary_tw.get('SF', {}):
                tw_fallback = dict(tw_ds)
                tw_fallback['t0'] = float(tw_ds['t0'] + delta_t0_median)
                tw_fallback['is_fallback'] = True
                summary_tw.setdefault('SF', {})[ch] = tw_fallback
                missing_sf_tw += 1

        missing_sf_align = 0
        for ch, al_ds in summary_alignment.get('DS', {}).items():
            if ch not in summary_alignment.get('SF', {}):
                al_fallback = [float(al_ds[0] + delta_t0_median), float(al_ds[1])]
                summary_alignment.setdefault('SF', {})[ch] = al_fallback
                missing_sf_align += 1

        if missing_sf_tw > 0:
            print(f"[INFO] Populated {missing_sf_tw} SF channels with DS time-walk fallback parameters (median offset: {delta_t0_median:.2f} ns).", flush=True)

        self.alignmentparameters_by_ref = summary_alignment
        consolidated_json = os.path.join(out_dir, f'run{self.runNr}_timing_calibration.json')
        combined_calib = {
            'DS': {
                'cscint': getattr(self, 'cscintvalues_by_ref', {}).get('DS', {}),
                'timewalk': summary_tw.get('DS', {}),
                'alignment': summary_alignment.get('DS', {})
            },
            'SF': {
                'cscint': getattr(self, 'cscintvalues_by_ref', {}).get('SF', {}),
                'timewalk': summary_tw.get('SF', {}),
                'alignment': summary_alignment.get('SF', {})
            }
        }
        with open(consolidated_json, 'w') as jf:
            json.dump(combined_calib, jf, indent=2)

        total_ds_tw = len(summary_tw.get('DS', {}))
        total_sf_tw = len(summary_tw.get('SF', {}))
        total_ds_align = len(summary_alignment.get('DS', {}))
        total_sf_align = len(summary_alignment.get('SF', {}))
        print(f'Fitted Time-Walk: DS={total_ds_tw} channels, SF={total_sf_tw} channels.')
        print(f'Fitted Alignment: DS={total_ds_align} channels, SF={total_sf_align} channels.')
        print(f'Consolidated timing calibration saved to {consolidated_json}')

    def LoadSystemObservables(self):
        systemobservablesfilename=f'{self.outputDir}/SystemObservables/run{self.TWCorrectionRun}/systemparameters.json'
        with open(systemobservablesfilename, 'r') as f:
            self.systemobservables=json.load(f)

    def yresidual3(self, detID):

        # self.MuFilter.GetPosition(detID,self.A,self.B)

        doca=self.muAna.Getyresidual(detID)
        s,p,b=self.muAna.parseDetID(detID)

        if s==3: return True
        key=10*s+p
        if key==12: return True # No hist for veto 3 yet

        dy_min, dy_max = self.dycut(self.cutdists[f'dy_{key}'])
        if doca>dy_max or doca<dy_min: return False
        else: return True

    def CheckForSFtrack(self):
        for idx, track in enumerate(self.M.Reco_MuonTracks):
            if track.GetUniqueID()==1:
                return idx
        return -1

    def nSiPMscut(self, hit, nLeft, nRight):
        s,p,b=self.muAna.parseDetID(hit.GetDetectorID())
        if s==1:
            if nLeft<6 or nRight<6: return False
        elif s==2:
            if nLeft<4 or nRight<4: return False
        elif s==3:
            pass
        return True

    def TimingDiscriminantCut(self):
        if self.timealignment=='old':
            if self.td<0: return False
            else: return True
        else:
            if not 'timingdiscriminant' in self.cutdists:
                print(f'No timing discriminant histogram')
                return
            hist=self.cutdists['timingdiscriminant']
            mean, stddev=hist.GetMean(), hist.GetStdDev()
            if self.td < mean-2*stddev or self.td > mean+2*stddev: return False
            else: return True

    def dycut(self, hist, nsig=1):
        dymin=hist.GetMean()-nsig*hist.GetStdDev()
        dymax=hist.GetMean()+nsig*hist.GetStdDev()
        return dymin, dymax

    def slopecut(self, slopecut=0.1):
        slopeX, slopeY = self.mom.x()/self.mom.z(), self.mom.y()/self.mom.z()
        if abs(slopeX)>slopecut or abs(slopeY)>slopecut: return 0
        else: return 1

    def tds0cut(self, val=None):
        if getattr(self.options, 'noTds0Cut', False):
            return 1
        t = self.reft if val is None else val
        if t is None or t == -999.:
            return 0
        if self.tds0_min is not None and self.tds0_max is not None:
            return 1 if self.tds0_min <= t <= self.tds0_max else 0
        time_align = getattr(self, 'timealignment', 'new')
        if time_align == 'old':
            default_min, default_max = (20.25, 24.15)
        elif time_align == 'new':
            default_min, default_max = (13.0, 15.5)
        else:
            default_min, default_max = (18.5, 20.5)
        return 1 if default_min <= t <= default_max else 0
