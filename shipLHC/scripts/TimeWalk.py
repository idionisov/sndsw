import ROOT, os, csv, json
from AnalysisFunctions import Analysis
from pathlib import Path
import numpy as np

deciles=[i/10 for i in range(11)]
SiPM_distance_cut=5

ROOT.gStyle.SetTitleXOffset(0.7)
ROOT.gStyle.SetTitleYOffset(0.6)
for i in ('x', 'y', 'z','t'): ROOT.gStyle.SetTitleSize(0.04, i)

class TimeWalk(ROOT.FairTask):

    def __init__(self, options, monitor):

        self.M=monitor
        self.simulation = options.simulation
        self.muAna = Analysis(options)
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

        if options.mode.find('1')>0:
            self.referencesystem=1
            self.mode=options.mode[:-1]
        elif options.mode.find('3')>0:
            self.referencesystem=3
            self.mode=options.mode[:-1]
        else:
            self.referencesystem=options.referencesystem
            self.mode=options.mode
        self.refsysname='DS' if self.referencesystem==3 else 'SF'
        self.path='TI18'
        if options.path.find('testbeam_June2023_H8')>0:
            self.path='H8'

        #### Time walk correction run is independent of time alignment settings.
        self.TWCorrectionRun=str(5408).zfill(6)

        if self.simulation:
            self.timealignment='sim'
            self.GetLumis() # loads the simulation luminosities dictionary into tw.simlumis
            self.simMode = options.simMode
            # self.cutdists=self.muAna.GetCutDistributions("005408", ('dy', 'timingdiscriminant'))
            self.cutdists=self.muAna.GetCutDistributions("005408", ('dy', ))

        elif not self.simulation:
            statedict={'zeroth':'uncorrected', 'tof':'uncorrected',
                        'tw':'corrected', 'res':'corrected',
                        'selectioncriteria':'corrected',
                        'systemalignment':'corrected',
                        'reconstructmuonposition':'corrected',
                        'showerprofiles':'corrected',
                        'numusignalevents':'corrected',
                        'tds0-studies':'uncorrected',
                        'extendedreconstruction':'corrected', 'struckquark':'corrected'}

            self.state=statedict[self.mode]

            self.timealignment=self.muAna.GetTimeAlignmentType(runNr=self.runNr)

            if options.debug: self.trackevents=[]

            ### If no time-walk correction run is provided. Set the default correction run depending on time alignment of the data set
            if options.AlignmentRun==None or self.muAna.GetTimeAlignmentType(runNr=options.AlignmentRun) != self.timealignment:
                if self.timealignment=='old': self.AlignmentRun=str(5097).zfill(6)
                elif self.timealignment=='new': self.AlignmentRun=str(5408).zfill(6)
                elif self.timealignment=='new+LHCsynch': self.AlignmentRun=str(5999).zfill(6)
            self.muAna.AlignmentRun=self.AlignmentRun

            if self.mode not in ['tds0-studies', 'zeroth']:
                self.cutdists = self.muAna.GetCutDistributions(self.TWCorrectionRun, ('dy', ))

                self.muAna.MakeAlignmentParameterDict()
                self.muAna.Makecscintdict(self.TWCorrectionRun, state=self.state)

                n=9
                self.muAna.MakeTWCorrectionDict(n=n)
            else:
                self.cutdists = {}

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
                    "corrected": [(100, -80, 20), (200, -25, 25)]
                },
                "dtvqdc": {
                    "uncorrected": [(200, 0, 200), (200, 0, 50)],
                    "corrected": [(200, 0, 200), (200, -25, 25)]
                }
            }

        self.notInDS=0

        self.trackRequired=True if self.mode in('zeroth', 'tof', 'tw', 'res', 'systemalignment', 'reconstructmuonposition', 'selectioncriteria', 'showerprofiles', 'tds0-studies') else False

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

        # Only accept events associated with colliding IP1 bunches
        if not self.simulation:
            # M.xing does not exist for simulation runs
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
            if track.GetUniqueID()==1:
                inScifi=True
                tracks[1].append(self.M.Reco_MuonTracks[i])
            if track.GetUniqueID()==3:
                inDS=True
                tracks[3].append(self.M.Reco_MuonTracks[i])

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

            # Hist just used for discussion in the thesis
            if self.mode=='zeroth': self.MakeExplantoryHist(event.Digi_MuFilterHits)

            for channel in channels_t:
                SiPM,clock=channel
                if self.muAna.IsSmallSiPMchannel(SiPM): continue
                if f'{detID}_{SiPM}' == '21005_15': continue # Dead SiPM
                qdc=self.muAna.GetChannelVal(SiPM, channels_qdc)
                if qdc==-999.: continue
                fixed_ch=f'{detID}_{SiPM}'
                if self.mode=='zeroth': self.zeroth(fixed_ch, clock, qdc)
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

    def zeroth(self,fixed_ch,clock,qdc):

        ReadableFixedCh=self.muAna.MakeHumanReadableFixedCh(fixed_ch)

        detID, SiPM = ( int(fixed_ch.split('_')[i]) for i in range(2) )
        s,p,b=self.muAna.parseDetID(detID)

        correctedtime=clock*self.TDC2ns
        t_rel=self.reft-correctedtime

        if not self.muAna.DSVcheck(detID):  dtvpred=f'dtvxpred_{fixed_ch}_{self.state}'
        else: dtvpred=f'dtvypred_{fixed_ch}_{self.state}'

        if not dtvpred in self.hists:
            coord='x' if not self.muAna.DSVcheck(detID) else 'y'
            title='{No time-walk correction t_{0}^{DS}-t^{uncorr}_{SiPM} v '+coord+'-position}'
            splittitle='#splitline{'+ReadableFixedCh+'}'+title
            axestitles=coord+'_{predicted} [cm];t_{0}^{DS}-t^{uncorr}_{SiPM} [ns]'
            fulltitle=splittitle+';'+axestitles
            histformat = self.histformatting["dtvxpred"][self.state]
            self.hists[dtvpred]=ROOT.TH2F(dtvpred,fulltitle,*histformat[0], *histformat[1])

        attlen=f'attlen_{fixed_ch}_{self.state}'
        if not attlen in self.hists:
            coord='x' if not self.muAna.DSVcheck(detID) else 'y'
            title='{Predicted position against QDC_{SiPM}}'
            splittitle='#splitline{'+ReadableFixedCh+'}'+title
            axestitles=coord+'_{predicted} [cm];QDC_{SiPM} [a.u]'
            fulltitle=splittitle+';'+axestitles
            self.hists[attlen]=ROOT.TH2F(attlen,fulltitle, 110, -100, 10, 200, 0., 200)

        self.hists[dtvpred].Fill(self.Ex[p].x(), t_rel)
        self.hists[attlen].Fill(self.Ex[p].x(), qdc)

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

    def tof(self, fixed_ch, clock, qdc):

        hists=self.hists
        ReadableFixedCh=self.muAna.MakeHumanReadableFixedCh(fixed_ch)

        detID, SiPM=fixed_ch.split('_')
        s,p,b = self.muAna.parseDetID(int(detID))

        ToFcorrectedtime = self.muAna.MuFilterCorrectedTime(self.MuFilter, fixed_ch, clock, -1, x=self.Ex[p].x())
        t_rel = self.reft - ToFcorrectedtime

        dtvqdc = f'dtvqdc_{fixed_ch}_{self.state}'

        for cut in (None, f'-{SiPM_distance_cut}cmFiducialCut'):

            if cut!=None: histname = dtvqdc + cut
            else: histname = dtvqdc

            if not histname in hists:
                if cut==f'-{SiPM_distance_cut}cmFiducialCut':
                    subtitle='{No tw correction t_{0}^{DS}-t^{uncorr}_{SiPM} v QDC_{SiPM}};QDC_{SiPM} [a.u];t_{0}^{DS}-t^{uncorr}_{SiPM} [ns]'
                else: subtitle='{'+str(SiPM_distance_cut)+' cm fiducial cut, no tw correction t_{0}^{DS}-t^{uncorr}_{SiPM} v QDC_{SiPM}};QDC_{SiPM} [a.u];t_{0}^{DS}-t^{uncorr}_{SiPM} [ns]'

                title='#splitline{'+ReadableFixedCh+'}'+subtitle
                histformat=self.histformatting["dtvqdc"][self.state]
                hists[histname]=ROOT.TH2F(histname, title, *histformat[0],*histformat[1])

        self.hists[dtvqdc].Fill(qdc,t_rel)
        if self.passfiducialvolumecut:
            self.hists[dtvqdc+f'-{SiPM_distance_cut}cmFiducialCut'].Fill(qdc,t_rel)

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
            histformat = self.histformatting['dtvxpred'][self.state]
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
            histformat = self.histformatting['dtvqdc'][self.state]
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

        if self.mode in ('zeroth', 'tof', 'tw', 'res'):
            if getattr(self.options, 'outFile', None):
                outfilename = self.options.outFile
            else:
                outfilename = f'{self.outpath}run{self.runNr}/timewalk_{self.mode}_{self.options.nStart}.root'

            os.makedirs(os.path.dirname(os.path.abspath(outfilename)), exist_ok=True)
            f = ROOT.TFile.Open(outfilename, 'recreate')

            # Global run/reference histograms written at the root level of the file
            global_written = set()
            if self.mode == 'zeroth':
                global_hist_keys = [
                    'reft-SF' if self.referencesystem == 1 else 'reft-DS',
                    'NDSHSiPMs',
                    'reducedchi2'
                ]
                if self.refsysname == 'DS':
                    global_hist_keys.append('tds0vxpred')

                for gkey in global_hist_keys:
                    if gkey in self.hists:
                        f.WriteObject(self.hists[gkey], gkey, 'kOverwrite')
                        global_written.add(gkey)

            # Explanatory histograms written in 24004_4 subdirectory
            if self.mode == 'zeroth':
                for bar_third in (1, 2, 3):
                    histname = f'qdcthird{bar_third}_24004_4'
                    if histname in self.hists:
                        d_exp = f.GetDirectory('24004_4')
                        if not d_exp:
                            d_exp = f.mkdir('24004_4')
                        d_exp.WriteObject(self.hists[histname], histname, 'kOverwrite')
                        global_written.add(histname)

            # Write channel-specific histograms into their dedicated TDirectory (<detID>_<SiPM>)
            for h in self.hists:
                if h in global_written:
                    continue

                parts = h.split('_')
                if len(parts) == 4:
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

            if self.mode == 'zeroth' and not getattr(self.options, 'noCscintFit', False):
                self.FitAndSaveCscint(f)

            f.Close()
            print(f'{len(self.M.h)} histograms saved to {outfilename}')

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

    def FitAndSaveCscint(self, root_file):
        min_entries = getattr(self.options, 'minEntries', 10000)
        fit_min = getattr(self.options, 'fitMin', -65.0)
        fit_max = getattr(self.options, 'fitMax', 5.0)
        min_slice_entries = getattr(self.options, 'minSliceEntries', 15)
        n_segments = getattr(self.options, 'nSegments', 5)

        summary_cscint = {}

        for hname, h2 in sorted(self.hists.items()):
            if not hname.startswith('dtvxpred_'):
                continue
            if h2.GetEntries() < min_entries:
                continue

            parts = hname.split('_')
            # format: dtvxpred_<detID>_<SiPM>_<state>
            if len(parts) >= 4:
                try:
                    detID = int(parts[1])
                    sipm = int(parts[2])
                except ValueError:
                    continue
                fixed_ch = f'{detID}_{sipm}'
            else:
                continue

            if self.muAna.IsSmallSiPMchannel(sipm) and not getattr(self.options, 'fitSmallSiPMs', False):
                continue
            if fixed_ch == '21005_15':
                continue

            gr = ROOT.TGraphErrors()
            gr.SetName(f'gr_mpv_{fixed_ch}')
            gr.SetTitle(f'MPV vs x_{{pred}} {fixed_ch};x_{{pred}} [cm];t_{{0}}^{{DS}} - t_{{SiPM}} [ns]')

            ip = 0
            for bx in range(1, h2.GetNbinsX() + 1):
                x = h2.GetXaxis().GetBinCenter(bx)
                if x < fit_min or x > fit_max:
                    continue
                py = h2.ProjectionY(f'py_{fixed_ch}_{bx}', bx, bx)
                if py.GetEntries() < min_slice_entries:
                    continue

                max_bin = py.GetMaximumBin()
                mode = py.GetBinCenter(max_bin)

                fl = ROOT.TF1(f'fl_{fixed_ch}_{bx}', 'landau', mode - 1.5, mode + 3.0)
                fit_res = py.Fit(fl, 'QNS0')
                mpv = None
                mpv_err = None
                if int(fit_res) == 0 and 0 < fl.GetParError(1) < 2.5:
                    p1 = fl.GetParameter(1)
                    if mode - 3.0 < p1 < mode + 4.0:
                        mpv = p1
                        mpv_err = fl.GetParError(1)

                if mpv is None:
                    fg = ROOT.TF1(f'fg_{fixed_ch}_{bx}', 'gaus', mode - 1.5, mode + 1.5)
                    res_g = py.Fit(fg, 'QNS0')
                    if int(res_g) == 0 and 0 < fg.GetParError(1) < 2.5:
                        mpv = fg.GetParameter(1)
                        mpv_err = fg.GetParError(1)
                    else:
                        mpv = mode
                        mpv_err = py.GetBinWidth(max_bin) / (12.0 ** 0.5)

                gr.SetPoint(ip, x, mpv)
                gr.SetPointError(ip, 0.0, mpv_err)
                ip += 1

            if gr.GetN() < 5:
                continue

            func = ROOT.TF1(f'fit_cscint_{fixed_ch}', 'pol1', fit_min, fit_max)
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

            seg_bounds = np.linspace(fit_min, fit_max, n_segments + 1)
            subrange_speeds = []
            for s_idx in range(n_segments):
                seg_func = ROOT.TF1(f'seg_{fixed_ch}_{s_idx}', 'pol1', seg_bounds[s_idx], seg_bounds[s_idx + 1])
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

            cname = f'c_cscint_{fixed_ch}'
            canv = ROOT.TCanvas(cname, f'Effective Speed Fit {fixed_ch}', 900, 700)
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
            pave.AddText(f'Channel: {fixed_ch}')
            pave.AddText(f'c_{{scint}} = {c_scint:.2f} #pm {delta_c_fit:.2f} (fit) #pm {delta_c_syst:.2f} (syst) cm/ns')
            pave.AddText(f'#chi^{{2}}/ndf = {chi2:.1f}/{ndf} ({chi2/ndf:.2f})')
            pave.Draw('SAME')
            canv.Update()

            ch_dir = root_file.GetDirectory(fixed_ch)
            if not ch_dir:
                ch_dir = root_file.mkdir(fixed_ch)
            ch_dir.WriteObject(canv, cname, 'kOverwrite')
            ch_dir.WriteObject(gr, gr.GetName(), 'kOverwrite')
            ch_dir.WriteObject(func, func.GetName(), 'kOverwrite')

            summary_cscint[fixed_ch] = [
                float(c_scint),
                float(delta_c_fit),
                float(delta_c_syst)
            ]

        out_dir = self.options.outputDir
        os.makedirs(out_dir, exist_ok=True)
        json_filename = os.path.join(out_dir, f'run{self.runNr}_cscintvalues.json')
        with open(json_filename, 'w') as jf:
            json.dump(summary_cscint, jf, indent=2)
        print(f'Fitted c_scint for {len(summary_cscint)} channels (minEntries={min_entries}). Calibration saved to {json_filename}')

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

    def tds0cut(self, nsig=1):
        if self.referencesystem==1: return 1

        if self.reft < 18 or self.reft > 21: return 0
        else: return 1
