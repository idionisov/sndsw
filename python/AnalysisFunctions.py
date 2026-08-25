#!/usr/bin/env python
import ROOT,os,csv,json
from datetime import datetime
import math as m
import numpy as np
from itertools import combinations

class Analysis(object):

	def __init__(self, options):

		self.options=options
		if options.numuStudy or options.nueStudy: self.nuStudy=True
		else: self.nuStudy=False

		# Adding flag for LaserMeasurements/ work to use some functions defined in here
		if not hasattr(options, "LaserMeasurements"):
			if options.runNumber==-1: self.runNr='005408'
			else: self.runNr = str(options.runNumber).zfill(6)
			self.TWCorrectionRun = str(5408).zfill(6)

			self.timealignment=self.GetTimeAlignmentType(self.runNr)
			self.state=options.state
			if hasattr(options, 'datafiletype'): self.fileext=options.datafiletype
			else: self.fileext='json'
			self.CorrectionType=options.CorrectionType

			afswork = getattr(options, 'afswork', '/afs/cern.ch/work/i/idioniso/sndVetoUS')
			afsuser = getattr(options, 'afsuser', '/afs/cern.ch/work/i/idioniso/sndVetoUS/twfiles')

			# if self.options.path.find('commissioning/TI18')>0:
			if options.datalocation=='commissioning':
				self.path=afswork+'-commissioning/'
			elif options.datalocation=='physics':
				self.path=afswork+'-physics2022/'
			elif options.datalocation=='H8':
				self.path=afswork+'-H8/'

			self.referencesystem=options.referencesystem
			self.refsysname='DS' if self.referencesystem==3 else 'SF'

			if options.numuStudy: self.Get_numuevents()
			elif options.nueStudy: self.Get_nueevents()

			self.simulation = options.simulation
			if self.simulation:
				self.simEngine = self.GetSimEngine

		self.correctionparams=lambda ps : [y for x,y in enumerate(ps) if x%2==0]

		# self.correctionfunction = lambda ps, qdc : ps[3]*(qdc-ps[0])/( ps[1] + ps[2]*(qdc-ps[0])*(qdc-ps[0]) ) + ps[4]*(qdc-ps[0])

		# ps = t_0, alpha, beta, QDC_0, gamma
		self.correctionfunction = lambda ps, qdc : ps[1] / (ps[2]*qdc - ps[3]) + ps[4]*qdc

		self.A, self.B = ROOT.TVector3(), ROOT.TVector3()
		self.systemAndPlanes = {1:2,2:5,3:7}
		self.systemAndBars = {1:7,2:10,3:60}
		self.systemAndChannels = {1:[0,8],2:[2,6],3:[0,1]}
		self.systemAndSiPMs={1:range(16),2:(0,1,3,4,6,7,8,9,11,12,14,15),3:(1,)}
		self.verticalBarDict={0:1, 1:3, 2:5, 3:6}
		self.gelsides={0:'right', 1:'left', 2:'right', 3:'left', 4:'left'}
		self.subsystemNames={1:'veto', 2:'upstream', 3:'downstream'}
		self.verbose=False

		if hasattr(options, 'datafiletype'): self.fileext=options.datafiletype
		else: self.fileext='csv'

		self.sigmatds0=0.263, 9.5E-5
		freq = 160.316E6
		self.TDC2ns = 1E9/freq

	def GetSimEngine(self):
		simEngine = self.options.geoFile.split('.')[1].split('-')[0]
		return simEngine

	def print_timestamp(self, message=""):
		print(f"{datetime.now().strftime('%Y-%m-%d %H:%M:%S')} - {message}")

	def SetTask(self, task):
		self.task=task

	def DSHcheck(self, detID): # True if detID is a horizontal bar
		s,p,b=self.parseDetID(detID)
		if s!=3: return False
		if p<3 and b<60:  return True
		else: return False

	def DSVcheck(self, detID): # True if detID is a vertical
		s,p,b=self.parseDetID(detID)
		if s!=3: return False
		if p<3 and b>59: return True
		elif p==3: return True
		else: return False

	def GetListOfChannels(self, subsystem): #Only returns horizontal channels for the DS (s==3)
		channels=[f'{self.MakeFixedCh((subsystem, plane, bar, SiPM))}' for plane in range(self.systemAndPlanes[subsystem]) for bar in range(self.systemAndBars[subsystem]) for SiPM in self.systemAndSiPMs[subsystem] ]
		return channels

	def GetGeoFile(self,runNumber):
		if type(runNumber)==str: runNumber=int(runNumber)
		year = self.GetRunYear(runNumber)
		if year==2022: geoFile =  "geofile_sndlhc_TI18_V4_2022.root"
		elif year==2023: geoFile =  "geofile_sndlhc_TI18_V3_2023.root"
		elif year==2024: geoFile =  "geofile_sndlhc_TI18_V0_2024.root"
		elif year==2025: geoFile =  "geofile_sndlhc_TI18_V0_2025.root"

		return geoFile

	def GetRunYear(self, runNr):
		if isinstance(runNr, str): runNr=int(runNr)

		if runNr < 5485: year='2022'
		elif 5485 <= runNr < 7656 : year='2023'
		else: year='2024'

		return year

	def BuildBarLengths(self, MuFilter):
		Vetobarlength = MuFilter.GetConfParF('MuFilter/VetoBarX')
		USbarlength = MuFilter.GetConfParF('MuFilter/UpstreamBarX')
		DSbarlength_hor = MuFilter.GetConfParF('MuFilter/DownstreamBarX')
		DSbarlength_vert = MuFilter.GetConfParF('MuFilter/UpstreamBarY_ver')

		self.barlengths={1:Vetobarlength, 2:USbarlength, 3:DSbarlength_hor, 4:DSbarlength_vert}

	def IsSmallSiPMchannel(self, i):
		if i==2 or i==5 or i==10 or i==13: return True
		else: return False

	def GetSide(self, fixed_ch):
		detID, SiPM = fixed_ch.split('_')
		s,p,b=self.parseDetID(int(detID))
		if s<3:
			if int(SiPM)<8: side='left'
			elif int(SiPM)>7: side='right'
			else: return
		elif s==3:
			s,p,b=self.parseDetID(int(detID))
			if p!=3 and b<60 and int(SiPM)==0: side='left'
			elif p!=3 and b<60 and int(SiPM)==1: side='right'
			elif p!=3 and b>59 and int(SiPM)==0: side='top'
			elif p==3 and int(SiPM)==0: side='top'
			else:
				print('huh?')
				return

		return side

	def IsGel(self, fixed_ch):
		detID, SiPM = fixed_ch.split('_')
		s,p,b = self.parseDetID(int(detID))
		if s!=2: return 0
		if int(SiPM)>7: side='right'
		else: side='left'
		if self.gelsides[p]==side: return 1
		else: return 0

	def parseDetID(self, detID):
		if not isinstance(detID, int): detID=int(detID)
		subsystem=detID//10000
		plane=detID%10000//1000
		bar=detID%1000
		return subsystem, plane, bar

	def MakeDetID(self, fixed):
		fixed_subsystem, fixed_plane, fixed_bar = fixed
		if fixed_subsystem in (1,2,3):
			return int(f'{str(fixed_subsystem)}{str(fixed_plane)}{str(fixed_bar).zfill(3)}')

	def MakeFixedCh(self, fixed):
		fixed_subsystem, fixed_plane, fixed_bar, fixed_SiPM = fixed
		return f'{str(fixed_subsystem)}{str(fixed_plane)}{str(fixed_bar).zfill(3)}_{str(fixed_SiPM)}'

	def GetDSPlaneNumber(self, detID):
		s,p,b=self.parseDetID(detID)
		if s<3: return p
		if b>59 and p<3: plane=p*2+1
		elif any( [b<60, b>59 and p==3] ): plane=p*2
		return plane

	def MakeHumanReadableFixedCh(self, fixed_ch):
		s,p,b=self.parseDetID(int(fixed_ch.split('_')[0]))
		SiPM=int(fixed_ch.split('_')[1])
  		# Not using f-string because ROOT can't cope with combining them and TLatex
		res='US plane '+str(p+1)+', bar '+str(b+1)+', SiPM '+str(SiPM+1)+', '+self.GetSide(fixed_ch)+' side' # +1 to make things more readable. I hope this doesn't complicate things
		return res

	def MakeHumanReadableDetID(self, detID):
		s,p,b=self.parseDetID(detID)
  		# Not using f-string because ROOT can't cope with combining them and TLatex
		res=str(self.subsystemNames[s])+', plane '+str(p+1)+', bar '+str(b+1)
		return res

	def GetScifiAverageTime(self, scifi, scifihits, mode='tsf0'):

		stations = {i.GetDetectorID():i for i in scifihits}
		times=[]

		if mode=='tsf0':
			for hit in scifihits:
				detID = hit.GetDetectorID()
				t = scifi.GetCorrectedTime(detID, hit.GetTime(), 0)
				times.append(t)
			if len(times)==0: return -999
			else: return sum(times) / len(times)

		elif mode=='deltastations':
			nstations = scifi.GetConfParI('Scifi/nscifi')
			times={st:[] for st in range(nstations+1)}
			res={f'scifi-delta{i+1}{i}':0 for i in range(nstations)}
			for hit in scifihits:
				detID = hit.GetDetectorID()
				station = hit.GetStation()
				t = scifi.GetCorrectedTime(detID, hit.GetTime(), 0)
				times[station].append(t)
			for station in range(nstations-1):
				if len(times[station+1]) == 0 or len(times[station]) == 0: continue
				res[f'scifi-delta{station+1}{station}'] = sum(times[station+1])/len(times[station+1]) - sum(times[station])/len(times[station])
			return res

	def GetScifiTrackAverageTime(self, scifi, scifihits):

		stations = {i.GetDetectorID():i for i in scifihits}
		times=[]

		track = self.task.track
		nM = track.getNumPointsWithMeasurement()

		for n in range(nM):
			M=track.getPointWithMeasurement(n)
			W=M.getRawMeasurement()
			detID=W.getDetId()

			# hkey=W.getHitId()
			if detID not in stations:
				print(stations)
				print(f'Event {self.task.M.EventNumber}')
			trackHit = stations[detID]
			time = scifi.GetCorrectedTime(detID, trackHit.GetTime(), 0)
			times.append(time)

		if len(times)==0: return
		averagetime = sum(times) / len(times)
		return averagetime

	def GetPlaneData(self, hits, mufilter, **kwargs):
		planewise_data = {i:{} for i in range(5)}
		if "hit_times" in kwargs: hit_times=kwargs.get("hit_times", None)

		# Group hits by plane
		skipped_hits=0
		for hit in hits:

			detID = hit.GetDetectorID()
			s,p,b = self.parseDetID(detID)
			if not s==2:
				continue

			if not hit.isValid():
				# print(f'Hit not valid')
				continue

			# Cutting delayed MC hits
			# if self.simulation:
			# 	if detID not in hit_times:
			# 		# print(f'No hit times for detID {detID} in simulation event')
			# 		continue
			# 	delay = hit_times[detID]['delay_vs_earliest_tofcorr']
			# 	if delay > 5:
			# 		# print(f'Skipping hit for delay {delay} ns')
			# 		skipped_hits += 1
			# 		continue

			# muAna knows the reference time needed to apply for !simulation
			if not self.simulation: alignedtimes=self.GetCorrectedTimes(hit, mode='aligned', MuFilter=mufilter)
			else: alignedtimes=hit.GetAllTimes()

			atimes_left = [i[1] for i in alignedtimes if i[0]<8]
			atimes_right = [i[1] for i in alignedtimes if i[0]>=8]

			if len(atimes_left)==0 or len(atimes_right)==0: continue

			# Check nSiPMs left and right. Skip hit if less than 4 SiPMs on either left or right fire
			if self.nuStudy and (len(atimes_left)<4 or len(atimes_right)<4):
				print(f'skipping event for nSiPMs criteria: {len(atimes_left)}, {len(atimes_right)}')
				continue

			atimes_left_mean = sum(atimes_left)/len(atimes_left)
			atimes_right_mean = sum(atimes_right)/len(atimes_right)

			atimes_left_median = np.median(atimes_left)
			atimes_right_median = np.median(atimes_right)

			# Skip hit if abs( atimes_left(right) ) > L/2 / cscint_L(R)
			averagecscint_left, averagecscint_right = self.GetBarAveragecscint(mufilter, detID)

			if self.nuStudy and (abs(atimes_left_mean) > self.barlengths[2]/2 / averagecscint_left[0] or abs(atimes_right_mean) > self.barlengths[2]/2 / averagecscint_right[0]):
				print(f'skipping event for times criteria')
				continue

			planewise_data[p][detID] = {}
			planewise_data[p][detID]['bar-QDC'] = self.GetTotalQDC(hit.GetAllSignals())

			planewise_data[p][detID]['atimes-left-mean'] = atimes_left_mean
			planewise_data[p][detID]['atimes-right-mean'] = atimes_right_mean

			planewise_data[p][detID]['atimes-left-median'] = atimes_left_median
			planewise_data[p][detID]['atimes-right-median'] = atimes_right_median

			planewise_data[p][detID]['cscint-left'] = averagecscint_left
			planewise_data[p][detID]['cscint-right'] = averagecscint_right

			# Find the time that corresponds to the minimum
			mufilter.GetPosition(detID, self.A, self.B)
			x_ref = 0.5 * (self.A.x() + self.B.x())
			left_edge = self.A.x()
			right_edge = self.B.x()

			# Find the distances to the left and right edges of the bar
			distance_to_left_edge = {i: (left_edge - (x_ref + i * averagecscint_left[0])) for i in atimes_left}
			distance_to_right_edge = {i: (right_edge - (x_ref - i * averagecscint_right[0])) for i in atimes_right}

			# Find the time in atimes_left that produces distance_to_left_edge closest to zero
			min_left = min(distance_to_left_edge, key=lambda k: abs(distance_to_left_edge[k]))

			# Find the time in atimes_right that produces distance_to_right_edge closest to zero
			min_right = min(distance_to_right_edge, key=lambda k: abs(distance_to_right_edge[k]))

			planewise_data[p][detID]['atimes-left-min'] = min_left
			planewise_data[p][detID]['atimes-right-min'] = min_right

		# if self.simulation: print(f'Skipped {skipped_hits}/{len(hits)}')
		return planewise_data

	def GetBarycentres(self, hits, **kwargs):

		"""
		Adding some kwargs for using when the analysis instance
		isn't connected to a FairTask
		"""

		mufilter=kwargs.get("MuFilter")
		if "hasTrack" in kwargs: hasTrack=kwargs.get("hasTrack")
		else: hasTrack=False

		if self.simulation:
			hit_times = kwargs.get("hit_times", None)
		# 	if hit_times == None:
		# 		print('Warning: no hit times provided for simulation event!')
		# 		return

		if not hasattr(self, "barlengths"): self.BuildBarLengths(mufilter)

		barycentres={i:{} for i in range(5)}

		# Plane data rejects invalid hits!
		# For simulation, also rejects delayed hits
		if self.simulation: planewise_data = self.GetPlaneData(hits, mufilter, hit_times=hit_times)
		else: planewise_data = self.GetPlaneData(hits, mufilter)

		for plane in planewise_data:
			pdata = planewise_data[plane]
			if len(pdata)==0: continue

			# Determine quantities for the bars now
			planeQDC = sum([ pdata[detID]['bar-QDC'] for detID in pdata ])

			weighted_ys = []
			y_weights = []
			x_barycentres = {}

			for detID in pdata:
				# Get weighted y-position for each hit that passes selection in this plane
				s,p,b = self.parseDetID(detID)

				if hasTrack:
					if self.GetExtrapolatedBarDetID(p) == detID: trackInBar = True
					else: trackInBar=False
				else: trackInBar=False

				barQDC = pdata[detID]['bar-QDC']
				mufilter.GetPosition(detID, self.A, self.B)
				y_pos = barQDC/planeQDC * 0.5 * (self.A.y() + self.B.y())
				x_midpoint = 0.5 * (self.A.x() + self.B.x())
				weighted_ys.append(y_pos) # Get weighted y-pos
				y_weights.append(barQDC/planeQDC)

				# Get x-barycentre
				cscint_left, cscint_right = pdata[detID]['cscint-left'], pdata[detID]['cscint-right']

				# Get the average time for the left and right sides
				atimes_left_mean, atimes_right_mean = pdata[detID]['atimes-left-mean'], pdata[detID]['atimes-right-mean']
				# Get the minimum time for the left and right sides
				atimes_left_min, atimes_right_min = pdata[detID]['atimes-left-min'], pdata[detID]['atimes-right-min']
				atimes_left_median, atimes_right_median = pdata[detID]['atimes-left-median'], pdata[detID]['atimes-right-median']

				# Calculate the x-positions in the physics FoR for the average time
				xLphys_mean, xRphys_mean = (x_midpoint+atimes_left_mean*cscint_left[0]), (x_midpoint-atimes_right_mean*cscint_right[0])
				# Calculate the x-positions in the physics FoR for the median time
				xLphys_median, xRphys_median = (x_midpoint+atimes_left_median*cscint_left[0]), (x_midpoint-atimes_right_median*cscint_right[0])
				# Calculate the x-positions in the physics FoR with the minimum time
				xLphys_min, xRphys_min = (x_midpoint+atimes_left_min*cscint_left[0]), (x_midpoint-atimes_right_min*cscint_right[0])

				# Calculate the uncertainties in the x-positions
				xLphys_err, xRphys_err = self.Getxuncertainty(detID, pdata[detID], 'left'), self.Getxuncertainty(detID, pdata[detID], 'right')

				# Store barycentre determined by each bar
				x_barycentres[detID] = {
					'xL-mean':(xLphys_mean, xLphys_err),
					'xR-mean':(xRphys_mean, xRphys_err),

					'xL-median':(xLphys_median, xLphys_err),
					'xR-median':(xRphys_median, xRphys_err),

					'xL-mintime':(xLphys_min, xLphys_err),
					'xR-mintime':(xRphys_min, xRphys_err),

					'barQDC':barQDC
					# "trackInBar":shower_side,
					}

			y_barycentre=sum(weighted_ys)
			d_y_barycentre = 6/np.sqrt(12) * np.sqrt(sum([w**2 for w in y_weights])) # Uncertainty on y-barycentre

			mufilter.GetPosition(max(pdata.keys()), self.A, self.B)
			max_y = 0.5*(self.A.y() + self.B.y())
			mufilter.GetPosition(min(pdata.keys()), self.A, self.B)
			min_y = 0.5*(self.A.y() + self.B.y())
			lambda_y = max_y - min_y

			y_barycentres = {
				'yB':y_barycentre,
				'dyB':d_y_barycentre,
				'lambda_y':lambda_y
				}

			barycentres[plane] = {'x-barycentres':x_barycentres, "y-barycentre":y_barycentres}

		return barycentres

	def GetOverallXBarycentre(self, barycentres, mode):
		xs = {p:{} for p in barycentres.keys()}

		if mode=='relQDC':
			for p,pdata in barycentres.items():
				if len(pdata)==0: continue
				xb_data = pdata['x-barycentres']

				# Calculate the relative QDCs for each bar in the plane
				planeQDC = sum([xb_data[detID]['barQDC'] for detID in xb_data])
				relQDCs={detID:xb_data[detID]['barQDC']/planeQDC for detID in xb_data.keys()}

				xLphys_mean = sum([relQDCs[detID]*xb_data[detID]['xL-mean'][0] for detID in xb_data.keys()])
				xLphys_median = sum([relQDCs[detID]*xb_data[detID]['xL-median'][0] for detID in xb_data.keys()])
				xLphys_mintime = sum([relQDCs[detID]*xb_data[detID]['xL-mintime'][0] for detID in xb_data.keys()])
				sigma_xLphys = np.sqrt(sum([relQDCs[detID]**2*xb_data[detID]['xL-mean'][1]**2 for detID in xb_data.keys()]))

				xRphys_mean = sum([relQDCs[detID]*xb_data[detID]['xR-mean'][0] for detID in xb_data.keys()])
				xRphys_median = sum([relQDCs[detID]*xb_data[detID]['xR-median'][0] for detID in xb_data.keys()])
				xRphys_mintime = sum([relQDCs[detID]*xb_data[detID]['xR-mintime'][0] for detID in xb_data.keys()])
				sigma_xRphys = np.sqrt(sum([relQDCs[detID]**2*xb_data[detID]['xR-mean'][1]**2 for detID in xb_data.keys()]))

				# xB = sum([relQDCs[detID]*xb_data[detID]['xB'] for detID in xb_data.keys()])

				# lambda_x = sum([relQDCs[detID]*xb_data[detID]['lambda_x'] for detID in xb_data.keys()])

				xs[p]['xL-mean']=(xLphys_mean, sigma_xLphys)
				xs[p]['xR-mean']=(xRphys_mean, sigma_xRphys)
				xs[p]['xL-median']=(xLphys_median, sigma_xLphys)
				xs[p]['xR-median']=(xRphys_median, sigma_xRphys)
				xs[p]['xL-mintime']=(xLphys_mintime, sigma_xLphys)
				xs[p]['xR-mintime']=(xRphys_mintime, sigma_xRphys)

				# xs[p]['xB-mean']=xB
				# xs[p]['lambda_x-mean']=lambda_x
			return xs

		elif mode=='maxQDC':
			for p,pdata in barycentres.items():
				if len(pdata)==0: continue
				xb_data=pdata['x-barycentres']

				max_key = max(xb_data, key=lambda k: xb_data[k]['barQDC'])

				xL_mean=xb_data[max_key]['xL-mean'][0]
				xR_mean=xb_data[max_key]['xR-mean'][0]
				xL_median=xb_data[max_key]['xL-median'][0]
				xR_median=xb_data[max_key]['xR-median'][0]
				xL_mintime=xb_data[max_key]['xL-mintime'][0]
				xR_mintime=xb_data[max_key]['xR-mintime'][0]
				sigma_xLphys = xb_data[max_key]['xL-mean'][1]
				sigma_xRphys = xb_data[max_key]['xR-mean'][1]

				xs[p]['xL-mean'] =(xL_mean, sigma_xLphys)
				xs[p]['xR-mean'] = (xR_mean, sigma_xRphys)
				xs[p]['xL-median'] = (xL_median, sigma_xLphys)
				xs[p]['xR-median'] = (xR_median, sigma_xRphys)
				xs[p]['xL-mintime'] = (xL_mintime, sigma_xLphys)
				xs[p]['xR-mintime'] = (xR_mintime, sigma_xRphys)
			return xs

	def Getxuncertainty(self, detID, bardata, side):

		barside_sigma = self.CalculateBarsideTimeresolution(self.runNr, detID, side)

		err_side  = np.sqrt( (bardata[f'atimes-{side}-mean'] * bardata[f'cscint-{side}'][1])**2 + ( bardata[f'cscint-{side}'][0] * barside_sigma)**2)

		return err_side

	def GetDSHaverage(self, hits, mode='tds0'):

		"""
		For tds0, I am only using DS2 and 3 due to the delta tds21 plot being broad and asymmetric... not implemented, should I?
		"""
		stations = {i.GetDetectorID():i for i in hits}

		if mode=='tds0':
			total={i:0 for i in range(4)}
			counter={i:0 for i in range(4)}
			theTrack = self.task.M.Reco_MuonTracks[0]
			nM = theTrack.getNumPointsWithMeasurement()

			for n in range(nM):
				M = theTrack.getPointWithMeasurement(n)
				W = M.getRawMeasurement()
				detID = W.getDetId()
				s,p,b = self.parseDetID(detID)
				hkey = W.getHitId()
				trackHit = stations[detID]
				tdcs = trackHit.GetAllTimes()
				if not all ( [p in (0, 1, 2), b<60, len(tdcs)==2] ): continue

				if len(tdcs) != 2: continue # pass

				for item in tdcs:
					SiPM, clock = item
					dscorrectedtime=self.task.M.MuFilter.GetCorrectedTime(detID, SiPM, clock*self.TDC2ns, 0)
					total[p]+=dscorrectedtime
					counter[p]+=1

			# for hit in hits:

			# 	detID = hit.GetDetectorID()
			# 	s,p,b = self.parseDetID(detID)
			# 	if not all ( [s==3, p in (0,1,2), b<60] ):
			# 		continue

			# 	if not hit.isValid():
			# 		# print(f'Hit not valid')
			# 		continue

			# 	tdcs = hit.GetAllTimes()
			# 	if len(tdcs) != 2: continue # pass

			# 	for item in tdcs:
			# 		SiPM, clock = item
			# 		dscorrectedtime=self.task.M.MuFilter.GetCorrectedTime(detID, SiPM, clock*self.TDC2ns, 0)
			# 		total[p]+=dscorrectedtime
			# 		counter[p]+=1

			# Only using DS3 (p=2) and DS2 (p=1)
			sum_total, counter_total = total[1] + total[2], counter[1] + counter[2]
			# print(f'total, counter: {sum_total}, {counter_total}')
			if counter_total==0: return -999

			dsh_average, nfired = sum_total/counter_total, counter_total

			return dsh_average, nfired

		elif mode=='timingdiscriminant': # Take k=2 for the 3rd horizontal plane

			total=[]

			theTrack = self.task.M.Reco_MuonTracks[0]
			nM = theTrack.getNumPointsWithMeasurement()
			# print(f'{nM} measurements')
			for n in range(nM):
				M=theTrack.getPointWithMeasurement(n)
				W=M.getRawMeasurement()
				detID=W.getDetId()
				s,p,b = self.parseDetID(detID)
				hkey=W.getHitId()
				trackHit = stations[detID]
				tdcs = trackHit.GetAllTimes()
				# print(f'DetID: {detID}, len(tdcs)={len(tdcs)}')
				if not all ( [p in (0,1,2), b<60, len(tdcs)==2] ): continue

				if len(tdcs) != 2: continue # pass

				for item in tdcs:
					SiPM, clock = item
					dscorrectedtime=self.task.M.MuFilter.GetCorrectedTime(detID, SiPM, clock*self.TDC2ns, 0)
					total.append(dscorrectedtime)

			if len(total)==0:
				return -420
			return sum(total)/len(total)

		elif mode=='deltastations':
			res={f'delta32':0, 'delta21':0}
			ts={k:[] for k in range(3)}

			# theTrack = self.task.M.Reco_MuonTracks[0]
			# nM = theTrack.getNumPointsWithMeasurement()
			# # print(f'{nM} measurements')
			# for n in range(nM):
			# 	M=theTrack.getPointWithMeasurement(n)
			# 	W=M.getRawMeasurement()
			# 	detID=W.getDetId()
			# 	s,p,b = self.parseDetID(detID)
			# 	hkey=W.getHitId()
			# 	trackHit = stations[detID]
			# 	tdcs = trackHit.GetAllTimes()
			# 	# print(f'DetID: {detID}, len(tdcs)={len(tdcs)}')

			for hit in hits:
				detID = hit.GetDetectorID()
				s,p,b = self.parseDetID(detID)
				tdcs = hit.GetAllTimes()

				# criteria to selection only horizontal hits with both SiPMs firing
				if s!=3: continue
				if not all ( [p in (0,1,2), b<60, len(tdcs)==2] ): continue

				for item in tdcs:
					SiPM, clock = item
					dscorrectedtime=self.task.M.MuFilter.GetCorrectedTime(detID, SiPM, clock*self.TDC2ns, 0)
					ts[p].append(dscorrectedtime)

			if any([len(ts[i])==0 for i in ts]): return False

			res['t3']=sum(ts[2])/len(ts[2])
			res['t2']=sum(ts[1])/len(ts[1])
			res['t1']=sum(ts[0])/len(ts[0])

			return res

		elif mode=='testing-tds0':
			total={i:0 for i in range(4)}
			counter={i:0 for i in range(4)}
			theTrack = self.task.M.Reco_MuonTracks[0]
			for nM in range(theTrack.getNumPointsWithMeasurement()):
				M=theTrack.getPointWithMeasurement(nM)
				W=M.getRawMeasurement()
				detID=W.getDetId()
				s,p,b = self.parseDetID(detID)
				hkey=W.getHitId()
				# trackHit=hits[hkey]
				trackHit = stations[detID]
				tdcs = trackHit.GetAllTimes()
				if not all ( [p in (0,1,2), b<60, len(tdcs)!=2] ): continue

				if len(tdcs) != 2: continue # pass

				for item in tdcs:
					SiPM, clock = item
					dscorrectedtime=self.task.M.MuFilter.GetCorrectedTime(detID, SiPM, clock*self.TDC2ns, 0)
					total[p]+=dscorrectedtime
					counter[p]+=1

			# Only using DS3 and DS2
			sum_total, counter_total = total[2] + total[1], counter[2] + counter[1]
			if counter_total==0: return -999

			dsh_average, nfired = sum_total/counter_total, counter_total

			return dsh_average, nfired

		elif mode=='testing-deltastations':
			res={f'delta32':0, 'delta21':0}
			ts={k:[] for k in range(3)}
			theTrack = self.task.M.Reco_MuonTracks[0]
			for nM in range(theTrack.getNumPointsWithMeasurement()):
				M=theTrack.getPointWithMeasurement(nM)
				W=M.getRawMeasurement()
				detID=W.getDetId()
				s,p,b = self.parseDetID(detID)
				if not all ( [p in (0,1,2), b<60] ): continue

				hkey=W.getHitId()
				trackHit = stations[detID]
				times=trackHit.GetAllTimes()
				if not len(times)==2: continue

				for idx in times:
					SiPM, clock=idx
					dscorrectedtime=self.task.M.MuFilter.GetCorrectedTime(detID, SiPM, clock*self.TDC2ns, 0)
					ts[p].append(dscorrectedtime)

			if any([len(ts[i])==0 for i in ts]):
				# print(f'ts: {ts}')
				return False
			res['delta32']=sum(ts[2])/len(ts[2]) - sum(ts[1])/len(ts[1])
			res['delta21']=sum(ts[1])/len(ts[1]) - sum(ts[0])/len(ts[0])
			return res

	def GetTimingDiscriminant(self, hits=None, mufilter=None):
		if hits is None:
			if hasattr(self, "task") and hasattr(self.task, "M"):
				hits = self.task.M.eventTree.Digi_MuFilterHits
			else:
				return -420.0

		if mufilter is None:
			if hasattr(self, "task") and hasattr(self.task, "M"):
				mufilter = self.task.M.MuFilter
			else:
				return -420.0

		# Find the hit in US Plane 1 (subsystem=2, plane=0)
		US1hit = None
		for hit in hits:
			detID = hit.GetDetectorID()
			s, p, b = self.parseDetID(detID)
			if s == 2 and p == 0:
				US1hit = hit
				break

		if US1hit is None:
			return -420.0

		us1averagetime = self.GetAverageTime(US1hit, mufilter, correctTW=False)
		if us1averagetime is None or us1averagetime == -420:
			return -420.0

		# DS3 horizontal average time
		ds3haverage = self.GetDSHaverage(hits, mode="timingdiscriminant")
		if ds3haverage is None or ds3haverage == -420:
			return -420.0

		return ds3haverage - us1averagetime

	def fit_langau(self, hist,o,bmin,bmax):
		params = {0:'Width(scale)',1:'mostProbable',2:'norm',3:'sigma'}
		F = ROOT.TF1('langau',self.langaufun,0,200,4)
		for p in params: F.SetParName(p,params[p])
		rc = hist.Fit('landau','S'+o,'',bmin,bmax)
		res = rc.Get()
		if not res: return res
		F.SetParameter(2,res.Parameter(0))
		F.SetParameter(1,res.Parameter(1))
		F.SetParameter(0,res.Parameter(2))
		F.SetParameter(3,res.Parameter(2))
		F.SetParLimits(0,0,10)
		F.SetParLimits(1,0,100)
		F.SetParLimits(3,0,10)

		rc = hist.Fit(F,'S'+o,'',bmin,bmax)
		res = rc.Get()
		return res

	def langaufun(self, x,par):
		#Fit parameters:
		#par[0]=Width (scale) parameter of Landau density
		#par[1]=Most Probable (MP, location) parameter of Landau density
		#par[2]=Total area (integral -inf to inf, normalization constant)
		#par[3]=Width (sigma) of convoluted Gaussian function
		#
		#In the Landau distribution (represented by the CERNLIB approximation),
		#the maximum is located at x=-0.22278298 with the location parameter=0.
		#This shift is corrected within this function, so that the actual
		#maximum is identical to the MP parameter.
		#
		# Numeric constants
		invsq2pi = 0.3989422804014   # (2 pi)^(-1/2)
		mpshift  = -0.22278298       # Landau maximum location
		#
		# Control constants
		np = 100.0      # number of convolution steps
		sc =   5.0      # convolution extends to +-sc Gaussian sigmas
		#
		# Variables
		summe = 0.0
		#
		# MP shift correction
		mpc = par[1] - mpshift * par[0]
		#
		# Range of convolution integral
		xlow = x[0] - sc * par[3]
		xupp = x[0] + sc * par[3]
		#
		step = (xupp-xlow) / np
		#
		# Convolution integral of Landau and Gaussian by sum
		i=1.0
		if par[0]==0 or par[3]==0: return 9999
		while i<=np/2:
			i+=1
			xx = xlow + (i-.5) * step
			fland = ROOT.TMath.Landau(xx,mpc,par[0]) / par[0]
			summe += fland * ROOT.TMath.Gaus(x[0],xx,par[3])
			#
			xx = xupp - (i-.5) * step
			fland = ROOT.TMath.Landau(xx,mpc,par[0]) / par[0]
			summe += fland * ROOT.TMath.Gaus(x[0],xx,par[3])

		return (par[2] * step * summe * invsq2pi / par[3])

	def DSAcceptanceRange(self, MuFilter):
		barNumbers={
			'horizontal':
				{'top':'059', 'bottom':'000', 'left':'029', 'right':'029'},
			'vertical':
				{'top':'089', 'bottom':'089', 'left':'119', 'right':'060'}
			}

		vals = {'horizontal':{}, 'vertical':{}}
		res={}

		# x for horizontal planes, y for vertical planes
		for i in barNumbers.keys():
			for plane in range(4):

				if i=='horizontal' and plane==3:continue # only 4th vertical plane

				for position in barNumbers[i]:
					bar = barNumbers[i][position]
					detID=f'3{plane}{bar}'

					MuFilter.GetPosition(int(detID), self.A, self.B)

					if i=='horizontal' and position in ['top', 'bottom']:
						vals[i][position] = 1/2 * (self.A.y()+self.B.y())
					elif i=='horizontal' and position == 'left': vals[i][position] = self.A.x()
					elif i=='horizontal' and position == 'right': vals[i][position] = self.B.x()

					elif i=='vertical' and position in ['left', 'right']:
						vals[i][position] = 1/2 * (self.A.x()+self.B.x())
					elif i=='vertical' and position == 'top': vals[i][position] = self.A.y()
					elif i=='vertical' and position == 'bottom': vals[i][position] = self.A.y()

		return vals

	def GetAverageTime(self, mufiHit, mufilter, side='both', correctTW=True):
		value=[0, 0]
		count=[0, 0]
		nSiPMs=mufiHit.GetnSiPMs()
		if correctTW: times = self.GetCorrectedTimes(mufiHit, mode='aligned', MuFilter=mufilter) # in nanoseconds
		else: times = mufiHit.GetAllTimes() # in clockcycles

		detID=mufiHit.GetDetectorID()
		s, p, b = self.parseDetID(detID)
		for element in times:
			SiPM, clock = element

			# Check if this channel has determined alignment parameters
			fixed_ch = f'{detID}_{SiPM}'
			if correctTW and not fixed_ch in self.alignmentparameters: continue
			elif correctTW and fixed_ch in self.alignmentparameters: d = self.alignmentparameters[fixed_ch]
			if s==2 and self.IsSmallSiPMchannel(SiPM):continue

			if SiPM<nSiPMs:
				if not correctTW: value[0]+=clock*self.TDC2ns # in nanoseconds.
				else: value[0]+= clock - d[0] # in nanoseconds ## I can't subtract d[0] here because then I have a sign error
				# else: value[0]+= clock # in nanoseconds
				count[0]+=1
			else:
				if not correctTW: value[0]+=clock*self.TDC2ns
				else: value[0]+= clock - d[0] # in nanoseconds ## I can't subtract d[0] here because then I have a sign error
				# else: value[0]+= clock
				count[1]+=1

		if s == 2:
			if count[0] != 0 and count[1] != 0:
				if side == 'both':
					# average = 0.5*(value[0]/count[0]+value[1]/count[1])
					average = sum(value) / sum(count)
					return average
				if side == 'L':
					return value[0]/count[0]
				elif side == 'R':
					return value[1]/count[1]
				else: return -999.
			else:
				return
				# print(f'hit does not have fired SiPMs on both sides.')
		elif s == 3 and b < 60: # I don't require a SiPM to fire on both sides here because that is required by the TDS0 determination
			if side=='both':
				average = 0.5*(value[0]/count[0]+value[1]/count[1])
				return average
			elif side == 'L': return value[0]/count[0]
			elif side == 'R': return value[1]/count[1]

	# def GetMedianTime(self, hit, mode=None, particleToF=None):
	def GetMedianTime(self, hit, mode=None):

		if not mode: mode=self.state
		# if self.timealignment=='old' and particleToF==None:
		# 	print(f'A particle ToF must be provided for old time alignment data')
		# 	return -999

		values={'left':[], 'right':[]}
		nSiPMs=hit.GetnSiPMs()
		clocks=hit.GetAllTimes()
		if mode != 'uncorrected': qdcs=hit.GetAllSignals()
		detID=hit.GetDetectorID()
		s, p, b = self.parseDetID(detID)

		# Should not be needed
		if s==3 and b>59:
			vals=[i[1] for i in clocks]
			return sum(vals)/len(vals)

		nLeft, nRight = self.GetnFiredSiPMs(hit)
		nSiPMconditions={1: any( [nLeft<6, nRight<6]),
						2: any( [nLeft<4, nRight<4]),
						3: all(( b<60, any( [nLeft!=1, nRight!=1]) ))
      					}
		if nSiPMconditions[s]:
			# print(f'nSiPM conditions not met for {detID}')
			return

		for element in clocks:
			SiPM, clock = element
			fixed_ch=self.MakeFixedCh((s,p,b,SiPM))
			side=self.GetSide(fixed_ch)

			if s==2 and self.IsSmallSiPMchannel(SiPM): continue

			if mode == 'uncorrected': values[side].append(clock*self.TDC2ns)

			elif mode == 'corrected':
				qdc=self.GetChannelVal(SiPM, qdcs)
				if not qdc: continue
				correctedtime=self.correct_TW(fixed_ch, qdc, clock)
				if not correctedtime: continue
				values[side].append(correctedtime)

			elif mode == 'alignment':
				qdc=self.GetChannelVal(SiPM, qdcs)
				if not qdc: continue
				correctedtime=self.MuFilterCorrectedTime(fixed_ch, clock, qdc)
				if not correctedtime: continue

				values[side].append(correctedtime)

		if not all( (len(values[i])>0 for i in values) ):
			print(f'Zero entries on one or both sides of bar {detID}')
			return
		medians={}

		if s==3 and b<60 and all( [nLeft==1, nRight==1] ):
			for x in values: medians[x]=values[x][0]
			return medians

		for x in values:
			if len(values[x])%2==0:
				if len(values[x])==2: print(detID, values)
				medians[x] = 0.5* ( values[x][int(len(values[x])/2)] + values[x][int( len(values[x])/2+1)] )
			else:
				medians[x] = values[x][int( 0.5*(len(values[x])+1))]
		return medians

	def GetTotalQDC(self, signals):
		return sum( [i[1] for i in signals if not self.IsSmallSiPMchannel(i[0])] )

	def GetChannelVal(self, SiPM, chs):
		for entry in chs:
			fSiPM, val = entry
			if fSiPM == SiPM:
				return val
		return

	def Get_numuevents(self):
		numusignalevent_filepath = '/afs/cern.ch/work/a/aconsnd/numusignalevents.csv'
		self.nu_mu_events = {}

		with open(numusignalevent_filepath, 'r') as f:
			reader=csv.reader(f)
			# next(reader) # skip first row with the headers

			for idx,x in enumerate(reader):
				if idx==0: continue

				self.nu_mu_events[int(x[0])] = [int(x[1]), int(x[2])] + [float(i) for i in x[3:]]

	def GetNeutrinoIntType(self, event):

		if not hasattr(event, "MCTrack"):
			print(f'No MCTrack branch. Is this real data?')
			return

		if event.MCTrack[0].GetPdgCode() == event.MCTrack[1].GetPdgCode():
			i_flav = 0 #NC
		elif abs(event.MCTrack[1].GetPdgCode()) == 11:
			i_flav = 1 #nueCC
		elif abs(event.MCTrack[1].GetPdgCode()) == 13:
			i_flav = 2 #numuCC
		elif abs(event.MCTrack[1].GetPdgCode()) == 15:
			is1Mu = False
			for j_track in range(2, len(event.MCTrack)):
				if event.MCTrack[j_track].GetMotherId() == 1 and abs(event.MCTrack[j_track].GetPdgCode()) == 13:
					is1Mu = True
					break
			if is1Mu:
				i_flav = 4 #nutauCC1mu
			else:
				i_flav = 3 #nutauCC0mu

	def OneHitPerSystem(self, hits, systems, Nfired=False):
		verbose=self.verbose

		hitdict={}
		for s in systems: # systems always includes US and DS. Veto is included if a Scifi track is also formed.
			for p in range(self.systemAndPlanes[s]):
				key=10*s+p
				hitdict[key]=[]
		for hit in hits:
			if not hit.isValid(): continue
			detID=hit.GetDetectorID()
			s, p, b = self.parseDetID(detID)
			if s not in systems: continue
			key=10*s+p
			hitdict[key].append(detID)

		#### For returning fraction of planes with 1 scintillator firing.
		if Nfired:
			counter=0
			totalplanes=12 if systems==(2,3) else 14

			for key in hitdict:
				if key<30 and len(hitdict[key])==1: counter+=1
				elif key>=30 and key<33:
					if len(hitdict[key])!=2: continue
					bars=sorted([ self.parseDetID(hit)[2] for hit in hitdict[key] ]) # List of bar numbers sorted by number.
					if bars[0]<60: counter+=1
					if bars[1]>59: counter+=1
				elif key==33 and len(hitdict[key])==1:counter+=1
			return counter/totalplanes

		#### For requiring an event to have exactly 1 scintillator firing.
		for key in hitdict:

			hits=hitdict[key]
			if key<30 and len(hitdict[key]) != 1: return False

			elif key>=30 and key<33: # DS0 -> DS2
				if len(hits)!=2:
					if verbose: print(f'{key}: {hits}')
					return False
				bars=sorted([self.parseDetID(hit)[2] for hit in hits]) # List of bar numbers sorted by number.
				# [horizontal, vertical]
				if not bars[0]<60 and bars[1]>59:
					if verbose: print(f'bars fired for plane {key}: {bars}')
					return False
			elif key==33: # DS3 only vertical bars
				if not len(hits)==1:
					return False
			# if not DSVcheck(hitdict[key]): return False
		return True

	def GetSubsystemZlimits(self, subsystem, zPos):
		if subsystem==0:
			pass
		elif subsystem==1:
			zmin, zmax=self.zPos['MuFilter'][10], self.zPos['MuFilter'][11]
		elif subsystem==2:
			zmin, zmax=self.zPos['MuFilter'][20], self.zPos['MuFilter'][24]
		elif subsystem==3:
			zmin, zmax=self.zPos['MuFilter'][30], self.zPos['MuFilter'][36]
		else: return 0
		return zmin, zmax

	def GetSubsystemXYlimits(self, subsystem,geoobject):
		if subsystem==0:
			pass
		elif subsystem==1:
			geoobject.GetPosition(10000, self.A, self.B)
			xmin, xmax = self.B.x(), self.A.x()
			ymin = self.A.y() - geoobject.GetConfParF('MuFilter/VetoBarY')/2.
			geoobject.GetPosition(10006, self.A, self.B)
			ymax = self.A.y() + geoobject.GetConfParF('MuFilter/VetoBarY')/2.
		elif subsystem==2:
			geoobject.GetPosition(20000, self.A, self.B)
			xmin, xmax = self.B.x(), self.A.x()
			ymin = self.A.y() - geoobject.GetConfParF('MuFilter/UpstreamBarY')/2.
			geoobject.GetPosition(20009, self.A, self.B)
			ymax = self.A.y() + geoobject.GetConfParF('MuFilter/UpstreamBarY')/2.
		elif subsystem==3:
			pass
		return xmin, xmax, ymin, ymax

	def getNUSPlanes(self, hits):

		res=0

		USPlanes={k:0 for k in range(5)}
		for i, hit in enumerate(hits):
			if not hit.isValid():continue
			detID=hit.GetDetectorID()
			s, p, b=self.parseDetID(detID)
			if s != 2: continue
			USPlanes[p]+=1
		for plane in range(5):
			if USPlanes[plane]==1: res+=1
		return res

	def delta_min_t(self, aHit):
		times = aHit.GetAllTimes()
		if len(times)==0: return -999.
		nSiPMs = aHit.GetnSiPMs()
		ts_L, ts_R = [], []
		for channel in times:
			SiPM, time = channel
			if SiPM<nSiPMs: ts_L.append(time)
			elif SiPM>=nSiPMs: ts_R.append(time)
		return min(ts_L)-min(ts_R)

	def AllLiveSiPMs(self, hit):
		detID = hit.GetDetectorID()
		if detID != 21000: left_max=right_max=6
		else: left_max, right_max = 6, 5

		nleft, nright = self.GetnFiredSiPMs(hit)
		if nleft==left_max and nright==right_max: return True
		else: return False

	def GetnFiredSiPMs(self, hit):

		s,p,b=self.parseDetID(hit.GetDetectorID())
		nSiPMs=hit.GetnSiPMs()

		nFiredSiPMs_left=0
		nFiredSiPMs_right=0
		channels=hit.GetAllSignals()
		for ch in channels:
			SiPM, qdc = ch
			if s==2 and self.IsSmallSiPMchannel(SiPM):continue
			if SiPM<nSiPMs: nFiredSiPMs_left+=1
			elif SiPM>=nSiPMs: nFiredSiPMs_right+=1
		return nFiredSiPMs_left, nFiredSiPMs_right

	def GetFiredSiPMsOnPCBs(self, hits):
		channels_on_PCB = {}
		for hit in hits:
			detID = hit.GetDetectorID()
			s,p,b = self.parseDetID(detID)
			Fired_left, Fired_right = self.GetnFiredSiPMs(hit)
			fired = {'left':Fired_left, 'right':Fired_right}

			for side in fired:
				key=f'{p}_{side}'
				if not key in channels_on_PCB: channels_on_PCB[key]=0
				channels_on_PCB[key] += fired[side]

		return channels_on_PCB

	def GetSiPMNumberInSystem_LandR(self, detID, SiPM): # 20000 SiPM 8 -> 8
		if not isinstance(SiPM, int): SiPM=int(SiPM)
		s, p, b = self.parseDetID(int(detID))
		if s==1:
			nSiPMs, SiPMs_plane=16, 112
			return int(SiPM)+nSiPMs*b+p*SiPMs_plane
		elif s==2:
			nSiPMs, SiPMs_plane=16, 160
			return SiPM+nSiPMs*b+p*SiPMs_plane

		elif s==3: # Count left and right horizontal SiPMs consecutively
			nSiPMs_bar_hor, nSiPMs_bar_ver=2, 1
			nSiPMs_plane_hor, nSiPMs_plane_ver=120, 60
			if b>59 and p<3: # First 3 vertical layers
				horizontalSiPMs=(p+1)*nSiPMs_plane_hor
				verticalSiPMs=p*nSiPMs_plane_ver
				total=horizontalSiPMs+verticalSiPMs+(b-60)*nSiPMs_bar_ver
			elif b<60: # All horizontal layers
				horizontalSiPMs=p*nSiPMs_plane_hor
				verticalSiPMs=p*nSiPMs_plane_ver
				total=horizontalSiPMs+verticalSiPMs+b*nSiPMs_bar_hor+SiPM
			elif b>59 and p==3: # Final vertical layer
				total=3*120+3*60+(bar-60)
		return total

	def GetSiPMNumberInSystem_PCBbyPCB(self, detID, SiPM):
		if not isinstance(detID, int): detID=int(detID)
		if not isinstance(SiPM, int): SiPM=int(SiPM)

		s,p,b=self.parseDetID(detID)
		if s==2:
			nSiPMs=8
			SiPMs_plane=160
			total=b*nSiPMs+p*SiPMs_plane
			if SiPM<nSiPMs:
				return int(total+SiPM )
			else:
				return int(total+SiPMs_plane/2+(SiPM-nSiPMs))
		elif s==1:
			nSiPMs=8
			SiPMs_plane=112
			total=b*nSiPMs+p*SiPMs_plane
			if SiPM<nSiPMs:
				return int(total+SiPM )
			else:
				return int(total+SiPMs_plane/2+(SiPM-nSiPMs))

		elif s==3:
			nSiPMs_bar_hor, nSiPMs_bar_ver=2,1
			nSiPMs_plane_hor, nSiPMs_plane_ver=120, 60
			if b>59 and p!=3: # First 3 vertical layers
				horizontalSiPMs=int((p+1)*nSiPMs_plane_hor)
				verticalSiPMs=int(p*nSiPMs_plane_ver)
				total=horizontalSiPMs+verticalSiPMs+int((b-60)*nSiPMs_bar_ver)
				return total
			elif b<60: # All horizontal layers
				horizontalSiPMs=int(p*nSiPMs_plane_hor)
				verticalSiPMs=int(p*nSiPMs_plane_ver)
				if SiPM==0:
					total=horizontalSiPMs+verticalSiPMs+b
					return total
				elif SiPM==1:
					total=horizontalSiPMs+verticalSiPMs+b+60
					return total
			elif b>59 and p==3:
				total=3*120+3*60+(b-60)
				return int(total)

	def GetSiPMNumberInPlane_LandR(self, detID, SiPM):
		s, p, b = self.parseDetID(detID)
		if s == 1: return SiPM*b*12
		if s == 2:  return SiPM+b*16

	def GetSiPMNumberInPlane_LTR(detID, SiPM):
		s, p, b = self.parseDetID(detID)
		if s != 2:
			print('AAAAAAAHHHHHH')
		return SiPM+b*8

	def GetSiPMNumberInSystem_LTR(self, detID, SiPM): # 20000 SiPM 8 -> 400
		if not isinstance(SiPM, int): SiPM=int(SiPM)
		s, p, b = self.parseDetID(detID)
		if s==1: nSiPMs, SiPMs_plane=6, 52 # is it?
		elif s==2: nSiPMs, SiPMs_plane=8, 80
		elif s==3: nSiPMs, SiPMs_plane=1, 60 # wrong

		if SiPM<nSiPMs:
			return SiPM+nSiPMs*b+p*SiPMs_plane
		elif SiPM>=nSiPMs:
			SiPM_r=400+SiPM%nSiPMs
			return SiPM_r+nSiPMs*b+p*SiPMs_plane

	def SiPM2BarAndPosition(self, SiPM):
		# Pass the SiPM number on the PCB
		# returns the bar number and SiPM number within the bar
		barNumber = (SiPM-1)//8+1
		SiPMNumber = (SiPM-1)%8+1
		if SiPMNumber==0: SiPMNumber==8
		return barNumber, SiPMNumber

	def GetDeltaT(self, times, one_channel=None):
		# nSiPMs=aHit.GetnSiPMs()
		nSiPMs=8
		mean = [0,0]
		count = [0,0]
		# channels = aHit.GetAllTimes()
		for ch in times:
			SiPM, val = ch
			if one_channel != None:
				if not (SiPM == one_channel or SiPM == one_channel+nSiPMs): continue
			if self.IsSmallSiPMchannel(SiPM): continue
			if SiPM < nSiPMs:
				mean[0]+=val
				count[0]+=1
			else:
				mean[1]+=val
				count[1]+=1
		if count[0] != 0 and count[1] != 0:
			return (mean[0]/count[0]-mean[1]/count[1])/2.
		else: return

	def GetExtrapolatedBarDetID(self, plane):
		xEx, yEx, zEx = self.GetExtrapolatedPosition(plane)
		try:
			node = self.task.nav.FindNode(xEx, yEx, zEx)
			# if node:
			# 	print(f"Found node: {node.GetName()}")
			# else:
			# 	print("Warning: No node found at that position!")
		except Exception as e:
			print(f"Exception occurred during FindNode({xEx}, {yEx}, {zEx}): {e}")
		detIDEx = self.task.nav.FindNode(xEx, yEx, zEx).GetName()
		if not detIDEx.split('_')[0] == 'volMuUpstreamBar':
			return
		else: return int(detIDEx.split('_')[1])

	def GetExtrapolatedPosition(self, plane):
		if not self.task.hasTrack:
			print(f'No track to evaluate expected position!')
			return
		zEx = self.zPos['MuFilter'][20+plane]
		lam = (zEx - self.task.pos.z())/self.task.mom.z()

		xEx = self.task.pos.x() + lam*self.task.mom.x()
		yEx = self.task.pos.y() + lam*self.task.mom.y()

		return (xEx, yEx, zEx)

	def Getcscint(self, runNr, fixed_ch, state):

		filename=f'{self.path}cscintvalues/run{runNr}/cscint_{fixed_ch}.json'

		if not os.path.exists(filename): return

		with open(filename, 'r') as x:
			d=json.load(x)

		if state not in d: return

		if len(d[state])==6: return d[state][0], d[state][1]
		elif len(d[state])==7: return d[state][0], d[state][1], d[state][-1]

	def Getcscint_offset(self, runNr, fixed_ch, state):

		filename=f'{self.path}cscintvalues/run{runNr}/cscint_{fixed_ch}.json'

		if not os.path.exists(filename): return

		with open(filename, 'r') as x:
			d=json.load(x)

		if state not in d: return

		return d[state][2], d[state][3]

	def GetBarAveragecscint(self, MuFilter, detID):

		subsystem, plane, bar = self.parseDetID(detID)
		cscintvalues={'left':[], 'right':[]}
		for SiPM in self.systemAndSiPMs[subsystem]:
			fixed_ch=self.MakeFixedCh((subsystem, plane, bar, SiPM))
			side='left' if SiPM<8 else 'right'
			if not self.simulation:
				if fixed_ch not in self.cscintvalues:continue
				cscint=self.cscintvalues[fixed_ch]
			else:
				cscint = MuFilter.GetConfParF(f'MuFilter/US_signalspeed_{fixed_ch}')
				if cscint==0:continue
				cscint=(cscint,0) # at the moment, no uncertainty for signal speed in simulation
			if not cscint: continue
			cscintvalues[side].append(cscint)

		average_left_cscint, average_right_cscint = sum( [ci[0] for ci in cscintvalues['left']] ) / len(cscintvalues['left']), sum( [ci[0] for ci in cscintvalues['right']] ) / len(cscintvalues['right'])
		uncertainty_sq_left, uncertainty_sq_right = sum( [ci[1]**2 for ci in cscintvalues['left']] ), sum( [ci[1]**2 for ci in cscintvalues['right']] )
		uncertainty_left, uncertainty_right = ROOT.TMath.Sqrt(uncertainty_sq_left), ROOT.TMath.Sqrt(uncertainty_sq_right)

		return (average_left_cscint, uncertainty_left), (average_right_cscint, uncertainty_right)

	def Getcscint_subranges(self, runNr, fixed_ch, state):

		filename=f'{self.path}cscintvalues/run{runNr}/subrange-cscints_{fixed_ch}.json'

		if not os.path.exists(filename): return

		with open(filename, 'r') as x:
			d=json.load(x)

		if state not in d: return
		return d[state]

	def GetBarAveragesigmat(self, detID):

		subsystem, plane, bar = self.parseDetID(detID)
		bartimeresolutionvalues={}
		for SiPM in self.systemAndSiPMs[subsystem]:
			fixed_ch=self.MakeFixedCh((subsystem, plane, bar, SiPM))
			if fixed_ch not in self.timeresolutionvalues:continue
			timeresolutionvalue=self.timeresolutionvalues[fixed_ch]
			if m.isnan(timeresolutionvalue[0]): continue

			bartimeresolutionvalues[fixed_ch]=timeresolutionvalue

		average_timeresolution = sum( [bartimeresolutionvalues[ch][0] for ch in bartimeresolutionvalues] ) / len(bartimeresolutionvalues.items())
		uncertainty_sq = sum( [bartimeresolutionvalues[ch][1]**2 for ch in bartimeresolutionvalues] )
		uncertainty=ROOT.TMath.Sqrt(uncertainty_sq)

		return average_timeresolution, uncertainty

	# def GetBarAverageQDC(self, qdcs):
		# for

	def Getcscint_chi2pNDF_info(self, runNr,fixed_ch,state):
		iteration=0 if state=='uncorrected' else 1
		with open(f'{self.path}cscintvalues/run{runNr}/cscint_{fixed_ch}.csv', 'r') as handle:
			reader=csv.reader(handle)
			alldata=[row for row in reader]
			if len(alldata)<iteration+1: return
			try:
				data=alldata[iteration]
			except IndexError:
				print(f'{fixed_ch} IndexError')
		return float(data[-2]),int(data[-1])

	def Getcscint_chi2pNDF(self, runNr,fixed_ch,state):
		iteration=0 if state=='uncorrected' else 1
		if not os.path.exists(f'{self.path}cscintvalues/run{runNr}/cscint_{fixed_ch}.csv'): return
		with open(f'{self.path}cscintvalues/run{runNr}/cscint_{fixed_ch}.csv', 'r') as handle:
			reader=csv.reader(handle)
			alldata=[row for row in reader]
			if len(alldata)<iteration+1: return
			try:
				data=alldata[iteration]
			except IndexError:
				print(f'{fixed_ch} IndexError')
		return float(data[-2])/int(data[-1])

	def Makecscintdict(self, runNr, state='corrected'):

		if self.timealignment=='old': run=str(5097).zfill(6)
		elif self.timealignment=='new': run=str(5408).zfill(6)
		elif self.timealignment=='new+LHCsynch': run=str(5999).zfill(6)

		d={}
		s=2
		for p in range(self.systemAndPlanes[s]):
			for b in range(self.systemAndBars[s]):
				for SiPM in self.systemAndSiPMs[s]:
					fixed_ch=self.MakeFixedCh((s,p,b,SiPM))
					cscint=self.Getcscint(runNr, fixed_ch=fixed_ch, state=state)
					if not cscint: continue
					else: d[fixed_ch]=cscint
		if len(d)==0: self.cscintvalues=None
		self.cscintvalues=d

	def Maketimeresolutiondict(self, runNr, state, fitmode='FWHM', histkey='dtvqdc'):

		path=f'{self.path}TimeResolution/run{runNr}/'
		res={}
		# if state=='corrected' or state=='aligned': state='corrected_tDS0-tSiPMcorrected'
		all_channels=self.GetListOfChannels(2)
		for fixed_ch in all_channels:

			filename=f'{path}timeresolution_{fixed_ch}.json'
			if not os.path.exists(filename): continue

			with open(filename, 'r') as x:
				d=json.load(x)

			if histkey not in d:
				print(f'histkey {histkey} for {fixed_ch} not in d')
				continue
			if state not in d[histkey]:
				print(f'state {state} for {fixed_ch} not in d[{histkey}]')
				continue
			if type(d[histkey][state])==list:
				print(f'wtf: {filename}')
			try: res[fixed_ch] = d[histkey][state][fitmode]
			except KeyError:
				print(f'{fixed_ch} has no {state} time res for mode {fitmode}')
				continue
		self.timeresolutiondict = res

	def Writetimeresolutiondict(self, runNr, state):
		self.Maketimeresolutiondict(runNr, state)

		afs_filename = f'{self.path}Results/run{runNr}/run{runNr}_timeresolution_{state}.json'
		eos_filename = f'/eos/home-a/aconsnd/SWAN_projects/Data analysis/data/run{runNr}_timeresolution_{state}.json'

		# Writing JSON data to a file
		for filename in [afs_filename, eos_filename]:
			# Ensure the directory exists
			os.makedirs(os.path.dirname(filename), exist_ok=True)
			with open(filename, 'w') as json_file:
				json.dump(self.timeresolutiondict, json_file, indent=4)
			print(f'Time resolution dict written to {filename}')

	def MakeDeltatimeresolutiondict(self, runNr, state, fitmode='FWHM', histkey='dtvxpred'):
		path=f'{self.path}TimeResolution/run{runNr}/'
		self.deltatimeresolutiondict={}
		# if state=='corrected' or state=='aligned': state='corrected_tDS0-tSiPMcorrected'
		all_channels=self.GetListOfChannels(2)
		for fixed_ch in all_channels:

			filename=f'{path}timeresolutionvx_{fixed_ch}.json'
			if not os.path.exists(filename):
				print(f'No file at {filename}')
				continue

			with open(filename, 'r') as x:
				d=json.load(x)

			if histkey not in d:
				print(f'histkey {histkey} for {fixed_ch} not in d')
				continue
			if state not in d[histkey]:
				print(f'state {state} for {fixed_ch} not in d[{histkey}]')
				continue
			if type(d[histkey][state])==list:
				print(f'wtf: {filename}')
			try: res = d[histkey][state][fitmode]
			except KeyError:
				print(f'{fixed_ch} has no {state} time res for mode {fitmode}')
				continue

			self.deltatimeresolutiondict[fixed_ch] = {float(k):v for k,v in res.items()}

	def Makedeltacscintdict(self, runNr):
		path=f'{self.path}cscintvalues/run{runNr}/'
		self.deltacscintdict={}
		all_channels=self.GetListOfChannels(2)
		for fixed_ch in all_channels:

			filename=f'{path}cscint_{fixed_ch}.json'
			if not os.path.exists(filename):
				print(f'No file at {filename}')
				continue

			with open(filename, 'r') as x:
				d=json.load(x)

			if 'corrected' not in d or 'uncorrected' not in d:
				print(f'No corrected or uncorrected state for {fixed_ch}')
				continue

			res = d['corrected'][0]-d['uncorrected'][0]
			d_res = np.sqrt(d['corrected'][-1]**2 + d['uncorrected'][-1]**2)

			self.deltacscintdict[fixed_ch] = res, d_res

	def GetPolyParams(self, fixed_ch, n, runNr='005408'):

		fname=f'{self.path}Polyparams/run{runNr}/polyparams{n}_{fixed_ch}.json'
		if not os.path.exists(fname):
			return

		if os.path.getsize(fname) == 0:
			return

		if n < 7:

			with open(fname, 'r') as f:
				alldata = json.load(f)
				data=alldata['uncorrected']
			if n==5:
				params=[float(i) for i in data[1:13]]
				limits=[float(i) for i in data[13:15]]
				return params, limits
			else: print(f'No tw correction parameters stored for n={n}')

		else:
			with open(fname, 'r') as f:
				alldata = json.load(f)
			params = list(alldata['params'].values())
			return params

	def GetRMSresidual(self, fixed_ch, n, runNr='005408'):

		fname=f'{self.path}Polyparams/run{runNr}/polyparams{n}_{fixed_ch}.json'
		if not os.path.exists(fname):
			return

		if os.path.getsize(fname) == 0:
			return

		if n < 7:

			with open(fname, 'r') as f:
				alldata = json.load(f)
				data=alldata['uncorrected']
			if n==5:
				return float(data[15])
			else: print(f'No tw correction parameters stored for n={n}')

		else:
			with open(fname, 'r') as f:
				alldata = json.load(f)
			return alldata['rms_residual']

	def Get_cscintRMSresidual(self, fixed_ch, state, runNr='005408'):

		fname = f'{self.path}rootfiles/run{runNr}/timewalk_{fixed_ch}.root'
		if not os.path.exists(fname):
			return
		f = ROOT.TFile.Open(fname, "READ")

		try: canv = f.Get(f"linearfitcanvas_{fixed_ch}_{state}")
		except AttributeError:
			print(f'No canvas found for {fixed_ch} {state} in {fname}')
			return

		try: pad = canv.GetPrimitive(f"linearfitresidualspad_{fixed_ch}")
		except AttributeError:
			print(f'No pad found for {fixed_ch} {state} in {fname}')
			return

		try: hist = pad.GetPrimitive(f"twfitresidualshist_{fixed_ch}")
		except AttributeError:
			print(f'No histogram found for {fixed_ch} {state} in {fname}')
			return

		nbins = hist.GetNbinsX()

		# ----- gather bin centres, contents and errors -----------------
		centres = []
		contents = []
		for i in range(1, nbins + 1):       # ROOT bins start at 1
			centres.append(hist.GetBinCenter(i))
			contents.append(hist.GetBinContent(i))

		# Convert to Python lists → easier to work with
		residuals = np.array(contents)

		# ----- calculate RMS ------------------------------------------
		rms = np.sqrt(1/(len(residuals)-1)*np.mean(residuals**2))

		return rms

	def GetErrorFloor(self, fixed_ch, n, runNr='005408'):

		fname=f'{self.path}Polyparams/run{runNr}/polyparams{n}_{fixed_ch}.json'
		if not os.path.exists(fname):
			return

		if os.path.getsize(fname) == 0:
			return

		if n < 7:

			with open(fname, 'r') as f:
				alldata = json.load(f)
				data=alldata['uncorrected']
			if n==5:
				return float(data[15])
			else: print(f'No tw correction parameters stored for n={n}')

		else:
			with open(fname, 'r') as f:
				alldata = json.load(f)
			return alldata['alpha']

	def GetNDF(self, fixed_ch, n, runNr='005408'):

		fname=f'{self.path}Polyparams/run{runNr}/polyparams{n}_{fixed_ch}.json'
		if not os.path.exists(fname):
			return

		if os.path.getsize(fname) == 0:
			return

		if n < 7:

			with open(fname, 'r') as f:
				alldata = json.load(f)
				data=alldata['uncorrected']
			if n==5:
				return float(data[15])
			else: print(f'No tw correction parameters stored for n={n}')

		else:
			with open(fname, 'r') as f:
				alldata = json.load(f)
			return alldata['ndf']

	def Gettds0relativetime(self, runNr, fixed_ch, mode='mean', state='uncorrected', n=5):

		iteration=0 if state=='uncorrected' else 1
		if not os.path.exists(f'{self.path}Polyparams/run{runNr}/polyparams{n}_{fixed_ch}.csv'): return
		with open(f'{self.path}Polyparams/run{runNr}/polyparams{n}_{fixed_ch}.csv', 'r') as f:
			reader=csv.reader(f)
			alldata=[r for r in reader]
			if len(alldata)==0:return
			data=alldata[iteration]

		if n==4:
			if len(data)!=15:
				print(f'No alignment constant saved for {fixed_ch}')
				return
			else: tds0tSiPMmean=[float(i) for i in data[13:]]
		elif n==5:
			if len(data)!=17:
				print(f'No alignment constant saved for {fixed_ch}')
				return
			else: tds0tSiPMmean=[float(i) for i in data[15:]]
		else: print(f'No tw correction parameters stored for n={n}')
		return tds0tSiPMmean

	def Gettimeresolution(self, runNr, fixed_ch, state):

		filename=f'{self.path}TimeResolution/run{runNr}/timeresolution_{fixed_ch}.json'

		if not os.path.exists(filename): return

		with open(filename, 'r') as x:
			d=json.load(x)

		if state not in d: return

		return d[state][0], d[state][1]

	######### THESE ARE HARDCODED TO JUST GET THE MODE THAT IS USED IN THE ANALYSIS ##########
	def GetBarsideTimeresolution(self, runNr, state, mode='FWHM'):

		fname = f'{self.path}TimeResolution/run{runNr}/barside-timeresolutions.json'
		with open(fname, 'r') as f:
			res = json.load(f)
		return res

	def GetBarTimeresolution(self, runNr, state, mode='FWHM'):

		fname = f'{self.path}TimeResolution/run{runNr}/bar-timeresolutions.json'
		with open(fname, 'r') as f:
			res = json.load(f)
		return res
	###########################################################################################

	def CalculateBarsideTimeresolution(self, runNr, detID, side, state='corrected'):

		if not hasattr(self, "timeresolutiondict"): self.Maketimeresolutiondict(runNr, state)

		covar_component = self.GetCovariances(detID, side)
		SiPMresolution_component,N = self.GetResolutions(detID, side)

		calc_barside_resolution = np.sqrt((1/N**2)*(SiPMresolution_component + covar_component))
		return calc_barside_resolution

	def GetResolutions(self, detID, side, state='corrected'):

		if side=='left':SiPMs = [0,1,3,4,6,7]
		else: SiPMs = [8,9,11,12,14,15]

		resolutions = []
		for SiPM in SiPMs:
			key = f'{detID}_{SiPM}'
			if not key in self.timeresolutiondict: continue
			resolution = self.timeresolutiondict[key][0]

			resolutions.append(resolution)
		SiPMresolution_component = sum([i**2 for i in resolutions])
		return SiPMresolution_component, len(SiPMs)

	def GetCovariances(self, detID, side):

		if not hasattr(self, "timingcovariance"):
			self.MakeTimingCovarianceDict()

		if side=='left':SiPMs = [0,1,3,4,6,7]
		else: SiPMs = [8,9,11,12,14,15]
		combs = list(combinations(SiPMs, 2))

		covars=[]
		for i,j in combs:

			key = f'timingxt_{detID}_SiPMs{i}-{j}'
			if not key in self.timingcovariance: continue
			covariance = self.timingcovariance[key]
			covars.append(covariance)
		final_component = sum([2*i for i in covars])
		return final_component

	def GetSkewness(self, hist, xlow, xhigh):
		"""
		Calculate the skewness of a ROOT histogram within a specified range.

		Parameters:
			hist (TH1): The ROOT histogram.
			range_min (float, optional): Lower bound of the range. Use the histogram minimum if None.
			range_max (float, optional): Upper bound of the range. Use the histogram maximum if None.

		Returns:
			float: The skewness of the histogram in the specified range.
		"""
		# Initialize variables
		sum_weights = 0  # Total weight
		sum_x = 0        # First moment (mean numerator)
		sum_x2 = 0       # Second moment
		sum_x3 = 0       # Third moment

		# Loop over bins
		for bin_idx in range(1, hist.GetNbinsX() + 1):
			bin_center = hist.GetBinCenter(bin_idx)
			bin_content = hist.GetBinContent(bin_idx)

			# Only consider bins within the range
			if xlow <= bin_center <= xhigh:
				sum_weights += bin_content
				sum_x += bin_content * bin_center
				sum_x2 += bin_content * bin_center**2
				sum_x3 += bin_content * bin_center**3

		# Calculate mean (mu1)
		if sum_weights == 0:
			raise ValueError("No data in the specified range.")

		mean = sum_x / sum_weights

		# Calculate central moments
		mu2 = (sum_x2 / sum_weights) - mean**2
		mu3 = (sum_x3 / sum_weights) - 3 * mean * mu2 - mean**3

		# Avoid division by zero
		if mu2 <= 0:
			raise ValueError("Variance (mu2) is zero or negative; skewness is undefined.")

		# Calculate skewness
		skewness = mu3 / mu2**1.5
		return skewness

	def FitForMPV(self, runNr, fixed_ch, state):
		fname=f'{self.path}rootfiles/run{runNr}/timewalk_{fixed_ch}.root'
		if not os.path.exists(fname): return -999.
		f=ROOT.TFile.Open(fname, 'READ')
		histname=f'dtvqdc_{fixed_ch}_{state}'
		if not histname in [k.GetName() for k in f.GetListOfKeys()]: return -999
		tmp=f.Get(histname)
		hist=tmp.Clone()
		hist.SetDirectory(ROOT.gROOT)
		f.Close()
		xproj=hist.ProjectionX()
		# xmin,xmax=xproj.GetXaxis().GetXmin(), xproj.GetXaxis().GetXmax()
		mode=xproj.GetBinCenter(xproj.GetMaximumBin())
		res=self.fit_langau(xproj, 'LQ',max(mode-2, 0))
		# 'LQ',0.8*tmp.GetBinCenter(bmin),1.5*tmp.GetBinCenter(bmax)
		MPV=res.Parameter(1)
		MPV_err=res.ParError(1)
		chi2, NDF= res.Chi2(), res.Ndf()
		return (MPV, MPV_err, chi2, NDF)

	def Getchi2pNDF(self, runNr, fixed_ch, state, n=5):

		if n<7:
			fname=f'{self.path}chi2s/run{runNr}/chi2s{n}_{fixed_ch}.json'
			if not os.path.exists(fname):
				print(f'no file {fname}')
				return
			with open(fname) as jf:
				data=json.load(jf)
			chi2pNDF = data[state][0]/data[state][1]

			return chi2pNDF

		else:
			fname=f'{self.path}Polyparams/run{runNr}/polyparams{n}_{fixed_ch}.json'

			if not os.path.exists(fname):
				return

			if os.path.getsize(fname) == 0:
				return

			with open(fname, 'r') as f:
				alldata = json.load(f)
			chi2 = alldata['chi2']
			ndf = alldata['ndf']
			return chi2/ndf

	def GetBadchi2pNDFdict(self, runNr, subsystem, state):

		iteration=0 if state=='uncorrected' else 1
		path=f'{self.path}/chi2s/run{runNr}/'
		res={}
		for filename in os.listdir(path):
			#fixed_ch=filename[filename.find(str(subsystem)):filename.find('.csv')]
			fixed_ch=f'{filename.split("_")[1]}_{filename.split("_")[2].split(".")[0]}'
			if fixed_ch[0]!=str(subsystem): continue
			with open(path+filename, 'r') as handle:
				reader=csv.reader(handle)
				alldata=[row for row in reader]
				if len(alldata)<iteration:
					print(filename)
					continue
				data=alldata[iteration-1]
				#chi2pNDF=data.pop()
				chi2pNDF=float(data[1])/int(data[2])
				res[fixed_ch]=chi2pNDF
		sorted_tuples=sorted(res.items(), key=lambda x:x[1])
		sorted_d={k:v for k,v in sorted_tuples}
		return sorted_d

	def Makechi2pNDFdict(self, runNr, subsystem, state, n=5):

		iteration=0 if state=='uncorrected' else 1
		path=f'{self.path}/chi2s/run{runNr}/'
		res={}
		for filename in os.listdir(path):
			if filename.split('_')[0].find(f'chi2s{n}') ==-1: continue
			#fixed_ch=filename[filename.find(str(subsystem)):filename.find('.csv')]
			fixed_ch=f'{filename.split("_")[1]}_{filename.split("_")[2].split(".")[0]}'
			if fixed_ch[0]!=str(subsystem): continue
			with open(path+filename, 'r') as handle:
				reader=csv.reader(handle)
				alldata=[row for row in reader]
				if len(alldata)==0 or len(alldata)<iteration: return -999.
				if len(alldata)<iteration:
					print(filename)
					continue
			data=alldata[iteration-1]
			chi2pNDF=float(data[1])/int(data[2])
			res[fixed_ch]=chi2pNDF
		sorted_tuples=sorted(res.items(), key=lambda x:x[1])
		sorted_d={k:v for k,v in sorted_tuples}
		return sorted_d

	def GetAttenuationLengthData(self, runNr, fixed_ch):
		fname=self.path+'attenuationlengths/run'+str(runNr)+'/csvfiles/attenuationlength_'+fixed_ch+'.csv'
		if not os.path.exists(fname): return -999.
		with open(fname, 'r') as handle:
			reader=csv.reader(handle)
			alldata=[row for row in reader]
			data=alldata[0]
			return data

	def GetMPV(self, runNr, fixed_ch, iteration):
		fname=f'{self.path}MPVs/run{runNr}/MPV_{fixed_ch}.csv'
		if not os.path.exists(fname): return -999.
		with open(fname, 'r') as h:
			reader=csv.reader(h)
			data=[row for row in reader]
			if data==[]: return -999. # To be investigated
			res=data[iteration-1]
		return float(res[0])

	"""
	Important note: the Analysis.correct_ToF function corrects
	the SiPM time to x=L/2 in the physics FoR. That is not equal to the bar centre!!!!
	"""

	def correct_ToF(self, MuFilter, fixed_ch, clock, xEx):
		detID=int(fixed_ch.split('_')[0])
		SiPM=int(fixed_ch.split('_')[-1])

		# print(f'time: {clock*6.25}, x: {xEx}')

		# Correct to the centre of the bar in the physics FoR
		MuFilter.GetPosition(detID, self.A, self.B)
		xref = 0.5 * (self.A.x() + self.B.x())
		if not self.simulation:
			cs = self.cscintvalues[fixed_ch]
			c_SiPM = float(cs[0])
			time = clock*self.TDC2ns
		# Simulation times are in nanoseconds already
		else:
			c_SiPM = MuFilter.GetConfParF(f'MuFilter/US_signalspeed_{fixed_ch}')
			time=clock

		ToFcorrection=(xEx-xref)/c_SiPM

		if SiPM<8: corrected_t = time + ToFcorrection
		else: corrected_t = time - ToFcorrection

		return (SiPM, corrected_t)

	def GetCorrectedTimes(self, hit, **kwargs):

		x=kwargs.get('x', 0)
		mufilter=kwargs.get('MuFilter', None)
		mode=kwargs.get('mode', 'aligned')

		if mufilter==None:
			print(f'Not passed MuFilter to GetCorrectedTimes!!')
			return

		detID=hit.GetDetectorID()

		alignedtimes=[]
		clocks, qdcs=hit.GetAllTimes(), hit.GetAllSignals()
		for i in clocks:
			SiPM, clock = i
			if self.IsSmallSiPMchannel(SiPM): continue

			fixed_ch=f'{detID}_{SiPM}'

			# Using qdc == -1 as a flag to only correct tof, for simulation data
			if mode == 'tof': qdc=-1
			else: qdc=self.GetChannelVal(SiPM, qdcs)

			if not self.simulation: correctedtime=self.MuFilterCorrectedTime(mufilter, fixed_ch, clock, qdc, x)
			# Apply ToF correction for the simulated data
			else: correctedtime = self.CorrectSimulatedToF(mufilter, fixed_ch, clock, x)

			if correctedtime==None:
				continue

			if mode=='aligned' and not self.simulation:
				d = self.alignmentparameters[f'{detID}_{SiPM}']
				correctedtime = self.task.reft - correctedtime - d[0]
			elif mode=='tof' and not self.simulation:
				correctedtime = self.task.reft - correctedtime

			alignedtimes.append((SiPM, correctedtime))
		return alignedtimes

	def MuFilterCorrectedTime(self, MuFilter, fixed_ch, clock, qdc=-1, x=0):

		time=clock*self.TDC2ns
		if not self.simulation:

			# To correct time: need tw params, alignment param and cscint value if correcting for signal ToF
			if not fixed_ch in self.twparameters:
				print(f'no tw params for {fixed_ch}')
				return
			if (not fixed_ch in self.alignmentparameters):
				if self.options.referencesystem==3 or (self.options.referencesystem==1 and int(fixed_ch[4])<8): print(f'no d param for {fixed_ch}')
				return
			if fixed_ch not in self.cscintvalues and x!=0:
				return

			#### Correct ToF if needed.
			if x==0:
				ToFcorrectedtime = time
			else:
				ToFcorrectedtime = self.correct_ToF(MuFilter, fixed_ch, clock, x)[1]

			# No timewalk correction if -1 passed as qdc (used for simulation data)
			if qdc==-1:
				ToFTWcorrectedtime = ToFcorrectedtime
			else:
				#### TW corrected time then ToF & TW corrected time
				twparams = self.twparameters[fixed_ch]
				twcorrection = self.correctionfunction(twparams, qdc)
				ToFTWcorrectedtime = ToFcorrectedtime+twcorrection

			return ToFTWcorrectedtime

		else:
			#### Correct ToF if needed.
			if x==0:
				ToFcorrectedtime=time
			else:
				# if not hasattr(self.task.M.MuFilter.)
				if self.task.M.MuFilter.GetConfParF(f'MuFilter/US_signalspeed_{fixed_ch}')==0:
					return
				ToFcorrectedtime=self.correct_ToF(self.task.M.MuFilter, fixed_ch, clock, x)[1]
			return ToFcorrectedtime

	def CorrectSimulatedToF(self, MuFilter, fixed_ch, time, xEx):

		detID, SiPM = int(fixed_ch.split('_')[0]), int(fixed_ch.split('_')[-1])
		signalspeed = MuFilter.GetConfParF(f'MuFilter/US_signalspeed_{fixed_ch}')
		if signalspeed == 0:
			print(f'No signalspeed for {fixed_ch}')
			return

		side = 'left' if SiPM < 8 else 'right'

		MuFilter.GetPosition(detID, self.A, self.B)
		x_m = 0.5 * (self.A.x() + self.B.x())
		x = time * signalspeed
		tof_corr = (xEx-x_m)/signalspeed

		if side == 'left':
			# tof_corr = time - dx/signalspeed
			# phys_frame_t = 1/signalspeed * (x + x_m)
			# tof_corr = phys_frame_t + tof_corr
			tof_corr = time - tof_corr

		elif side == 'right':
			# phys_frame_t = 1/signalspeed * (x_m - x)
			# tof_corr = phys_frame_t - tof_corr
			tof_corr = time + tof_corr

		return tof_corr

	def GetQDCpeak(self, runNr, fixed_ch):

		filename=f'{self.path}rootfiles/run{runNr}/timewalk_{fixed_ch}.root'
		if not os.path.exists(filename):return

		f=ROOT.TFile.Open(filename, 'READ')
		histname=f'dtvqdc_{fixed_ch}_corrected'
		if not hasattr(f, histname):
			f.Close()
			return
		hist=f.Get(histname)
		xproj=hist.ProjectionX()
		qdcpeak=xproj.GetBinLowEdge(xproj.GetMaximumBin())
		error=xproj.GetStdDev()/ROOT.TMath.Sqrt(xproj.GetEntries())
		f.Close()
		return qdcpeak, error

	def GetCutDistributions(self, runNr, distmodes=('dy', 'slopes', 'nSiPMs', 'timingdiscriminant'), task='TimeWalk'):
		Allmodes=('dy', 'nSiPMs', 'slopes', 'timingdiscriminant')
		filename=f'{self.path}rootfiles/run{runNr}/SelectionCriteria.root'
		if not os.path.exists(filename):
			if self.timealignment=='old': filename=f'{self.path}rootfiles/run005097/SelectionCriteria.root'
			elif self.timealignment=='new': filename=f'{self.path}rootfiles/run005408/SelectionCriteria.root'
			elif self.timealignment=='new+LHCsynch': filename=f'{self.path}rootfiles/run005999/SelectionCriteria.root'

		if isinstance(distmodes, str):
			distmodes=(distmodes,)

		for distmode in distmodes:
			if distmode not in Allmodes:
				print('Use a valid mode.')
				return -999

		f=ROOT.TFile.Open(filename, 'READ')
		dists={}

		for distmode in distmodes:
			if distmode=='dy' or distmode=='nSiPMs':
				for p in range(self.systemAndPlanes[2]):
					key=str(20+p)
					name=f'{distmode}_{key}'
					if not hasattr(f, name):
						print(f'No hist {name}')
						f.Close()
						return -999
					hist=f.Get(name)

					# Clone so it’s independent of the file
					hist_clone = hist.Clone()
					hist_clone.SetDirectory(0)   # detach from file

					# Rename after clone + detach
					new_name = hist_clone.GetName()
					if task == 'SelectionCriteria':
						hist_clone.SetName(f'sc-{new_name}')

					dists[new_name] = hist_clone

					# hist=f.Get(name).Clone()
					# name=hist.GetName()
					# if task=='SelectionCriteria': hist.SetName(f'sc-{name}')
					# hist.SetDirectory(ROOT.gROOT)
					# dists[name]=hist

			elif distmode=='slopes':
				if not hasattr(f, distmode):
					f.Close()
					print(f'No {distmode} hist')
					return -999
				hist=f.Get(f'{distmode}').Clone()
				hist.SetDirectory(ROOT.gROOT)
				dists[distmode]=hist

			elif distmode=='timingdiscriminant':
				if not hasattr(f, distmode):
					f.Close()
					print(f'No {distmode} hist')
					return -999
				hist=f.Get(f'{distmode}').Clone()
				name=hist.GetName()
				if task=='SelectionCriteria': hist.SetName(f'sc-{name}')
				hist.SetDirectory(ROOT.gROOT)
				dists[distmode]=hist

		f.Close()
		return dists

	def GetTimeAlignmentType(self, runNr):
		if not isinstance(runNr, int): runNr=int(runNr)

		if any( [int(runNr) < 5116, int(runNr) > 5174 and int(runNr) < 5193] ): return 'old'
		elif int(runNr) < 5413: return 'new'
		else: return 'new+LHCsynch'

	def isLHCsynch(self, runNr):
		if int(runNr)<5431: return False
		else: return True

	def GetSiPMtime(self, runNr, fixed_ch, state, mode='mean'):

		filename=f'{self.path}rootfiles/run{runNr}/timewalk_{fixed_ch}.root'
		if not os.path.exists(filename): return

		f=ROOT.TFile.Open(filename, 'READ')
		histname=f'tSiPM_{fixed_ch}_{state}'
		if not hasattr(f, histname):
			f.Close()
			return
		hist=f.Get(histname)
		if mode=='mean':
			tSiPM=hist.GetMean()
			uncertainty=hist.GetStdDev()
		elif mode=='median':
			pass
		f.Close()
		return tSiPM, uncertainty

	def GetAlignmentParameters(self, runNr, fixed_ch):
		fname=f'{self.path}Alignmentparams/run{runNr}/alignmentparams_{self.refsysname}_{fixed_ch}.csv'
		if not os.path.exists(fname):
			return
		with open(fname, 'r') as handle:
			reader=csv.reader(handle)
			alldata=[row for row in reader]
			if len(alldata) == 0: return
			try:
				data=alldata[0]
			except IndexError:
				print(f'{fixed_ch} IndexError')
		return (float(data[0]), float(data[1]))

	def Gettds0mean(self, runNr, fixed_ch, mode='mean', state='aligned'):

		filename=f'{self.path}rootfiles/run{runNr}/timewalk_{fixed_ch}.root'
		if not os.path.exists(filename):
			print(f'No timewalk file for {fixed_ch}')
			return

		f=ROOT.TFile.Open(filename, 'READ')

		if state=='aligned': histname=f'dt_{fixed_ch}_{self.refsysname}aligned'
		else: print(f'No conditions met for finding histname:\n')

		if not hasattr(f, histname):
			f.Close()
			print(f'No hist {histname} for {fixed_ch}')
			return
		hist=f.Get(histname)

		if mode=='mean':
			mean=hist.GetMean()
			uncertainty = mean/ROOT.TMath.Sqrt(hist.GetEntries())
			f.Close()
			return mean, uncertainty

		elif mode=='truncated':
			mode = hist.GetBinCenter(hist.GetMaximumBin())
			timing_range = (mode - 1*hist.GetStdDev(), mode + 1*hist.GetStdDev())
			# print(timing_range)
			# Select bins within the timing range
			bin_centres = [hist.GetBinCenter(i) for i in range(1, hist.GetNbinsX() + 1) if timing_range[0] <= hist.GetBinCenter(i) <= timing_range[1]]
			frequencies = [hist.GetBinContent(i) for i in range(1, hist.GetNbinsX() + 1) if timing_range[0] <= hist.GetBinCenter(i) <= timing_range[1]]

			mean_weighted = np.average(bin_centres, weights=frequencies)
			variance_weighted = np.average((np.array(bin_centres) - mean_weighted)**2, weights=frequencies)
			std_dev_weighted = np.sqrt(variance_weighted)
			return mean_weighted, std_dev_weighted/np.sqrt(sum(frequencies)) # Standard error
			# modaltime=hist.GetBinCenter(hist.GetMaximumBin())
			# low, high=hist.FindBin(modaltime-2*hist.GetStdDev()), hist.FindBin(modaltime+2*hist.GetStdDev()) # mode(hist) + 2*stddev, mode(hist) - 2*std dev
			# tmp=[[hist.GetBinCenter(i),hist.GetBinContent(i)] for i in range(low, high+1)] # Data within 2 standard deviations of the mode

			# bin_centers = np.array([x[0] for x in tmp])
			# frequencies = np.array([x[1] for x in tmp])

			# # Treat as weighted data
			# mean_weighted = np.average(bin_centers, weights=frequencies)
			# variance_weighted = np.average((bin_centers - mean_weighted)**2, weights=frequencies)
			# std_dev_weighted = np.sqrt(variance_weighted)

			# # Effective sample size for weighted data
			# # Using Kish's approximation formula
			# effective_n = np.sum(frequencies)**2 / np.sum(frequencies**2) if np.sum(frequencies**2) > 0 else len(frequencies)

			# standard_error = std_dev_weighted / np.sqrt(effective_n)
			# f.Close()
			# return mean_weighted, standard_error

		elif mode=='mode':
			correction=hist.GetBinCenter(hist.GetMaximumBin())
			uncertainty=hist.GetStdDev()/ROOT.TMath.Sqrt(hist.GetEntries())
			return correction, uncertainty
		f.Close()

	def GettSiPMcorrectedmean(self, runNr, fixed_ch, mode='mean'):
		filename=f'{self.path}rootfiles/run{runNr}/timewalk_{fixed_ch}.root'
		if not os.path.exists(filename): return
		f=ROOT.TFile.Open(filename, 'READ')
		# if self.timealignment=='new': histname=f'tSiPMToFcorrected_{fixed_ch}_corrected_tDS0-tSiPMcorrected'
		histname=f'tSiPMToFcorrected_{fixed_ch}_corrected_tDS0-tSiPMcorrected'
		# elif self.timealignment=='old': histname=f'ScifiAlignedToFcorrectedtSiPM_{fixed_ch}_corrected'
		if not hasattr(f, histname):
			f.Close()
			return
		hist=f.Get(histname)
		if mode=='mean':
			correction=hist.GetBinCenter(hist.GetMaximumBin())
			error=hist.GetStdDev()/2
			# for
			# correction=hist.GetMean()

			# uncertainty=hist.GetStdDev()

		elif mode=='median':
			correction, error =0, 0
		else:
			correction, error =0, 0
		f.Close()
		return correction, error

	def GetEntriesInHist(self, runNr, fixed_ch, mode, state):

		filename=f'{self.path}rootfiles/run{runNr}/timewalk_{fixed_ch}.root'
		if not os.path.exists(filename):return

		f=ROOT.TFile.Open(filename, 'READ')
		histname=f'{mode}_{fixed_ch}_{state}'
		if not hasattr(f, histname):
			f.Close()
			return
		hist=f.Get(histname)
		entries=hist.GetEntries()
		f.Close()
		return entries

	def MakeTWCorrectionDict(self, n, withErrors=False):

		d={}
		for s in (2,): # Only make dict for US
			for p in range(self.systemAndPlanes[s]):
				for b in range(self.systemAndBars[s]):
					for SiPM in self.systemAndSiPMs[s]:
						fixed_ch=self.MakeFixedCh((s,p,b,SiPM))
						tmp=self.GetPolyParams(fixed_ch, n, '005999')
						if not tmp:
							print(f'No tw params for {fixed_ch}')
							continue

						if n<7:
							paramsAndErrors=tmp[0]
							if not withErrors: params=self.correctionparams(paramsAndErrors)
							else: params=paramsAndErrors

						else:
							params = tmp
						d[fixed_ch]=params
		if len(d)==0: self.twparameters=None
		self.twparameters=d

	def WriteTWParamDict(self):
		if not hasattr(self, "twparameters"): self.MakeTWCorrectionDict(self.timealignment)

		twparamsdir = f'{self.path}/Polyparams/run{self.runNr}/'
		twparamsfilename=twparamsdir+f'twparams.json'
		with open(twparamsfilename, 'w') as jf:
			json.dump(self.twparameters, jf)
		print(f'TW param dict written to {twparamsfilename}')

	### Make dictionary of the alignment parameter determined as the truncated y-mean of tw-corr (tds0 - tSiPM)
	def MakeAlignmentParameterDict(self):

		print(f'Making alignment parameter dict for {self.timealignment} alignment.')

		if self.timealignment=='old': run=str(5097).zfill(6)
		elif self.timealignment=='new': run=str(5408).zfill(6)
		elif self.timealignment=='new+LHCsynch': run=str(5999).zfill(6)

		d={}
		for s in (2,): # Just US
			for p in range(self.systemAndPlanes[s]):
				for b in range(self.systemAndBars[s]):
					for SiPM in self.systemAndSiPMs[s]:
						fixed_ch=self.MakeFixedCh((s,p,b,SiPM))
						correction=self.GetAlignmentParameters(run, fixed_ch)
						if not correction: continue
						d[fixed_ch]=correction
		if len(d)==0: self.alignmentparameters=None
		self.alignmentparameters=d

	def WriteAlignmentParamDict(self):
		if not hasattr(self, "alignmentparameters"): self.MakeAlignmentParameterDict()

		afs_alignmentparamsfilename=f'{self.path}/Alignmentparams/run{self.runNr}/alignmentparameter{self.refsysname}.json'
		eos_alignmentparamsfilename=f'/eos/home-a/aconsnd/SWAN_projects/Data analysis/data/run{self.runNr}_alignmentparameter{self.refsysname}.json'
		for alignmentparamsfilename in (afs_alignmentparamsfilename, eos_alignmentparamsfilename):
			with open(alignmentparamsfilename, 'w') as jf:
				json.dump(self.alignmentparameters, jf)
			print(f'Alignment param dict written to {alignmentparamsfilename}')

	def MakeTimingCovarianceDict(self, mode='truncated'):

		if self.timealignment=='old': run=str(5097).zfill(6)
		elif self.timealignment=='new': run=str(5408).zfill(6)
		elif self.timealignment=='new+LHCsynch': run=str(5999).zfill(6)

		if mode=='truncated': filename=f'{self.path}TimingCovariance/run{run}/run{run}_truncatedcovariance.json'
		elif mode=='full': filename=f'{self.path}TimingCovariance/run{run}/timingcovariance.json'
		else: return

		with open(filename, 'r') as x:
			self.timingcovariance = json.load(x)

	def SaveTimingCovarianceDict(self, mode='truncated'):

		if self.timealignment=='old': run=str(5097).zfill(6)
		elif self.timealignment=='new': run=str(5408).zfill(6)
		elif self.timealignment=='new+LHCsynch': run=str(5999).zfill(6)

		if mode=='truncated': filename=f'{self.path}TimingCovariance/run{run}/run{run}_truncatedcovariance.json'
		elif mode=='full': filename=f'{self.path}TimingCovariance/run{run}/timingcovariance.json'
		else: return

		if not hasattr(self, 'timingcovariance'): self.MakeTimingCovarianceDict(mode)

		afs_filename=f'{self.path}TimingCovariance/run{self.runNr}/run{self.runNr}_{mode}covariance.json'
		eos_filename=f'/eos/home-a/aconsnd/SWAN_projects/Data analysis/data/run{self.runNr}_{mode}covariance.json'

		for filename in (afs_filename, eos_filename):
			with open(filename, 'w') as jf:
				json.dump(self.timingcovariance, jf)
			print(f'Timing covariance dict written to {filename}')

	def MakeTimingCorrelationDict(self, runNr):
		filename=f'{self.path}TimingCovariance/run{self.runNr}/timingcorrelation.json'
		with open(filename) as jsonfile:
			self.timingcorrelation=json.load(jsonfile)

	def MakeQDCMIPJson(self, runNr):
		### Currently
		self.langaufun()
		d={}
		for s in (1,2):
			for p in range(self.systemAndPlanes[s]):
				for b in range(self.systemAndBars[s]):
					print(f'{s}, {p}, {b}')
					for SiPM in self.systemAndSiPMs[s]:
						fixed_ch=self.MakeFixedCh((s,p,b,SiPM))
						MPV, MPV_err, chi2, NDF = self.FitForMPV(runNr, fixed_ch, 'corrected')

						d[fixed_ch] = MPV, MPV_err

		mpvpath=f'{self.afswork}MPVs/run{self.runNr}/'
		mpvfilename=mpvpath+f'MPVs.json'
		if not os.path.exists(mpvfilename):
			os.makedirs(mpvpath, exist_ok=True)
		with open(mpvfilename, 'w') as outfile:
			json.dump(d, outfile, indent=4)
		print(f'mpv dictionary written to {mpvfilename}')

	def GetAverageSiPMTimingCovariance(self, runNr):
		r={}

		if not hasattr(self, 'timingcovariance'): self.MakeTimingCovarianceDict()

		for key in self.timingcovariance:
			if key=='timingxt_reft_US': continue
			x, detID, tmp=key.split('_')

			SiPMs=tmp[len('SiPMs'):].split('-')

			#### Special case for last SiPM on bar end due to how combinations algorithm works
			if SiPMs[1] in ('7', '15'):
				fixed_ch=f'{detID}_{SiPMs[1]}'
				if not fixed_ch in r: r[fixed_ch]=[]
				r[fixed_ch].append(self.timingcovariance[key])

			fixed_ch = f'{detID}_{SiPMs[0]}'
			if not fixed_ch in r: r[fixed_ch]=[]
			r[fixed_ch].append(self.timingcovariance[key])

		for fixed_ch in r:
			avg=1/len(r[fixed_ch]) * sum(r[fixed_ch])
			r[fixed_ch]=avg

		self.SiPMaveragetimingcovariance=r

	def GetCovariance(self, runNr, detID, SiPMs):
		if not hasattr(self, 'timingcovariance'):
			self.MakeTimingCovarianceDict()

		key=f'timingxt_{detID}_SiPMs{SiPMs[0]}-{SiPMs[1]}'
		if not key in self.timingcovariance: return
		else: return self.timingcovariance[key]

	def MakePlaneWiseBarMapThing(self, MuFilter):

		d={}
		for plane in range(5):
			for bar in range(10):
				detID = int(f'2{plane}00{bar}')
				MuFilter.GetPosition(detID, self.A, self.B)
				d[detID] = [self.A.x(), self.B.x(), self.A.y(), self.B.y()]
		for plane in range(4):
			if plane!=3:bars=range(120)
			else: bars=range(60,120)
			for bar in bars:
				bar=str(bar).zfill(3)
				detID=int(f'3{plane}{bar}')
				MuFilter.GetPosition(detID, self.A, self.B)
				d[detID] = [self.A.x(), self.B.x(), self.A.y(), self.B.y()]

		for f in (f'{self.path}BarPositions.json', f'/eos/user/a/aconsnd/SWAN_projects/numuInvestigation/data/BarPositions.json', f'/eos/user/a/aconsnd/SWAN_projects/Simulation/data/BarPositions.json'):
			with open(f, 'w') as jf:
				json.dump(d, jf)
		print(f'Bar positions written to three locations')
		return d

	def GetMCPointBarycentres(self, hit, MuFiPointsInHit):

		detID = hit.GetDetectorID()

		# Get end coordinates of bar
		self.task.MuFilter.GetPosition(detID, self.A, self.B)
		left_x, right_x = self.A[0], self.B[0]

		# Get bar side average signal speeds and time resolutions
		left_signalspeed = self.task.MuFilter.GetBarSideSignalSpeed(detID, "left")
		right_signalspeed = self.task.MuFilter.GetBarSideSignalSpeed(detID, "right")

		left_timeresolution = self.task.MuFilter.GetBarSideTimeResolution(detID, "left")
		right_timeresolution = self.task.MuFilter.GetBarSideTimeResolution(detID, "right")

		# Just doing x-barycentre with MCPoints
		detID = hit.GetDetectorID()

		# Get inter-quartile range of times
		times = np.array([mcp.GetTime() for mcp in MuFiPointsInHit])
		q1 = np.percentile(times, 25)
		q3 = np.percentile(times, 75)
		iqr_range = (q3-q1)

		# Get minimum MC point time
		min_time = min(mcp.GetTime() for mcp in MuFiPointsInHit)

		hit_z = self.A.z() # Assuming z-coordinate is the same for both ends of the bar

		time2leftside = lambda mcp_time, mcp_position : np.random.normal(mcp_time + (self.A-mcp_position).Mag()/left_signalspeed, left_timeresolution)
		time2rightside = lambda mcp_time, mcp_position : np.random.normal(mcp_time + (self.B-mcp_position).Mag()/right_signalspeed, right_timeresolution)

		# dictionaries of mc points and the time of arrival of the associated photons
		times_to_leftside = {idx:time2leftside(mcp.GetTime(), ROOT.TVector3(mcp.GetX(), mcp.GetY(), mcp.GetZ())) for idx, mcp in enumerate(MuFiPointsInHit)}
		times_to_rightside = {idx:time2rightside(mcp.GetTime(), ROOT.TVector3(mcp.GetX(), mcp.GetY(), mcp.GetZ())) for idx, mcp in enumerate(MuFiPointsInHit)}

		min_lefttime_idx = min(times_to_leftside, key=times_to_leftside.get)
		min_righttime_idx = min(times_to_rightside, key=times_to_rightside.get)

		# Take smeared time and remove the time of light creation
		smeared_timeL = times_to_leftside[min_lefttime_idx] - MuFiPointsInHit[min_lefttime_idx].GetTime()
		xL_reco = -smeared_timeL*left_signalspeed + left_x

		# Take smeared time and remove the time of light creation
		smeared_timeR = times_to_rightside[min_righttime_idx] - MuFiPointsInHit[min_righttime_idx].GetTime()
		xR_reco = smeared_timeR*right_signalspeed + right_x

		# Get total energy of the MC points in the hit
		total_energy = sum(mcp.GetEnergyLoss() for mcp in MuFiPointsInHit)

		return {'xL':xL_reco, 'xR':xR_reco, 'tL':smeared_timeL, 'tR':smeared_timeR, 'iqr_range':iqr_range, 'hit_z':hit_z, 'min_time':min_time, 'n_points':len(MuFiPointsInHit), 'total_energy':total_energy}


	# def MCP_times2leftside(mcp, left_signalspeed, left_timeresolution, left_pos):
	# 	mcp_time = mcp.GetTime()
	# 	mcp_position = ROOT.TVectmcp.

	# def MCP_times2rightside(mcp):

	# 	time2leftside = lambda mcp_time, mcp_position : np.random.normal(mcp_time + (A-mcp_position).Mag()/left_signalspeed, left_timeresolution)
	# 	time2rightside = lambda mcp_time, mcp_position : np.random.normal(mcp_time + (B-mcp_position).Mag()/right_signalspeed, right_timeresolution)
