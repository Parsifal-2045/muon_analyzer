# Associators
Some changed from the ipython dump to the running code highlited with comments:
- Purity and efficiency cuts should be "safe" to change, meaning that the new values propagate correctly
- Every Tracks tag gain the additional input tag Pixel
- useGEMs and usePhase2Tracker change from False to True 

## <span style="color:red">L3 from L1tk</span>
```python
cms.EDProducer("MuonAssociatorEDProducer",
    AbsoluteNumberOfHits_muon = cms.bool(False),
    AbsoluteNumberOfHits_track = cms.bool(False),
    CSClinksTag = cms.InputTag("simMuonCSCDigis","MuonCSCStripDigiSimLinks"),
    CSCsimHitsTag = cms.InputTag("g4SimHits","MuonCSCHits"),
    CSCsimHitsXFTag = cms.InputTag("mix","g4SimHitsMuonCSCHits"),
    CSCwireLinksTag = cms.InputTag("simMuonCSCDigis","MuonCSCWireDigiSimLinks"),
    DTdigiTag = cms.InputTag("simMuonDTDigis"),
    DTdigisimlinkTag = cms.InputTag("simMuonDTDigis"),
    DTrechitTag = cms.InputTag("hltDt1DRecHits"),
    DTsimhitsTag = cms.InputTag("g4SimHits","MuonDTHits"),
    DTsimhitsXFTag = cms.InputTag("mix","g4SimHitsMuonDTHits"),
    EfficiencyCut_muon = cms.double(0.0),
    EfficiencyCut_track = cms.double(0.0),
    GEMdigisimlinkTag = cms.InputTag("simMuonGEMDigis","GEM"),
    GEMsimhitsTag = cms.InputTag("g4SimHits","MuonGEMHits"),
    GEMsimhitsXFTag = cms.InputTag("mix","g4SimHitsMuonGEMHits"),
    NHitCut_muon = cms.uint32(0),
    NHitCut_track = cms.uint32(0),
    PurityCut_muon = cms.double(0.25), 
    PurityCut_track = cms.double(0.25), # purity cut untouched 
    ROUList = cms.vstring(
        'TrackerHitsTIBLowTof',
        'TrackerHitsTIBHighTof',
        'TrackerHitsTIDLowTof',
        'TrackerHitsTIDHighTof',
        'TrackerHitsTOBLowTof',
        'TrackerHitsTOBHighTof',
        'TrackerHitsTECLowTof',
        'TrackerHitsTECHighTof',
        'TrackerHitsPixelBarrelLowTof',
        'TrackerHitsPixelBarrelHighTof',
        'TrackerHitsPixelEndcapLowTof',
        'TrackerHitsPixelEndcapHighTof'
    ),
    RPCdigisimlinkTag = cms.InputTag("simMuonRPCDigis","RPCDigiSimLink"),
    RPCsimhitsTag = cms.InputTag("g4SimHits","MuonRPCHits"),
    RPCsimhitsXFTag = cms.InputTag("mix","g4SimHitsMuonRPCHits"),
    ThreeHitTracksAreSpecial = cms.bool(False),
    UseGrouped = cms.bool(True),
    UseMuon = cms.bool(False),
    UsePixels = cms.bool(True),
    UseSplitting = cms.bool(True),
    UseTracker = cms.bool(True),
    acceptOneStubMatchings = cms.bool(False),
    associatePixel = cms.bool(True),
    associateRecoTracks = cms.bool(True),
    associateStrip = cms.bool(True),
    associatorByWire = cms.bool(False),
    crossingframe = cms.bool(False),
    dumpDT = cms.bool(False),
    dumpInputCollections = cms.untracked.bool(False),
    ignoreMissingTrackCollection = cms.untracked.bool(True),
    includeZeroHitMuons = cms.bool(True),
    inputCSCSegmentCollection = cms.InputTag("cscSegments"),
    inputDTRecSegment4DCollection = cms.InputTag("dt4DSegments"),
    links_exist = cms.bool(True),
    phase2TrackerSimLinkSrc = cms.InputTag("simSiPixelDigis","Tracker"),
    pixelSimLinkSrc = cms.InputTag("simSiPixelDigis","Pixel"), # pixel added
    rejectBadGlobal = cms.bool(True),
    simtracksTag = cms.InputTag("g4SimHits"),
    simtracksXFTag = cms.InputTag("mix","g4SimHits"),
    stripSimLinkSrc = cms.InputTag("simSiStripDigis"),
    tpRefVector = cms.bool(True),
    tpTag = cms.InputTag("TPmu"),
    tracksTag = cms.InputTag("hltIter2Phase2L3FromL1TkMuonMerged"),
    useGEMs = cms.bool(True), # changed from False to True
    usePhase2Tracker = cms.bool(True) # changed from False to True
)
```
##  <span style="color:purple">L2 muons from L1 tk</span>
```python
cms.EDProducer("MuonAssociatorEDProducer",
    AbsoluteNumberOfHits_muon = cms.bool(False),
    AbsoluteNumberOfHits_track = cms.bool(False),
    CSClinksTag = cms.InputTag("simMuonCSCDigis","MuonCSCStripDigiSimLinks"),
    CSCsimHitsTag = cms.InputTag("g4SimHits","MuonCSCHits"),
    CSCsimHitsXFTag = cms.InputTag("mix","g4SimHitsMuonCSCHits"),
    CSCwireLinksTag = cms.InputTag("simMuonCSCDigis","MuonCSCWireDigiSimLinks"),
    DTdigiTag = cms.InputTag("simMuonDTDigis"),
    DTdigisimlinkTag = cms.InputTag("simMuonDTDigis"),
    DTrechitTag = cms.InputTag("hltDt1DRecHits"),
    DTsimhitsTag = cms.InputTag("g4SimHits","MuonDTHits"),
    DTsimhitsXFTag = cms.InputTag("mix","g4SimHitsMuonDTHits"),
    EfficiencyCut_muon = cms.double(0.0),
    EfficiencyCut_track = cms.double(0.0),
    GEMdigisimlinkTag = cms.InputTag("simMuonGEMDigis","GEM"),
    GEMsimhitsTag = cms.InputTag("g4SimHits","MuonGEMHits"),
    GEMsimhitsXFTag = cms.InputTag("mix","g4SimHitsMuonGEMHits"),
    NHitCut_muon = cms.uint32(0),
    NHitCut_track = cms.uint32(0),
    PurityCut_muon = cms.double(0.25),
    PurityCut_track = cms.double(0.25), # purity cut untouched
    ROUList = cms.vstring(
        'TrackerHitsTIBLowTof',
        'TrackerHitsTIBHighTof',
        'TrackerHitsTIDLowTof',
        'TrackerHitsTIDHighTof',
        'TrackerHitsTOBLowTof',
        'TrackerHitsTOBHighTof',
        'TrackerHitsTECLowTof',
        'TrackerHitsTECHighTof',
        'TrackerHitsPixelBarrelLowTof',
        'TrackerHitsPixelBarrelHighTof',
        'TrackerHitsPixelEndcapLowTof',
        'TrackerHitsPixelEndcapHighTof'
    ),
    RPCdigisimlinkTag = cms.InputTag("simMuonRPCDigis","RPCDigiSimLink"),
    RPCsimhitsTag = cms.InputTag("g4SimHits","MuonRPCHits"),
    RPCsimhitsXFTag = cms.InputTag("mix","g4SimHitsMuonRPCHits"),
    ThreeHitTracksAreSpecial = cms.bool(False),
    UseGrouped = cms.bool(True),
    UseMuon = cms.bool(True),
    UsePixels = cms.bool(True),
    UseSplitting = cms.bool(True),
    UseTracker = cms.bool(False),
    acceptOneStubMatchings = cms.bool(False),
    associatePixel = cms.bool(True),
    associateRecoTracks = cms.bool(True),
    associateStrip = cms.bool(True),
    associatorByWire = cms.bool(False),
    crossingframe = cms.bool(False),
    dumpDT = cms.bool(False),
    dumpInputCollections = cms.untracked.bool(False),
    ignoreMissingTrackCollection = cms.untracked.bool(True),
    includeZeroHitMuons = cms.bool(True),
    inputCSCSegmentCollection = cms.InputTag("cscSegments"),
    inputDTRecSegment4DCollection = cms.InputTag("dt4DSegments"),
    links_exist = cms.bool(True),
    phase2TrackerSimLinkSrc = cms.InputTag("simSiPixelDigis","Tracker"),
    pixelSimLinkSrc = cms.InputTag("simSiPixelDigis","Pixel"), # pixel added
    rejectBadGlobal = cms.bool(True),
    simtracksTag = cms.InputTag("g4SimHits"),
    simtracksXFTag = cms.InputTag("mix","g4SimHits"),
    stripSimLinkSrc = cms.InputTag("simSiStripDigis"),
    tpRefVector = cms.bool(True),
    tpTag = cms.InputTag("TPmu"),
    tracksTag = cms.InputTag("hltL2MuonsFromL1TkMuon","UpdatedAtVtx"),
    useGEMs = cms.bool(True), # changed from False to True
    usePhase2Tracker = cms.bool(True) # changed from False to True
)
```
## <span style="color:lightblue">L3 muons OI</span>
```python
cms.EDProducer("MuonAssociatorEDProducer",
    AbsoluteNumberOfHits_muon = cms.bool(False),
    AbsoluteNumberOfHits_track = cms.bool(False),
    CSClinksTag = cms.InputTag("simMuonCSCDigis","MuonCSCStripDigiSimLinks"),
    CSCsimHitsTag = cms.InputTag("g4SimHits","MuonCSCHits"),
    CSCsimHitsXFTag = cms.InputTag("mix","g4SimHitsMuonCSCHits"),
    CSCwireLinksTag = cms.InputTag("simMuonCSCDigis","MuonCSCWireDigiSimLinks"),
    DTdigiTag = cms.InputTag("simMuonDTDigis"),
    DTdigisimlinkTag = cms.InputTag("simMuonDTDigis"),
    DTrechitTag = cms.InputTag("hltDt1DRecHits"),
    DTsimhitsTag = cms.InputTag("g4SimHits","MuonDTHits"),
    DTsimhitsXFTag = cms.InputTag("mix","g4SimHitsMuonDTHits"),
    EfficiencyCut_muon = cms.double(0.0),
    EfficiencyCut_track = cms.double(0.0),
    GEMdigisimlinkTag = cms.InputTag("simMuonGEMDigis","GEM"),
    GEMsimhitsTag = cms.InputTag("g4SimHits","MuonGEMHits"),
    GEMsimhitsXFTag = cms.InputTag("mix","g4SimHitsMuonGEMHits"),
    NHitCut_muon = cms.uint32(0),
    NHitCut_track = cms.uint32(0),
    PurityCut_muon = cms.double(0.25),
    PurityCut_track = cms.double(0.25), # purity cut untouched
    ROUList = cms.vstring(
        'TrackerHitsTIBLowTof',
        'TrackerHitsTIBHighTof',
        'TrackerHitsTIDLowTof',
        'TrackerHitsTIDHighTof',
        'TrackerHitsTOBLowTof',
        'TrackerHitsTOBHighTof',
        'TrackerHitsTECLowTof',
        'TrackerHitsTECHighTof',
        'TrackerHitsPixelBarrelLowTof',
        'TrackerHitsPixelBarrelHighTof',
        'TrackerHitsPixelEndcapLowTof',
        'TrackerHitsPixelEndcapHighTof'
    ),
    RPCdigisimlinkTag = cms.InputTag("simMuonRPCDigis","RPCDigiSimLink"),
    RPCsimhitsTag = cms.InputTag("g4SimHits","MuonRPCHits"),
    RPCsimhitsXFTag = cms.InputTag("mix","g4SimHitsMuonRPCHits"),
    ThreeHitTracksAreSpecial = cms.bool(False),
    UseGrouped = cms.bool(True),
    UseMuon = cms.bool(True),
    UsePixels = cms.bool(True),
    UseSplitting = cms.bool(True),
    UseTracker = cms.bool(False),
    acceptOneStubMatchings = cms.bool(False),
    associatePixel = cms.bool(True),
    associateRecoTracks = cms.bool(True),
    associateStrip = cms.bool(True),
    associatorByWire = cms.bool(False),
    crossingframe = cms.bool(False),
    dumpDT = cms.bool(False),
    dumpInputCollections = cms.untracked.bool(False),
    ignoreMissingTrackCollection = cms.untracked.bool(True),
    includeZeroHitMuons = cms.bool(True),
    inputCSCSegmentCollection = cms.InputTag("cscSegments"),
    inputDTRecSegment4DCollection = cms.InputTag("dt4DSegments"),
    links_exist = cms.bool(True),
    phase2TrackerSimLinkSrc = cms.InputTag("simSiPixelDigis","Tracker"),
    pixelSimLinkSrc = cms.InputTag("simSiPixelDigis","Pixel"), # Pixel added
    rejectBadGlobal = cms.bool(True),
    simtracksTag = cms.InputTag("g4SimHits"),
    simtracksXFTag = cms.InputTag("mix","g4SimHits"),
    stripSimLinkSrc = cms.InputTag("simSiStripDigis"),
    tpRefVector = cms.bool(True),
    tpTag = cms.InputTag("TPmu"),
    tracksTag = cms.InputTag("hltL3MuonsPhase2L3OI"),
    useGEMs = cms.bool(True), # changed from False to True
    usePhase2Tracker = cms.bool(True) # changed from False to True
)
```
## <span style="color:Gold">L3 muon merged</span>
```python
cms.EDProducer("MuonAssociatorEDProducer",
    AbsoluteNumberOfHits_muon = cms.bool(False),
    AbsoluteNumberOfHits_track = cms.bool(False),
    CSClinksTag = cms.InputTag("simMuonCSCDigis","MuonCSCStripDigiSimLinks"),
    CSCsimHitsTag = cms.InputTag("g4SimHits","MuonCSCHits"),
    CSCsimHitsXFTag = cms.InputTag("mix","g4SimHitsMuonCSCHits"),
    CSCwireLinksTag = cms.InputTag("simMuonCSCDigis","MuonCSCWireDigiSimLinks"),
    DTdigiTag = cms.InputTag("simMuonDTDigis"),
    DTdigisimlinkTag = cms.InputTag("simMuonDTDigis"),
    DTrechitTag = cms.InputTag("hltDt1DRecHits"),
    DTsimhitsTag = cms.InputTag("g4SimHits","MuonDTHits"),
    DTsimhitsXFTag = cms.InputTag("mix","g4SimHitsMuonDTHits"),
    EfficiencyCut_muon = cms.double(0.0),
    EfficiencyCut_track = cms.double(0.0),
    GEMdigisimlinkTag = cms.InputTag("simMuonGEMDigis","GEM"),
    GEMsimhitsTag = cms.InputTag("g4SimHits","MuonGEMHits"),
    GEMsimhitsXFTag = cms.InputTag("mix","g4SimHitsMuonGEMHits"),
    NHitCut_muon = cms.uint32(0),
    NHitCut_track = cms.uint32(0),
    PurityCut_muon = cms.double(0.25), 
    PurityCut_track = cms.double(0.25), # Purity cuts untouched 
    ROUList = cms.vstring(
        'TrackerHitsTIBLowTof',
        'TrackerHitsTIBHighTof',
        'TrackerHitsTIDLowTof',
        'TrackerHitsTIDHighTof',
        'TrackerHitsTOBLowTof',
        'TrackerHitsTOBHighTof',
        'TrackerHitsTECLowTof',
        'TrackerHitsTECHighTof',
        'TrackerHitsPixelBarrelLowTof',
        'TrackerHitsPixelBarrelHighTof',
        'TrackerHitsPixelEndcapLowTof',
        'TrackerHitsPixelEndcapHighTof'
    ),
    RPCdigisimlinkTag = cms.InputTag("simMuonRPCDigis","RPCDigiSimLink"),
    RPCsimhitsTag = cms.InputTag("g4SimHits","MuonRPCHits"),
    RPCsimhitsXFTag = cms.InputTag("mix","g4SimHitsMuonRPCHits"),
    ThreeHitTracksAreSpecial = cms.bool(False),
    UseGrouped = cms.bool(True),
    UseMuon = cms.bool(False),
    UsePixels = cms.bool(True),
    UseSplitting = cms.bool(True),
    UseTracker = cms.bool(True),
    acceptOneStubMatchings = cms.bool(False),
    associatePixel = cms.bool(True),
    associateRecoTracks = cms.bool(True),
    associateStrip = cms.bool(True),
    associatorByWire = cms.bool(False),
    crossingframe = cms.bool(False),
    dumpDT = cms.bool(False),
    dumpInputCollections = cms.untracked.bool(False),
    ignoreMissingTrackCollection = cms.untracked.bool(True),
    includeZeroHitMuons = cms.bool(True),
    inputCSCSegmentCollection = cms.InputTag("cscSegments"),
    inputDTRecSegment4DCollection = cms.InputTag("dt4DSegments"),
    links_exist = cms.bool(True),
    phase2TrackerSimLinkSrc = cms.InputTag("simSiPixelDigis","Tracker"),
    pixelSimLinkSrc = cms.InputTag("simSiPixelDigis","Pixel"), # Pixel added
    rejectBadGlobal = cms.bool(True),
    simtracksTag = cms.InputTag("g4SimHits"),
    simtracksXFTag = cms.InputTag("mix","g4SimHits"),
    stripSimLinkSrc = cms.InputTag("simSiStripDigis"),
    tpRefVector = cms.bool(True),
    tpTag = cms.InputTag("TPmu"),
    tracksTag = cms.InputTag("hltPhase2L3MuonMerged"),
    useGEMs = cms.bool(True), # changed from False to True
    usePhase2Tracker = cms.bool(True) # changed from False to True
)
```
## <span style="color:ForestGreen">L3 Muon Tracks</span>
```python
cms.EDProducer("MuonAssociatorEDProducer",
    AbsoluteNumberOfHits_muon = cms.bool(False),
    AbsoluteNumberOfHits_track = cms.bool(False),
    CSClinksTag = cms.InputTag("simMuonCSCDigis","MuonCSCStripDigiSimLinks"),
    CSCsimHitsTag = cms.InputTag("g4SimHits","MuonCSCHits"),
    CSCsimHitsXFTag = cms.InputTag("mix","g4SimHitsMuonCSCHits"),
    CSCwireLinksTag = cms.InputTag("simMuonCSCDigis","MuonCSCWireDigiSimLinks"),
    DTdigiTag = cms.InputTag("simMuonDTDigis"),
    DTdigisimlinkTag = cms.InputTag("simMuonDTDigis"),
    DTrechitTag = cms.InputTag("hltDt1DRecHits"),
    DTsimhitsTag = cms.InputTag("g4SimHits","MuonDTHits"),
    DTsimhitsXFTag = cms.InputTag("mix","g4SimHitsMuonDTHits"),
    EfficiencyCut_muon = cms.double(0.0),
    EfficiencyCut_track = cms.double(0.0),
    GEMdigisimlinkTag = cms.InputTag("simMuonGEMDigis","GEM"),
    GEMsimhitsTag = cms.InputTag("g4SimHits","MuonGEMHits"),
    GEMsimhitsXFTag = cms.InputTag("mix","g4SimHitsMuonGEMHits"),
    NHitCut_muon = cms.uint32(0),
    NHitCut_track = cms.uint32(0),
    PurityCut_muon = cms.double(0.25),
    PurityCut_track = cms.double(0.25), # Purity cuts untouched
    ROUList = cms.vstring(
        'TrackerHitsTIBLowTof',
        'TrackerHitsTIBHighTof',
        'TrackerHitsTIDLowTof',
        'TrackerHitsTIDHighTof',
        'TrackerHitsTOBLowTof',
        'TrackerHitsTOBHighTof',
        'TrackerHitsTECLowTof',
        'TrackerHitsTECHighTof',
        'TrackerHitsPixelBarrelLowTof',
        'TrackerHitsPixelBarrelHighTof',
        'TrackerHitsPixelEndcapLowTof',
        'TrackerHitsPixelEndcapHighTof'
    ),
    RPCdigisimlinkTag = cms.InputTag("simMuonRPCDigis","RPCDigiSimLink"),
    RPCsimhitsTag = cms.InputTag("g4SimHits","MuonRPCHits"),
    RPCsimhitsXFTag = cms.InputTag("mix","g4SimHitsMuonRPCHits"),
    ThreeHitTracksAreSpecial = cms.bool(False),
    UseGrouped = cms.bool(True),
    UseMuon = cms.bool(False),
    UsePixels = cms.bool(True),
    UseSplitting = cms.bool(True),
    UseTracker = cms.bool(True),
    acceptOneStubMatchings = cms.bool(False),
    associatePixel = cms.bool(True),
    associateRecoTracks = cms.bool(True),
    associateStrip = cms.bool(True),
    associatorByWire = cms.bool(False),
    crossingframe = cms.bool(False),
    dumpDT = cms.bool(False),
    dumpInputCollections = cms.untracked.bool(False),
    ignoreMissingTrackCollection = cms.untracked.bool(True),
    includeZeroHitMuons = cms.bool(True),
    inputCSCSegmentCollection = cms.InputTag("cscSegments"),
    inputDTRecSegment4DCollection = cms.InputTag("dt4DSegments"),
    links_exist = cms.bool(True),
    phase2TrackerSimLinkSrc = cms.InputTag("simSiPixelDigis","Tracker"),
    pixelSimLinkSrc = cms.InputTag("simSiPixelDigis","Pixel"), # Pixel added
    rejectBadGlobal = cms.bool(True),
    simtracksTag = cms.InputTag("g4SimHits"),
    simtracksXFTag = cms.InputTag("mix","g4SimHits"),
    stripSimLinkSrc = cms.InputTag("simSiStripDigis"),
    tpRefVector = cms.bool(True),
    tpTag = cms.InputTag("TPmu"),
    tracksTag = cms.InputTag("hltPhase2L3MuonTracks"),
    useGEMs = cms.bool(True), # changed from False to True
    usePhase2Tracker = cms.bool(True) # changed from False to True
)
```
# Validators

# <span style="color:red">L3 from L1tk</span>
```python
cms.EDProducer("MuonTrackValidator",
    BiDirectional_RecoToSim_association = cms.bool(True),
    UseAssociators = cms.bool(False),
    associatormap = cms.InputTag("tpToL1TkMergedMuonAssociation"),
    associators = cms.vstring('MuonAssociationByHits'),
    beamSpot = cms.InputTag("offlineBeamSpot"),
    dirName = cms.string('HLT/Muon/MuonTrack/'),
    ignoremissingtrackcollection = cms.untracked.bool(True),
    label = cms.VInputTag("hltIter2Phase2L3FromL1TkMuonMerged"),
    label_pileupinfo = cms.InputTag("addPileupInfo"),
    label_tp = cms.InputTag("TPmu"),
    label_tp_refvector = cms.bool(True),
    muonHistoParameters = cms.PSet(
        cotThetaRes_nbin = cms.int32(100),
        cotThetaRes_rangeMax = cms.double(0.01),
        cotThetaRes_rangeMin = cms.double(-0.01),
        do_MUOhitsPlots = cms.bool(False),
        do_TRKhitsPlots = cms.bool(True),
        dxyRes_nbin = cms.int32(100),
        dxyRes_rangeMax = cms.double(0.02),
        dxyRes_rangeMin = cms.double(-0.02),
        dzRes_nbin = cms.int32(100),
        dzRes_rangeMax = cms.double(0.05),
        dzRes_rangeMin = cms.double(-0.05),
        etaRes_nbin = cms.int32(100),
        etaRes_rangeMax = cms.double(0.02),
        etaRes_rangeMin = cms.double(-0.02),
        maxCSCHit = cms.double(50.5),
        maxDTHit = cms.double(50.5),
        maxDxy = cms.double(2.0),
        maxDz = cms.double(30.0),
        maxEta = cms.double(2.5),
        maxFTracks = cms.int32(20),
        maxLayers = cms.double(20.5),
        maxNHit = cms.double(40.5),
        maxNTracks = cms.int32(100),
        maxPU = cms.double(199.5),
        maxPhi = cms.double(3.1416),
        maxPixels = cms.double(5.5),
        maxPt = cms.double(2000.0),
        maxRPCHit = cms.double(10.5),
        maxRpos = cms.double(4.0),
        maxZpos = cms.double(30.0),
        minCSCHit = cms.double(-0.5),
        minDTHit = cms.double(-0.5),
        minDxy = cms.double(-2.0),
        minDz = cms.double(-30.0),
        minEta = cms.double(-2.5),
        minFTracks = cms.int32(0),
        minLayers = cms.double(-0.5),
        minNHit = cms.double(-0.5),
        minNTracks = cms.int32(0),
        minPU = cms.double(-0.5),
        minPhi = cms.double(-3.1416),
        minPixels = cms.double(-0.5),
        minPt = cms.double(0.9),
        minRPCHit = cms.double(-0.5),
        minRpos = cms.double(0.0),
        minZpos = cms.double(-30.0),
        nintCSCHit = cms.int32(51),
        nintDTHit = cms.int32(51),
        nintDxy = cms.int32(40),
        nintDz = cms.int32(60),
        nintEta = cms.int32(50),
        nintFTracks = cms.int32(20),
        nintLayers = cms.int32(21),
        nintNHit = cms.int32(41),
        nintNTracks = cms.int32(100),
        nintPU = cms.int32(100),
        nintPhi = cms.int32(36),
        nintPixels = cms.int32(6),
        nintPt = cms.int32(50),
        nintRPCHit = cms.int32(11),
        nintRpos = cms.int32(40),
        nintZpos = cms.int32(60),
        phiRes_nbin = cms.int32(200),
        phiRes_rangeMax = cms.double(0.01),
        phiRes_rangeMin = cms.double(-0.01),
        ptRes_nbin = cms.int32(200),
        ptRes_rangeMax = cms.double(0.5),
        ptRes_rangeMin = cms.double(-0.5),
        useFabsEta = cms.bool(False),
        useInvPt = cms.bool(False),
        useLogPt = cms.untracked.bool(True),
        usemuon = cms.bool(False),
        usetracker = cms.bool(True)
    ),
    muonTPSelector = cms.PSet(
        chargedOnly = cms.bool(True),
        intimeOnly = cms.bool(True),
        lip = cms.double(30.0),
        maxRapidity = cms.double(2.4),
        minHit = cms.int32(0),
        minRapidity = cms.double(-2.4),
        pdgId = cms.vint32(13, -13),
        ptMax = cms.double(1e+100),
        ptMin = cms.double(0.9),
        signalOnly = cms.bool(True),
        src = cms.InputTag("TPmu"),
        stableOnly = cms.bool(True),
        tip = cms.double(3.5)
    ),
    outputFile = cms.string(''),
    parametersDefiner = cms.string('LhcParametersDefinerForTP'),
    simHitTpMapTag = cms.InputTag("simHitTPAssocProducer"),
    useGEMs = cms.bool(True),
    useME0 = cms.bool(False) # Validation doesn't use ME0
)
```
#  <span style="color:purple">L2 muons from L1 tk</span>
```python
cms.EDProducer("MuonTrackValidator",
    BiDirectional_RecoToSim_association = cms.bool(True),
    UseAssociators = cms.bool(False),
    associatormap = cms.InputTag("tpToL2MuonAssociation"),
    associators = cms.vstring('MuonAssociationByHits'),
    beamSpot = cms.InputTag("offlineBeamSpot"),
    dirName = cms.string('HLT/Muon/MuonTrack/'),
    ignoremissingtrackcollection = cms.untracked.bool(True),
    label = cms.VInputTag("hltL2MuonsFromL1TkMuon:UpdatedAtVtx"),
    label_pileupinfo = cms.InputTag("addPileupInfo"),
    label_tp = cms.InputTag("TPmu"),
    label_tp_refvector = cms.bool(True),
    muonHistoParameters = cms.PSet(
        cotThetaRes_nbin = cms.int32(100),
        cotThetaRes_rangeMax = cms.double(0.1),
        cotThetaRes_rangeMin = cms.double(-0.1),
        do_MUOhitsPlots = cms.bool(True),
        do_TRKhitsPlots = cms.bool(False),
        dxyRes_nbin = cms.int32(100),
        dxyRes_rangeMax = cms.double(10.0),
        dxyRes_rangeMin = cms.double(-10.0),
        dzRes_nbin = cms.int32(100),
        dzRes_rangeMax = cms.double(25.0),
        dzRes_rangeMin = cms.double(-25.0),
        etaRes_nbin = cms.int32(100),
        etaRes_rangeMax = cms.double(0.1),
        etaRes_rangeMin = cms.double(-0.1),
        maxCSCHit = cms.double(50.5),
        maxDTHit = cms.double(50.5),
        maxDxy = cms.double(10.0),
        maxDz = cms.double(30.0),
        maxEta = cms.double(2.5),
        maxFTracks = cms.int32(20),
        maxLayers = cms.double(20.5),
        maxNHit = cms.double(60.5),
        maxNTracks = cms.int32(100),
        maxPU = cms.double(199.5),
        maxPhi = cms.double(3.1416),
        maxPixels = cms.double(5.5),
        maxPt = cms.double(2000.0),
        maxRPCHit = cms.double(10.5),
        maxRpos = cms.double(4.0),
        maxZpos = cms.double(30.0),
        minCSCHit = cms.double(-0.5),
        minDTHit = cms.double(-0.5),
        minDxy = cms.double(-10.0),
        minDz = cms.double(-30.0),
        minEta = cms.double(-2.5),
        minFTracks = cms.int32(0),
        minLayers = cms.double(-0.5),
        minNHit = cms.double(-0.5),
        minNTracks = cms.int32(0),
        minPU = cms.double(-0.5),
        minPhi = cms.double(-3.1416),
        minPixels = cms.double(-0.5),
        minPt = cms.double(0.9),
        minRPCHit = cms.double(-0.5),
        minRpos = cms.double(0.0),
        minZpos = cms.double(-30.0),
        nintCSCHit = cms.int32(51),
        nintDTHit = cms.int32(51),
        nintDxy = cms.int32(40),
        nintDz = cms.int32(60),
        nintEta = cms.int32(50),
        nintFTracks = cms.int32(20),
        nintLayers = cms.int32(21),
        nintNHit = cms.int32(61),
        nintNTracks = cms.int32(100),
        nintPU = cms.int32(100),
        nintPhi = cms.int32(36),
        nintPixels = cms.int32(6),
        nintPt = cms.int32(50),
        nintRPCHit = cms.int32(11),
        nintRpos = cms.int32(40),
        nintZpos = cms.int32(60),
        phiRes_nbin = cms.int32(200),
        phiRes_rangeMax = cms.double(0.1),
        phiRes_rangeMin = cms.double(-0.1),
        ptRes_nbin = cms.int32(200),
        ptRes_rangeMax = cms.double(5.0),
        ptRes_rangeMin = cms.double(-1.0),
        useFabsEta = cms.bool(False),
        useInvPt = cms.bool(False),
        useLogPt = cms.untracked.bool(True),
        usemuon = cms.bool(True),
        usetracker = cms.bool(False)
    ),
    muonTPSelector = cms.PSet(
        chargedOnly = cms.bool(True),
        intimeOnly = cms.bool(True),
        lip = cms.double(30.0),
        maxRapidity = cms.double(2.4),
        minHit = cms.int32(0),
        minRapidity = cms.double(-2.4),
        pdgId = cms.vint32(13, -13),
        ptMax = cms.double(1e+100),
        ptMin = cms.double(0.9),
        signalOnly = cms.bool(True),
        src = cms.InputTag("TPmu"),
        stableOnly = cms.bool(True),
        tip = cms.double(3.5)
    ),
    outputFile = cms.string(''),
    parametersDefiner = cms.string('LhcParametersDefinerForTP'),
    simHitTpMapTag = cms.InputTag("simHitTPAssocProducer"),
    useGEMs = cms.bool(True),
    useME0 = cms.bool(False)
)
```
# <span style="color:lightblue">L3 muons OI</span>
```python
cms.EDProducer("MuonTrackValidator",
    BiDirectional_RecoToSim_association = cms.bool(True),
    UseAssociators = cms.bool(False),
    associatormap = cms.InputTag("tpToL3OIMuonAssociation"),
    associators = cms.vstring('MuonAssociationByHits'),
    beamSpot = cms.InputTag("offlineBeamSpot"),
    dirName = cms.string('HLT/Muon/MuonTrack/'),
    ignoremissingtrackcollection = cms.untracked.bool(True),
    label = cms.VInputTag("hltL3MuonsPhase2L3OI"),
    label_pileupinfo = cms.InputTag("addPileupInfo"),
    label_tp = cms.InputTag("TPmu"),
    label_tp_refvector = cms.bool(True),
    muonHistoParameters = cms.PSet(
        cotThetaRes_nbin = cms.int32(100),
        cotThetaRes_rangeMax = cms.double(0.01),
        cotThetaRes_rangeMin = cms.double(-0.01),
        do_MUOhitsPlots = cms.bool(True),
        do_TRKhitsPlots = cms.bool(True),
        dxyRes_nbin = cms.int32(100),
        dxyRes_rangeMax = cms.double(0.02),
        dxyRes_rangeMin = cms.double(-0.02),
        dzRes_nbin = cms.int32(100),
        dzRes_rangeMax = cms.double(0.05),
        dzRes_rangeMin = cms.double(-0.05),
        etaRes_nbin = cms.int32(100),
        etaRes_rangeMax = cms.double(0.02),
        etaRes_rangeMin = cms.double(-0.02),
        maxCSCHit = cms.double(50.5),
        maxDTHit = cms.double(50.5),
        maxDxy = cms.double(2.0),
        maxDz = cms.double(30.0),
        maxEta = cms.double(2.5),
        maxFTracks = cms.int32(20),
        maxLayers = cms.double(20.5),
        maxNHit = cms.double(80.5),
        maxNTracks = cms.int32(100),
        maxPU = cms.double(199.5),
        maxPhi = cms.double(3.1416),
        maxPixels = cms.double(5.5),
        maxPt = cms.double(2000.0),
        maxRPCHit = cms.double(10.5),
        maxRpos = cms.double(4.0),
        maxZpos = cms.double(30.0),
        minCSCHit = cms.double(-0.5),
        minDTHit = cms.double(-0.5),
        minDxy = cms.double(-2.0),
        minDz = cms.double(-30.0),
        minEta = cms.double(-2.5),
        minFTracks = cms.int32(0),
        minLayers = cms.double(-0.5),
        minNHit = cms.double(-0.5),
        minNTracks = cms.int32(0),
        minPU = cms.double(-0.5),
        minPhi = cms.double(-3.1416),
        minPixels = cms.double(-0.5),
        minPt = cms.double(0.9),
        minRPCHit = cms.double(-0.5),
        minRpos = cms.double(0.0),
        minZpos = cms.double(-30.0),
        nintCSCHit = cms.int32(51),
        nintDTHit = cms.int32(51),
        nintDxy = cms.int32(40),
        nintDz = cms.int32(60),
        nintEta = cms.int32(50),
        nintFTracks = cms.int32(20),
        nintLayers = cms.int32(21),
        nintNHit = cms.int32(81),
        nintNTracks = cms.int32(100),
        nintPU = cms.int32(100),
        nintPhi = cms.int32(36),
        nintPixels = cms.int32(6),
        nintPt = cms.int32(50),
        nintRPCHit = cms.int32(11),
        nintRpos = cms.int32(40),
        nintZpos = cms.int32(60),
        phiRes_nbin = cms.int32(200),
        phiRes_rangeMax = cms.double(0.01),
        phiRes_rangeMin = cms.double(-0.01),
        ptRes_nbin = cms.int32(200),
        ptRes_rangeMax = cms.double(0.5),
        ptRes_rangeMin = cms.double(-0.5),
        useFabsEta = cms.bool(False),
        useInvPt = cms.bool(False),
        useLogPt = cms.untracked.bool(True),
        usemuon = cms.bool(True),
        usetracker = cms.bool(True)
    ),
    muonTPSelector = cms.PSet(
        chargedOnly = cms.bool(True),
        intimeOnly = cms.bool(True),
        lip = cms.double(30.0),
        maxRapidity = cms.double(2.4),
        minHit = cms.int32(0),
        minRapidity = cms.double(-2.4),
        pdgId = cms.vint32(13, -13),
        ptMax = cms.double(1e+100),
        ptMin = cms.double(0.9),
        signalOnly = cms.bool(True),
        src = cms.InputTag("TPmu"),
        stableOnly = cms.bool(True),
        tip = cms.double(3.5)
    ),
    outputFile = cms.string(''),
    parametersDefiner = cms.string('LhcParametersDefinerForTP'),
    simHitTpMapTag = cms.InputTag("simHitTPAssocProducer"),
    useGEMs = cms.bool(True),
    useME0 = cms.bool(False)
)
```
# <span style="color:Gold">L3 muon merged</span>
```python
cms.EDProducer("MuonTrackValidator",
    BiDirectional_RecoToSim_association = cms.bool(True),
    UseAssociators = cms.bool(False),
    associatormap = cms.InputTag("tpToL3MuonMergedAssociation"),
    associators = cms.vstring('MuonAssociationByHits'),
    beamSpot = cms.InputTag("offlineBeamSpot"),
    dirName = cms.string('HLT/Muon/MuonTrack/'),
    ignoremissingtrackcollection = cms.untracked.bool(True),
    label = cms.VInputTag("hltPhase2L3MuonMerged"),
    label_pileupinfo = cms.InputTag("addPileupInfo"),
    label_tp = cms.InputTag("TPmu"),
    label_tp_refvector = cms.bool(True),
    muonHistoParameters = cms.PSet(
        cotThetaRes_nbin = cms.int32(100),
        cotThetaRes_rangeMax = cms.double(0.01),
        cotThetaRes_rangeMin = cms.double(-0.01),
        do_MUOhitsPlots = cms.bool(False),
        do_TRKhitsPlots = cms.bool(True),
        dxyRes_nbin = cms.int32(100),
        dxyRes_rangeMax = cms.double(0.02),
        dxyRes_rangeMin = cms.double(-0.02),
        dzRes_nbin = cms.int32(100),
        dzRes_rangeMax = cms.double(0.05),
        dzRes_rangeMin = cms.double(-0.05),
        etaRes_nbin = cms.int32(100),
        etaRes_rangeMax = cms.double(0.02),
        etaRes_rangeMin = cms.double(-0.02),
        maxCSCHit = cms.double(50.5),
        maxDTHit = cms.double(50.5),
        maxDxy = cms.double(2.0),
        maxDz = cms.double(30.0),
        maxEta = cms.double(2.5),
        maxFTracks = cms.int32(20),
        maxLayers = cms.double(20.5),
        maxNHit = cms.double(40.5),
        maxNTracks = cms.int32(100),
        maxPU = cms.double(199.5),
        maxPhi = cms.double(3.1416),
        maxPixels = cms.double(5.5),
        maxPt = cms.double(2000.0),
        maxRPCHit = cms.double(10.5),
        maxRpos = cms.double(4.0),
        maxZpos = cms.double(30.0),
        minCSCHit = cms.double(-0.5),
        minDTHit = cms.double(-0.5),
        minDxy = cms.double(-2.0),
        minDz = cms.double(-30.0),
        minEta = cms.double(-2.5),
        minFTracks = cms.int32(0),
        minLayers = cms.double(-0.5),
        minNHit = cms.double(-0.5),
        minNTracks = cms.int32(0),
        minPU = cms.double(-0.5),
        minPhi = cms.double(-3.1416),
        minPixels = cms.double(-0.5),
        minPt = cms.double(0.9),
        minRPCHit = cms.double(-0.5),
        minRpos = cms.double(0.0),
        minZpos = cms.double(-30.0),
        nintCSCHit = cms.int32(51),
        nintDTHit = cms.int32(51),
        nintDxy = cms.int32(40),
        nintDz = cms.int32(60),
        nintEta = cms.int32(50),
        nintFTracks = cms.int32(20),
        nintLayers = cms.int32(21),
        nintNHit = cms.int32(41),
        nintNTracks = cms.int32(100),
        nintPU = cms.int32(100),
        nintPhi = cms.int32(36),
        nintPixels = cms.int32(6),
        nintPt = cms.int32(50),
        nintRPCHit = cms.int32(11),
        nintRpos = cms.int32(40),
        nintZpos = cms.int32(60),
        phiRes_nbin = cms.int32(200),
        phiRes_rangeMax = cms.double(0.01),
        phiRes_rangeMin = cms.double(-0.01),
        ptRes_nbin = cms.int32(200),
        ptRes_rangeMax = cms.double(0.5),
        ptRes_rangeMin = cms.double(-0.5),
        useFabsEta = cms.bool(False),
        useInvPt = cms.bool(False),
        useLogPt = cms.untracked.bool(True),
        usemuon = cms.bool(False),
        usetracker = cms.bool(True)
    ),
    muonTPSelector = cms.PSet(
        chargedOnly = cms.bool(True),
        intimeOnly = cms.bool(True),
        lip = cms.double(30.0),
        maxRapidity = cms.double(2.4),
        minHit = cms.int32(0),
        minRapidity = cms.double(-2.4),
        pdgId = cms.vint32(13, -13),
        ptMax = cms.double(1e+100),
        ptMin = cms.double(0.9),
        signalOnly = cms.bool(True),
        src = cms.InputTag("TPmu"),
        stableOnly = cms.bool(True),
        tip = cms.double(3.5)
    ),
    outputFile = cms.string(''),
    parametersDefiner = cms.string('LhcParametersDefinerForTP'),
    simHitTpMapTag = cms.InputTag("simHitTPAssocProducer"),
    useGEMs = cms.bool(True),
    useME0 = cms.bool(False)
)
```
# <span style="color:ForestGreen">L3 Muon Tracks</span>
```python
cms.EDProducer("MuonTrackValidator",
    BiDirectional_RecoToSim_association = cms.bool(True),
    UseAssociators = cms.bool(False),
    associatormap = cms.InputTag("tpToL3TkMuonAssociation"),
    associators = cms.vstring('MuonAssociationByHits'),
    beamSpot = cms.InputTag("offlineBeamSpot"),
    dirName = cms.string('HLT/Muon/MuonTrack/'),
    ignoremissingtrackcollection = cms.untracked.bool(True),
    label = cms.VInputTag("hltPhase2L3MuonTracks"),
    label_pileupinfo = cms.InputTag("addPileupInfo"),
    label_tp = cms.InputTag("TPmu"),
    label_tp_refvector = cms.bool(True),
    muonHistoParameters = cms.PSet(
        cotThetaRes_nbin = cms.int32(100),
        cotThetaRes_rangeMax = cms.double(0.01),
        cotThetaRes_rangeMin = cms.double(-0.01),
        do_MUOhitsPlots = cms.bool(False),
        do_TRKhitsPlots = cms.bool(True),
        dxyRes_nbin = cms.int32(100),
        dxyRes_rangeMax = cms.double(0.02),
        dxyRes_rangeMin = cms.double(-0.02),
        dzRes_nbin = cms.int32(100),
        dzRes_rangeMax = cms.double(0.05),
        dzRes_rangeMin = cms.double(-0.05),
        etaRes_nbin = cms.int32(100),
        etaRes_rangeMax = cms.double(0.02),
        etaRes_rangeMin = cms.double(-0.02),
        maxCSCHit = cms.double(50.5),
        maxDTHit = cms.double(50.5),
        maxDxy = cms.double(2.0),
        maxDz = cms.double(30.0),
        maxEta = cms.double(2.5),
        maxFTracks = cms.int32(20),
        maxLayers = cms.double(20.5),
        maxNHit = cms.double(40.5),
        maxNTracks = cms.int32(100),
        maxPU = cms.double(199.5),
        maxPhi = cms.double(3.1416),
        maxPixels = cms.double(5.5),
        maxPt = cms.double(2000.0),
        maxRPCHit = cms.double(10.5),
        maxRpos = cms.double(4.0),
        maxZpos = cms.double(30.0),
        minCSCHit = cms.double(-0.5),
        minDTHit = cms.double(-0.5),
        minDxy = cms.double(-2.0),
        minDz = cms.double(-30.0),
        minEta = cms.double(-2.5),
        minFTracks = cms.int32(0),
        minLayers = cms.double(-0.5),
        minNHit = cms.double(-0.5),
        minNTracks = cms.int32(0),
        minPU = cms.double(-0.5),
        minPhi = cms.double(-3.1416),
        minPixels = cms.double(-0.5),
        minPt = cms.double(0.9),
        minRPCHit = cms.double(-0.5),
        minRpos = cms.double(0.0),
        minZpos = cms.double(-30.0),
        nintCSCHit = cms.int32(51),
        nintDTHit = cms.int32(51),
        nintDxy = cms.int32(40),
        nintDz = cms.int32(60),
        nintEta = cms.int32(50),
        nintFTracks = cms.int32(20),
        nintLayers = cms.int32(21),
        nintNHit = cms.int32(41),
        nintNTracks = cms.int32(100),
        nintPU = cms.int32(100),
        nintPhi = cms.int32(36),
        nintPixels = cms.int32(6),
        nintPt = cms.int32(50),
        nintRPCHit = cms.int32(11),
        nintRpos = cms.int32(40),
        nintZpos = cms.int32(60),
        phiRes_nbin = cms.int32(200),
        phiRes_rangeMax = cms.double(0.01),
        phiRes_rangeMin = cms.double(-0.01),
        ptRes_nbin = cms.int32(200),
        ptRes_rangeMax = cms.double(0.5),
        ptRes_rangeMin = cms.double(-0.5),
        useFabsEta = cms.bool(False),
        useInvPt = cms.bool(False),
        useLogPt = cms.untracked.bool(True),
        usemuon = cms.bool(False),
        usetracker = cms.bool(True)
    ),
    muonTPSelector = cms.PSet(
        chargedOnly = cms.bool(True),
        intimeOnly = cms.bool(True),
        lip = cms.double(30.0),
        maxRapidity = cms.double(2.4),
        minHit = cms.int32(0),
        minRapidity = cms.double(-2.4),
        pdgId = cms.vint32(13, -13),
        ptMax = cms.double(1e+100),
        ptMin = cms.double(0.9),
        signalOnly = cms.bool(True),
        src = cms.InputTag("TPmu"),
        stableOnly = cms.bool(True),
        tip = cms.double(3.5)
    ),
    outputFile = cms.string(''),
    parametersDefiner = cms.string('LhcParametersDefinerForTP'),
    simHitTpMapTag = cms.InputTag("simHitTPAssocProducer"),
    useGEMs = cms.bool(True),
    useME0 = cms.bool(False)
)
```