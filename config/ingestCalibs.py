from lsst.obs.pfs.ingest import PfsCalibsParseTask
config.parse.retarget(PfsCalibsParseTask)
config.register.validityUntilSuperseded = [
    'bias', 'dark', 'detectorMap', 'flat', 'fiberProfiles', 'fiberNorms' 'defects'
]

config.register.columns = {'arm': 'text',
                           'ccd': 'int',
                           'spectrograph': 'int',
                           'calibDate': 'text',  # Date of calib; used for filenames
                           'calibTime': 'text',  # Date+time of calib; used for setting validity ranges
                           'visit0': 'int',
                           'validStart': 'text',
                           'validEnd': 'text',
                           }

config.parse.translators = {'arm': 'translate_arm',
                            'ccd': 'translate_ccd',
                            'spectrograph': 'translate_spectrograph',
                            'calibDate': 'translate_calibDate',
                            'calibTime': 'translate_calibTime',
                            'visit0': 'translate_visit0',
                            }

config.register.detector = ['arm', 'spectrograph']
config.register.unique = ['arm', 'spectrograph', 'calibTime']
config.register.tables = ['bias', 'dark', 'detectorMap', 'flat', 'fiberProfiles', 'fiberNorms', 'ipc']
config.register.visit = ['calibDate', 'arm', 'spectrograph']
config.register.calibDate = "calibTime"  # date+time, for intra-day resolution
