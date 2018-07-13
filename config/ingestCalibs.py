from lsst.obs.pfs.ingest import PfsCalibsParseTask
config.parse.retarget(PfsCalibsParseTask)

config.register.columns = {'arm': 'text',
                           'ccd': 'int',
                           'spectrograph': 'int',
                           'calibDate': 'text',
                           'visit0': 'int',
                           'validStart': 'text',
                           'validEnd': 'text',
                           }

config.parse.translators = {'arm': 'translate_arm',
                            'ccd': 'translate_ccd',
                            'spectrograph': 'translate_spectrograph',
                            'calibDate': 'translate_calibDate',
                            'visit0': 'translate_visit0',
                            }

config.register.detector = ['arm', 'spectrograph']
config.register.unique = ['arm', 'spectrograph', 'calibDate']
config.register.tables = ['bias', 'dark', 'detectormap', 'flat', 'fibertrace']
config.register.visit = ['calibDate', 'arm', 'spectrograph']
