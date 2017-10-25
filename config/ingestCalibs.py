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

config.parse.translation = {'calibDate': 'calibDate',
                            }

config.parse.translators = {#'ccd': 'translate_ccd',
                            #'filter': 'translate_filter',
                            #'calibDate': 'translate_calibDate',
                            }

config.register.detector = ['arm', 'spectrograph']
config.register.unique = ['arm', 'spectrograph', 'calibDate']
config.register.tables = ['bias', 'dark', 'detectormap', 'flat', 'fibertrace'] # n.b. lower case
config.register.visit = ['calibDate', 'arm', 'spectrograph']
