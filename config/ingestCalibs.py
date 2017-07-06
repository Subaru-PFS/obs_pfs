from lsst.obs.pfs.ingest import PfsCalibsParseTask
config.parse.retarget(PfsCalibsParseTask)

config.register.columns = {'filter': 'text',
                           'ccd': 'int',
                           'calibDate': 'text',
                           'validStart': 'text',
                           'validEnd': 'text',
                           }

config.parse.translators = {#'ccd': 'translate_ccd',
                            #'filter': 'translate_filter',
                            #'calibDate': 'translate_calibDate',
                            }

config.register.unique = ['filter', 'ccd', 'calibDate']
config.register.tables = ['arc', 'bias', 'dark', 'flat', 'fibertrace']
config.register.visit = ['calibDate', 'filter']
