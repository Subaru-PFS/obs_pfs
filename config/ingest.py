from lsst.obs.pfs.ingest import PfsParseTask
config.parse.retarget(PfsParseTask)
config.register.columns = {'site': 'text', #J: JHU, L: LAM, X: Subaru offline, I: IPMU, A: ASIAA, S: Summit, P: Princeton, F: simulation (fake)
                           'category': 'text', #A: science, B: NTR, C: Meterology, D: HG
                           'field': 'text', # IMAGETYP
                           'visit': 'int',
                           'ccd': 'int', #[0-11]
                           'filter': 'text', #b: blue, r: red, n: nir, m: medium resolution red
                           'arm': 'text', #b: blue, r: red, n: nir, m: medium resolution red
                           'spectrograph': 'int', #1-4
                           'dateObs': 'text',
                           'expTime': 'double',
                           'dataType': 'text',
                           'taiObs': 'text',
                           'pfsConfigId': 'int',
                          }
config.register.unique = ['site', 'category', 'visit', 'filter', 'arm', 'spectrograph', 'pfsConfigId']
config.register.visit = ['visit', 'field', 'filter', 'spectrograph', 'arm', 'dateObs', 'taiObs', 'pfsConfigId']

config.parse.translation = {'dataType': 'IMAGETYP',
                            'expTime': 'EXPTIME',
                            'dateObs': 'DATE-OBS',
                            'taiObs': 'DATE-OBS',
}
config.parse.defaults = {'ccdTemp': "0", # Added in commissioning run 3
                       }
config.parse.translators = {'field': 'translate_field',
                            'dateObs': 'translate_date',
                            'taiObs': 'translate_date',
                            'pfsConfigId': 'translate_pfsConfigId',
}
