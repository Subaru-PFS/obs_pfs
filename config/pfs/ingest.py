from lsst.obs.subaru.ingest import PfsParseTask
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
                           #IMAGETYP, same as field
                           #'src': 'text' #this is new
                          }
config.register.unique = ['site', 'category', 'visit', 'filter', 'arm', 'spectrograph']
config.register.visit = ['visit', 'field', 'filter', 'arm', 'dateObs', 'taiObs']

config.parse.translation = {'dataType': 'IMAGETYP',
                            'expTime': 'EXPTIME',
                          #'ccd': 'DET-ID',
                          #'pa': 'INST-PA',
                          #'autoguider': 'T_AG',
                          #'ccdTemp': 'T_CCDTV',
                          #'config': 'T_CFGFIL',
                          #'frameId': 'FRAMEID',
                          #'expId': 'EXP-ID',
                            'dateObs': 'DATE-OBS',
                            'taiObs': 'DATE-OBS',
}
config.parse.defaults = {'ccdTemp': "0", # Added in commissioning run 3
                       }
config.parse.translators = {'field': 'translate_field',
                            'dateObs': 'translate_date',
                            'taiObs': 'translate_date'
                          #'visit': 'translate_visit',
                          #'pointing': 'translate_pointing',
                          #'filter': 'translate_filter',
}

