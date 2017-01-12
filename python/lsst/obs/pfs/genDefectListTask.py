#!/usr/bin/env python
#
# LSST Data Management System
# Copyright 2008, 2009, 2010, 2011 LSST Corporation.
#
# This product includes software developed by the
# LSST Project (http://www.lsst.org/).
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the LSST License Statement and
# the GNU General Public License along with this program.  If not,
# see <http://www.lsstcorp.org/LegalNotices/>.
#
"""
This module defines the CmdLineTask class GenDefectListTask
@author Andreas Ritter, Princeton University
"""
import lsst.afw.display as afwDisplay
import lsst.afw.image as afwImage
import lsst.afw.math as afwMath
from lsst.log import Log
from lsst.pex.config import Config, Field, ListField
from lsst.pipe.base import Struct, TaskRunner, ArgumentParser, CmdLineTask
import lsst.utils
import numpy as np
import os

class GenDefectListConfig(Config):
    badCols = ListField(doc="List of bad columns",
                        dtype=int,
                        default=[1023,1024,2045,2046,2047,2048,2049,2050,2051,2052,3071,3072])
    display = Field(doc="Set to display outputs",
                    dtype=bool,
                    default=False)
    expTimeHigh = Field(doc="Exposure time for stronger flats",
                        dtype=int,
                        default=30)
    expTimeLow = Field(doc="Exposure time for weaker flats",
                       dtype=int,
                       default=2)
    maxSigma = Field(doc="Maximum sigma for good pixels",
                     dtype=float,
                     default=6.)
    nPixCut = Field(doc="Number of rows/columns around the edge of images to ignore",
                    dtype=int,
                    default=40)
    date = Field(doc="date for output defect list",
                 dtype=str,
                 default="2015-12-01")

class GenDefectListTaskRunner(TaskRunner):
    """Get parsed values into the GenDefectListTask.run"""
    @staticmethod
    def getTargetList(parsedCmd, **kwargs):
        return [dict(expRefList=parsedCmd.id.refList, butler=parsedCmd.butler)]

    def __call__(self, args):
        task = self.TaskClass(config=self.config, log=self.log)
        if self.doRaise:
            self.log.info('args = %s' % args)
            result = task.run(**args)
        else:
            try:
                result = task.run(**args)
            except Exception, e:
                task.log.fatal("Failed with error message <%s>" % e)

        if self.doReturnResults:
            return Struct(
                args = args,
                metadata = task.metadata,
                result = result,
            )

class GenDefectListTask(CmdLineTask):
    """
    @brief Task to generate the defect list for one detector
    This task writes the file $OBS_PFS_DIR/pfs/defects/<date>/defects.dat.
    All pixels which are not linear are written to the defect list.
    Non-linearity is assumed if a pixel value is less/greater than the median
    value of the quotient of the post ISR images from all Flats with the
    longer exposure time and all Flats with the shorter exposure time, multiplied
    with -/+ the argument <maxSigma>.
    Bad pixels from other detectors are read via butler.get('defects')
    and included in the output file.
    Note that currently only individual pixels are written to the defects list
    (width=height=1), no clustering is done.

    Public methods:

    __init__(self, *args, **kwargs)
    getMedianFlats(self)
    run(self, expRefList, butler)
    """
    ConfigClass = GenDefectListConfig
    RunnerClass = GenDefectListTaskRunner
    _DefaultName = "genDefectListTask"

    def __init__(self, *args, **kwargs):
        """create and initialize the shared data"""

        """Initialize super class"""
        CmdLineTask.__init__(self, *args, **kwargs)

        """Initialize logger and set log level"""
        self.logger = Log.getLogger("genDefectList")
        self.logger.setLevel(self.log.getLevel())

        """Define exposures which the user might want to display or save,
           returned by the function self.getMedianFlats()"""
        self.medianLowExp = 0
        self.medianHighExp = 0
        self.divFlatExp = 0

        """Define variable to contain the detector ID for which to find bad pixels"""
        self.ccdNum = 0

    @classmethod
    def _makeArgumentParser(cls, *args, **kwargs):
        parser = ArgumentParser(name=cls._DefaultName)

        """TODO: If 'output' or 'rerun' are specified, write postISRCCD images there"""

        """We don't require the 'output' or 'rerun' parameters specified as we don't need
           to write output files except for the file 'defects.dat', the directory in which
           to write is solely specified by the 'date' parameter."""
        parser.requireOutput = False

        """Add ID argument to identify flats from which to identify the bad pixels"""
        parser.add_id_argument("--id", datasetType="postISRCCD",
                               help="input identifiers, e.g., --id visit=123 ccd=4")
        return parser

    def _getOutputString(self, ccdNum, col, row, width, height):
        """ Return a string which contains one defect in the output file 'defects.dat'
        @param ccdNum: detector ID for which to create the output string
        @param col: x position of defect start
        @param row: y position of defect start
        @param width: width of defect
        @param height: height of defect
        @returns a string which contains one defect in the output file 'defects.dat'
        """
        return (repr(ccdNum).ljust(5) + repr(col).ljust(6) +
                repr(row).ljust(6) + repr(width).ljust(6) +
                repr(height).ljust(6) + '\n')

    def _getDetectorId(self, butler, root, spectrograph, arm):
        """ Calculate and return the detector ID from 'arm' and 'spectrograph'
            Should just call pfsMapper._extractDetectorId() but I haven't figured
            out yet how to do this
        @param spectrograph: Spectrograph Number [1..4]
        @param arm: Arm ['b', 'r', 'm', 'n']
        @returns Spectrograph Id from lsst/obs/pfs/camera.py
        """
        mapper = butler.getMapperClass(root)(root=root)
        dataId = dict(spectrograph=spectrograph, arm=arm)
        return mapper._extractDetectorId(dataId)

    def getMedianFlats(self):
        """ Return exposures containing the median of flats with lower and higher
            exposure time, and the quotient medianHigh / medianLow"""
        return {'medianFlatsLow': self.medianLowExp,
                'medianFlatsHigh': self.medianHighExp,
                'divFlat': self.divFlatExp}

    def _getFlatsMedian(self, butler, spectrograph, arm, visits, expTime):
        """ Calculate the median of the images with the visit numbers listed in the
            parameter 'visits'
        @param butler: Butler to use for retrieving the images
        @param spectrograph: Spectrograph to use [1..4]
        @param arm: Arm to use ['b', 'r', 'm', 'n']
        @param visits: list of visit numbers to calculate median from
        @param expTime: exposure time of flats with visit numbers in visits
        @returns Median of flats visit visit numbers in visits
        """
        dataId = dict(visit=visits[0],
                      spectrograph=spectrograph,
                      arm=arm)
        flat = butler.get('postISRCCD', dataId).getMaskedImage().getImage().getArray()
        allFlats = np.ndarray(shape=[flat.shape[0], flat.shape[1], len(visits)], dtype='float32')
        allFlats[:,:,0] = flat
        for i in np.arange(1, len(visits)):
            dataId['visit'] = visits[i]
            allFlats[:,:,i] = butler.get('postISRCCD', dataId).getMaskedImage().getImage().getArray()

        nFlats = len(visits)
        visitsStr = ', '.join(str(e) for e in visits)
        if nFlats > 2:
            medianFlats = np.median(allFlats, 2)
            self.logger.info('%d Flats with %ds exposure time found with visits %s' % (nFlats,
                                                                                       expTime,
                                                                                       visitsStr))
        elif nFlats == 2:
            medianFlats = np.min(allFlats, 2)
            self.logger.warn('Only %d Flats with %ds exposure time found with visits %s' % (nFlats,
                                                                                            expTime,
                                                                                            visitsStr))
        elif nFlats == 0:
            raise RunTimeError("No Flats with %d s exposure time found" % self.config.expTimeLow)
        else:
            self.logger.warn('Only 1 Flat with %ds exposure time found with visit %s' % (expTime,
                                                                                         visitsStr))
            medianFlats = allFlats[:,:,0]
        medianFlats = np.where( medianFlats == 0., 0.000001, medianFlats)
        return medianFlats

    def run(self, expRefList, butler):
        """ Create defect list
            The output will be written to '$OBS_PFS_DIR/pfs/defects/<self.config.date>/defects.dat'.
            Pre-existing bad pixels from other detectors will be added to the file.
        @param expRefList: List of exposure references to Flats from which to find the defects
        @param butler: Butler to use
        """

        if len(expRefList) == 0:
            raise RuntimeError("no exposures found matching the ID")

        """ Quiet down lots of loggers, so we can see genDefectList logs better """
        loggerNames = ["daf.persistence.butler",
                       "afw.image.Mask",
                       "afw.ExposureFormatter",
                       "daf.persistence.LogicalLocation",
                       "CameraMapper",
                       "afw.image.ExposureInfo"]
        for loggerName in loggerNames:
            Log.getLogger(loggerName).setLevel(Log.WARN)

        dataRef = expRefList[0]
        root = os.path.dirname(os.path.dirname(butler.get('raw_filename', dataRef.dataId)[0]))
        print 'root = ',root
        self.ccdNum = self._getDetectorId(butler, root, dataRef.dataId['spectrograph'], dataRef.dataId['arm'])

        rawMaskedImage = butler.get('postISRCCD', dataRef.dataId).getMaskedImage()
        nRows = rawMaskedImage.getHeight()
        nCols = rawMaskedImage.getWidth()
        self.logger.debug("nRows = %d, nCols = %d" % (nRows, nCols))

        """check arguments for validity"""
        if self.config.expTimeLow <= 0:
            raise RuntimeError("expTimeLow needs to be greater than 0")
        if self.config.expTimeHigh <= self.config.expTimeLow:
            raise RuntimeError("expTimeHigh needs to be greater than expTimeLow")
        if self.config.nPixCut < 0:
            raise RuntimeError("nPixCut needs to be greater than or equal to 0")

        outFileName = os.path.join(lsst.utils.getPackageDir("obs_pfs"),
                                   "pfs/defects",
                                   self.config.date,
                                   "defects.dat")
        self.logger.debug('outFileName = <%s>' % outFileName)

        flatsLow = list()
        flatsHigh = list()
        for dataRefTemp in expRefList:
            visit = dataRefTemp.dataId['visit']
            md = butler.get('raw_md', dataRefTemp.dataId)
            if md.get("IMAGETYP") == "flat":
                if md.get('EXPTIME') == self.config.expTimeLow:
                    flatsLow.append(visit)
                elif md.get('EXPTIME') == self.config.expTimeHigh:
                    flatsHigh.append(visit)

        self.logger.debug("len(flatsLow) = %d" % len(flatsLow))
        self.logger.debug("len(flatsHigh) = %d" % len(flatsHigh))
        self.logger.debug("detector Id = %d" % self.ccdNum)

        medianLow = self._getFlatsMedian(butler,
                                         dataRef.dataId['spectrograph'],
                                         dataRef.dataId['arm'],
                                         flatsLow,
                                         self.config.expTimeLow)
        medianHigh = self._getFlatsMedian(butler,
                                          dataRef.dataId['spectrograph'],
                                          dataRef.dataId['arm'],
                                          flatsHigh,
                                          self.config.expTimeHigh)

        divFlat = medianHigh / medianLow

        sctrl = afwMath.StatisticsControl()
        sctrl.setWeighted(False)

        divFlatIm = afwImage.ImageF(divFlat)
        divFlatMIm = afwImage.makeMaskedImage(divFlatIm)

        stat = afwMath.makeStatistics(divFlatMIm, afwMath.MEANCLIP | afwMath.STDEVCLIP)
        mean = stat.getValue(afwMath.MEANCLIP)
        stddev = stat.getValue(afwMath.STDEVCLIP)

        mask = divFlatMIm.getMask()
        maskVal = 1 << divFlatMIm.getMask().getMaskPlane('BAD')

        """Write defects file"""
        if self.config.date != '':
            with open(outFileName, "w") as text_file:
                text_file.write("""#\n""")
                text_file.write("""# Defects file\n""")
                text_file.write("""#\n""")
                text_file.write("""# Convert to a fits table using getDefectFits.py;  """)
                text_file.write("""this is done for you by running scons\n""")
                text_file.write("""#\n""")
                text_file.write("""#CCD x0    y0    width height\n""")

                for spectrograph in np.arange(1,5):
                    for arm in ['b', 'r', 'n']:
                        ccdNum = self._getDetectorId(butler, root, spectrograph, arm)
                        self.log.debug("spectrograph = %d, arm = %s: ccdNum = %d" % (spectrograph,
                                                                                     arm,
                                                                                     ccdNum))
                        """ if detector ID is equal to detector ID for which to find the bad pixels,
                            find bad pixels and write them to the output file, otherwise write pre-existing
                            bad pixels for the other detectors """
                        if ccdNum == self.ccdNum:
                            for badCol in self.config.badCols:
                                s = self._getOutputString(ccdNum,
                                                          badCol,
                                                          0,
                                                          1,
                                                          nRows)
                                text_file.write(s)
                                mask[badCol,:] |= maskVal
                            nBad = 0
                            for iRow in np.arange(self.config.nPixCut,nRows-self.config.nPixCut):
                                for iCol in np.arange(self.config.nPixCut,nCols-self.config.nPixCut):
                                    if iCol not in self.config.badCols:
                                        if (np.absolute(divFlat[iRow, iCol] - mean)
                                            > self.config.maxSigma * stddev):
                                            nBad = nBad + 1
                                            s = self._getOutputString(ccdNum,
                                                                      iCol,
                                                                      iRow,
                                                                      1,
                                                                      1)
                                            text_file.write(s)
                                            mask[iCol, iRow] |= maskVal;
                        else:
                            defectList = butler.get("defects",
                                                    spectrograph=spectrograph,
                                                    arm=arm,
                                                    visit=dataRef.dataId['visit'],
                                                    )
                            self.log.debug("spectrograph = %d, arm = %s: len(defectList) = %d"
                                           % (spectrograph,
                                              arm,
                                              len(defectList)))
                            for defect in defectList:
                                bbox = defect.getBBox()
                                s = self._getOutputString(ccdNum,
                                                          bbox.getBeginX(),
                                                          bbox.getBeginY(),
                                                          bbox.getWidth(),
                                                          bbox.getHeight())
                                text_file.write(s)

                text_file.write("# Dead amps")

            self.logger.info('%d bad pixels found' % nBad)

        self.medianLowExp = afwImage.makeExposure(afwImage.makeMaskedImage(afwImage.ImageF(medianLow)))
        self.medianLowExp.getMaskedImage().getMask()[:,:] = mask[:,:]
        self.medianHighExp = afwImage.makeExposure(afwImage.makeMaskedImage(afwImage.ImageF(medianHigh)))
        self.medianHighExp.getMaskedImage().getMask()[:,:] = mask[:,:]
        self.divFlatExp = afwImage.makeExposure(divFlatMIm)

        if self.config.display:
            afwDisplay.setDefaultMaskTransparency(75)
            display0 = afwDisplay.getDisplay()
            display0.mtv(afwImage.makeExposure(divFlatMIm), title="parent")
