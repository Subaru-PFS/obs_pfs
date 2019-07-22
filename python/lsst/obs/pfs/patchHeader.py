#!/usr/bin/env python3

import inspect
import logging
import pathlib

import astropy.io.fits as pyfits

# ruaml_yaml is a bit more careful. In particular, it at least
# _complains_ about duplicate top-level keys.
try:
    import ruamel_yaml as yaml
except ImportError:
    import yaml

logger = logging.getLogger('patchHeader')

def _cardWithComment(value, comment=None):
    # Stupid recent astropy.io.fits change got rid of comment=comment
    if comment is None:
        return value
    else:
        return value, comment

def addCard(hdr, *, name, value,
            comment=None, overwriteExisting=False):
    """ Add a card to a FITS header

    Args
    ----
    hdr : pyfits header dictionary
       The header to modify. This is modified in place.
    name : str
       The card name to add
    value : object
       The value to add. Typing is done by pyfits.
    comment : str
       A comment for the new card.
    overwriteExisting : bool
       If True, overwrite existing cards, else raise.

    Returns
    -------
    hdr : fits header dictionary
       Returned just for convenience: the hdr is modified in place.

    If you want to add or overwrite, use modifyCard instead.
    """

    if not overwriteExisting and name in hdr:
        raise KeyError(f"card {name} already exists")

    hdr[name] = _cardWithComment(value, comment)

    return hdr

def deleteCard(hdr, *, name):
    """ Delete a card from a FITS header

    Args
    ----
    hdr : pyfits header dictionary
       The header to modify. This is modified in place.
    name : str
       The card name to delete
    """

    hdr.remove(name)
    return hdr

def modifyCard(hdr, *, name, newValue,
               comment=None,
               unlessValueIs=None, onlyIfValueIs=None,
               checkTypes=True,
               addIfMissing=False):
    """ Modify an existing FITS card.

    Args
    ----
    hdr : pyfits header dictionary
       The header to modify. This is modified in place.
    name : str
       The card name to modify
    newValue : object
       The value to replace any existing value with. Typing is done by pyfits.
    unlessValueIs: object
       If not None, only replace the value if it is NOT this.
    onlyIfValueIs: object
       If not None, only replace the value if it IS this.
    addIfMissing: False
       Do addCard(hdr, name, newValue) if the card does not exist yet.

    Returns
    -------
    hdr : fitsio.FITSHDR
       Returned just for convenience: the hdr is modified in place.

    """
    if name not in hdr:
        if addIfMissing:
            hdr[name] = _cardWithComment(newValue, comment)
        return hdr

    if unlessValueIs is not None and onlyIfValueIs is not None:
        raise RuntimeError("unlessValueIs and onlyIfValueIs conflict")

    if unlessValueIs is not None:
        if unlessValueIs != hdr[name]:
            hdr[name] = _cardWithComment(newValue, comment)
    elif onlyIfValueIs is not None:
        if onlyIfValueIs == hdr[name]:
            hdr[name] = _cardWithComment(newValue, comment)
    else:
        hdr[name] = _cardWithComment(newValue, comment)

__allRules = {addCard, deleteCard, modifyCard}
__allRulesByName = {r.__name__: r for r in __allRules}

def _checkOrRunRule(hdr, rule, doRun=True):
    """ Check a rule for validity. Raises an exception if not. """

    matchedRule = None
    if 'action' in rule:
        matchedRule = rule
    else:
        # Syntactic sugar: convert "actionName: foo" into
        #   action: actionName
        #   name: foo
        # and allow "foo" to be optional.
        for r in __allRules:
            if r.__name__ in rule:
                if 'action' in rule:
                    raise RuntimeError(f"cannot have action key as well as {r.__name__}")
                rule['action'] = r.__name__
                sugarValue = rule.pop(r.__name__)
                if sugarValue is not None:
                    raise RuntimeError(f"both action: {r.__name__} and {r.__name__}: {sugarValue} are set")

                raise NotImplementedError("rule: cardName not implemented, sorry.")

                if 'name' in rule:
                    raise RuntimeError(f"both {r.__name__} and {r.__name__}: {sugarValue} are set")
                matchedRule = r
                break

    if matchedRule is None or matchedRule['action'] not in __allRulesByName:
        raise KeyError(f"no rule matching {rule}")

    # OK, sanity check using the signature
    matchedRule = matchedRule.copy()
    ruleFunc = __allRulesByName[matchedRule['action']]
    action = matchedRule.pop('action')
    ruleSig = inspect.signature(ruleFunc)

    boundRule = ruleSig.bind(hdr, **matchedRule)
    boundRule.apply_defaults()

    if doRun:
        logger.info("running rule %s %s", action, matchedRule)
        ruleFunc(hdr, **boundRule.kwargs)
    else:
        logger.info("NOT running rule %s %s", action, matchedRule)

    return hdr

def _findPatchesForFitsName(fitsName, patchPath):
    """ Return the patches for a given FITS file.

    Args
    ----
    fitsName : str or pathlib.Path
      The name of a FITS file to find patches for.
    patchPath : str or pathlib.Path
      The name of an existing YAML file containing patches.

    Returns
    -------
    patches : a possibly empty tuple of patches.
    """

    logger.debug(f'looking for patches in {patchPath}')
    with open(patchPath, 'r') as patchFile:
        patchDict = yaml.safe_load(patchFile)

    fitsName = pathlib.PurePath(fitsName)
    testName = pathlib.PurePath(fitsName.name)
    patches = []
    for k,v in patchDict.items():
        logger.debug(f'checking (testName) vs. {k}')
        if testName.match(k):
            logger.debug(f'matched (testName) vs. {k}')
            patches.extend(v)

    return tuple(patches)

def _getPatchPaths(fitsName, header, yamlRootDir=None):
    """ Return the files which _could_ contain patches for the given FITS file.

    Args
    ----
    fitsName : path-like
      The name of the fits file, no path. No currently used.
    header : pyfits header
      We just used DATE-OBS from this.
    yamlRootDir : path-like
      We search in this directory for patch files

    Returns
    -------
    patchFiles : list of full pathlib.Paths
      All the existing files which _might_ contain patches. Should be searched in order.

    Notes
    -----
    The test scheme is to look in a directory of patches for all of:
      - a global file (all.yaml)
      - date-based files (2019-04-10, 2019-04, 2019)

    Obviously any selection scheme would work, but this one is easy to test with.
    """

    # Generate all possible file names
    yamlFileNames = ['all.yaml']
    try:
        dateObs = header['DATE-OBS']

        # We want the date only, no time.
        TisAt = dateObs.index('T')
        if TisAt > 0:
            dateObs = dateObs[:TisAt]
        Y, M, D = dateObs.split('-')
        yamlFileNames.append(f"{Y}.yaml")
        yamlFileNames.append(f"{Y}-{M}.yaml")

        yamlFileNames.append(f"{dateObs}.yaml")
    except Exception as e:
        logging.warning(f'failed to get DATE-OBS card: {e}')
        # patch for that better be in the global patch file....

    if yamlRootDir is None:
        yamlRootDir = './patches'  # Or some well-known directory

    # Filter down to _existing_ paths.
    yamlDir = pathlib.Path(yamlRootDir)
    yamlPaths = []
    for f in yamlFileNames:
        yamlPath = yamlDir / f
        if yamlPath.is_file():
            yamlPaths.append(yamlPath)

    return yamlPaths

def getPatchedHeader(fitsPath, hdu=0, yamlRootDir=None, doRun=True):
    """ Return the FITS header, possibly patched.

    Args
    ----
    fitsPath : path-like
      The full path of an existing fits file.
    hdu : int/str
      The HDU to return
    yamlRootDir : path-like
      We search in this directory for patch files
    doRun: bool
      Whether to actually do anything but print out what we would do.

    Returns
    -------
    hdr : pyfits.Header

    """

    fitsPath = pathlib.Path(fitsPath)
    header = pyfits.getheader(fitsPath, ext=hdu)
    fitsName = fitsPath.name

    # Get all possible patch file paths for the given FITS _name_.
    patchPaths = _getPatchPaths(fitsName, header, yamlRootDir=yamlRootDir)

    # Filter down to all actual patches which apply to this file
    allPatches = []
    for patchPath in patchPaths:
        patches = _findPatchesForFitsName(fitsName, patchPath)
        if patches:
            allPatches.extend(patches)

    # Apply the patches
    for p in allPatches:
        header = _checkOrRunRule(header, p, doRun=doRun)

    # Return (possibly) patched header.
    return header

def main(argv=None):
    if isinstance(argv, str):
        import shlex
        argv = shlex.split(argv)
    elif argv is None:
        import sys
        argv = sys.argv[1:]

    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('--yamlDir', type=str, default=None, help='directory to search for patch files')
    parser.add_argument('--dryrun', '-n', action='store_true', help='whether to actually do anything')
    parser.add_argument('--hdu', type=str, default='0', help='which HDU to process')
    parser.add_argument('--verbose', action='store_true', help='logging level')
    parser.add_argument('--printHeader', action='store_true', help='whether to print final header out')

    parser.add_argument('fitsPath', type=str, help='the FITS file to process')

    opts = parser.parse_args(args=argv)
    try:
        hdu = int(opts.hdu)
    except:
        hdu = opts.hdu

    logLevel = logging.DEBUG if opts.verbose else logging.INFO
    logging.basicConfig(level=logLevel)
    logger.setLevel(logLevel)

    hdr = getPatchedHeader(opts.fitsPath, hdu,
                           opts.yamlDir,
                           doRun=not opts.dryrun)

    if opts.printHeader:
        print(hdr.tostring(sep='\n', padding=False))
        return hdr
    else:
        print("SHUSH")
        return None

if __name__ == "__main__":
    main()
