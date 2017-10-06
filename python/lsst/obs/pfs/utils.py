"""
Utilities for working with PFS data
"""

def getLampElements(md):
    """Return a list of the elements found in the lamps that are on

    @param md: dafBase.PropertyList containing the header
    """

    md = md.toDict()
    elements = []
    for lamp in ['Ne', 'HgAr', 'Xe']:
        if md.get('W_AIT_SRC_%s' % (lamp,), False):
            if lamp == 'HgAr':
                elements.append('Ar')
                elements.append('Hg')
            elif lamp == 'HgCd':
                elements.append('Cd')
                elements.append('Hg')
                # the carrier gas is argon
                elements.append('Ar')
            else:
                elements.append(lamp)

    return elements
