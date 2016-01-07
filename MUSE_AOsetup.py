# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
import pdb
from astropy import units
from astropy.vo.client import conesearch
from astropy.coordinates import ICRS, FK5
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
def MUSE_AOregion_DeepCall(ra,dec,color='magenta',verbose=True):
    """
    Generating region file for AO setup for Deep Field MUSE call
    http://www.eso.org/sci/activities/docs/Call_MUSE_DF.pdf

    --- INPUT ---

    --- EXAMPLE OF USE ---


    """


    regionstring = """
# Region file format: DS9 version 4.1
global color=green dashlist=8 3 width=1 font="helvetica 10 normal roman" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1
fk5
# composite(11:49:35.808,+22:23:55.85,0) || composite=1 color=%s width=3
box(21:29:25.876,-7:41:31.94,50",50",310) || # color=yellow width=2
box(21:29:27.461,-7:41:03.85,30",8",40) || # width=2
box(21:29:27.783,-7:41:55.72,30",8",310) || # width=2
box(21:29:23.969,-7:41:08.16,30",8",310) || # width=2
box(21:29:24.290,-7:42:00.03,30",8",40) || # width=2
box(21:29:25.876,-7:41:31.94,60",60",310) || # color=red width=2
circle(21:29:25.876,-7:41:31.94,52") || # color=blue width=4
circle(21:29:25.876,-7:41:31.94,101.5") || # color=blue width=4
circle(21:29:25.901,-7:41:31.57,2") # color=magenta width=4
"""  % (color)

    return regionstring
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
def perform_conesearch(ra,dec,search_radius,catalog=None,verbose=True):
    """
    Perform a conesearch using astropy's VO capabilities


    --- INPUT ---
    ra                  richt ascension [deg]
    dec                 declination     [deg]
    search_radius       radius around ra and dec to search [deg]
    catalo              name of catalog to search; set to None to get option printed
    --- EXAMPLE OF USE ---
    import MUSE_AOsetup as mao
    catalog  = 'guide*star'
    objtable = mao.perform_conesearch(177.3992,22.39885,0.5,catalog=catalog)
    """

    if catalog == None:
        print ' - No catalog selected, please choose from the following list (wildcards ok, e.g., "guide*star":'
        for cat in conesearch.list_catalogs(): print '  ',cat
        return None
    else:
        if verbose: print ' - Performing cone search in radius ',search_radius,\
            ' deg around (ra,dec) = (',ra,',',dec,')'

    cat_list = conesearch.list_catalogs(pattern=catalog)

    obj_coord = ICRS(ra=177.39920*units.deg, dec=22.39885*units.deg)

    search_result = conesearch.search_all(obj_coord, search_radius*units.deg, catalog_db=cat_list)

    if verbose:
        print ' - Searched ',len(search_result),' catalogs. Each catalog contains the following number of entries:'
        print ' '
        for url, tab in search_result.items():
            print '   CatUrl: ',url
            print '   TabCol: ',tab.nrows
            print '   Nrows : ',tab.array.data.dtype.names
            print ' '

    return search_result
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =