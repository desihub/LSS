'''
Functions written by Claude and Arnaud de Mattia to help with LSS tasks.

Requires the cosmodesi environment

Tables: the one in-memory type the single-process stages of :mod:`mockfactory.desi` share.

:class:`astropy.table.Table`, the type desitarget and fiberassign already speak. The stages that
run over MPI, :mod:`~mockfactory.desi.base` and the building of target catalogs, keep
:class:`mpytools.Catalog`; :func:`as_table` is where the two meet, and it wraps a catalog's
columns without copying them.

Two habits keep a table as cheap as the structured arrays it replaced, which copied every
column each time one was added:

- a column is added with :func:`set_column`, which does not copy it, rather than by plain
  assignment, which does;
- columns are dropped or kept with ``remove_columns`` or ``keep_columns`` on
  ``table.copy(copy_data=False)``, since ``table[names]`` copies every column it keeps.

And one for speed: arithmetic on a :class:`~astropy.table.Column` is about 1.8 times slower
than on the array behind it, so a hot expression reads ``table[name].value``.

'''

import logging
import os

import numpy as np

from astropy.table import Table, vstack

logger = logging.getLogger('lsscat')

#: Value the survey pipeline uses for a column that a left join left empty.
NULL = 999999


def as_table(array):
    """
    Return ``array`` as an :class:`astropy.table.Table`, the table every stage takes and returns.

    A table stores its columns apart, so adding one costs that column and nothing else; the
    structured arrays this replaces copied every column each time one was added, and a random
    catalog gains about fifteen on its way through the stages.

    A :class:`~astropy.table.Table` is returned as is. An ``mpytools`` catalog, as the rest of
    ``mockfactory`` produces, is wrapped without copying its columns. A structured array, as
    :func:`fitsio.read` gives, has each column copied out, so that the table holds contiguous
    columns and does not keep the whole record array alive.
    """
    import mpytools
    if isinstance(array, Table):
        return array
    if isinstance(array, mpytools.Catalog):
        return Table({name: array.get(name, return_type=None) for name in array.columns()},
                     copy=False)
    array = np.asarray(array)
    return Table({name: np.ascontiguousarray(array[name]) for name in array.dtype.names},
                 copy=False)


def count_tiles_claude(array, tilelocids=False):
    """
    Return, per target, how many tiles could have reached it and which ones.

    ``NTILE`` is the number of distinct tiles and ``TILES`` stands for the set of them, so that
    a later stage can group targets by the overlap of tiles they sit under and measure how
    complete each such overlap is. The survey pipeline names the set by writing the tile
    identifiers out sorted and joined by ``-``; here the set is carried as a code instead, for
    the memory that name costs at the scale of a random catalog. See :func:`_group_code`.

    Parameters
    ----------
    array : array
        Potential assignments, with ``TARGETID``, ``TILEID`` and, for ``tilelocids``,
        ``TILELOCID``.
    tilelocids : bool, default=True
        Whether to also code the fiber locations, as ``TILELOCIDS``.
    """
    array = as_table(array)
    targetid, ntile, tiles = _group_code(array['TARGETID'], array['TILEID'])
    toret = Table({'TARGETID': targetid.astype(array['TARGETID'].dtype, copy=False),
                   'NTILE': ntile.astype('i8', copy=False), 'TILES': tiles}, copy=False)
    if tilelocids:
        set_column(toret, 'TILELOCIDS', _group_code(
            array['TARGETID'], array['TILELOCID'])[2])
    logger.info('counted tiles for {:d} targets, up to {:d} tiles each'
                .format(len(toret), int(ntile.max()) if len(ntile) else 0))
    return toret
