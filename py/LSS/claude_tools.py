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


def set_column(table, name, value, dtype=None):
    """
    Set column ``name`` of ``table`` to ``value``, typed as a structured array field would be.

    The column keeps its own type if it exists, and takes ``dtype`` if it does not; ``dtype``
    of None keeps the type of ``value``. A scalar fills every row. ``value`` is not copied,
    which plain assignment to a table would do.
    """
    if name in table.colnames:
        dtype = table[name].dtype
    value = np.asarray(value)
    if dtype is not None:
        value = value.astype(dtype, copy=False)
    if not value.ndim:
        value = np.full(len(table), value, dtype=value.dtype)
    if name in table.colnames:
        table.replace_column(name, value, copy=False)
    else:
        table.add_column(value, name=name, copy=False)
    return table


def encode_keys(*keys):
    """
    Turn several key columns into one integer label per distinct combination.

    The catalog stages join and group on keys that are either wide (``TARGETID``) or compound
    (``TARGETID``, ``LOCATION``, ``TILEID``), and packing them into a single integer by
    arithmetic overflows or collides. Each column is replaced by a dense code instead, and the
    codes are combined in a width that is known to fit.
    """
    keys = [np.asarray(key) for key in keys]
    code = np.zeros(len(keys[0]), dtype='i8')
    for key in keys:
        _, dense = np.unique(key, return_inverse=True)
        code = code * (dense.max() + 1 if len(dense) else 1) + dense
    return code


def last_of_each(key, sort=None, tie=None):
    """
    Return the index of one row per distinct value of ``key``: the one with the largest
    ``sort``, and among those the one with the smallest ``tie``.

    Ties are common and they matter. The rows a target can be kept at are ranked by what they
    say about it, and at a full survey about two in five targets have several rows at the top
    of that ranking, equally good by every criterion the ranking uses. Which one is kept still
    decides the fiber location the target is charged to, and so the completeness its
    neighbours are weighted by. The survey pipeline leaves the choice to an unstable sort,
    which makes its catalogs irreproducible at that level; here it is settled by ``tie``, so
    that the same inputs give the same catalog however the rows were assembled.

    The indices come back sorted, so the result is in the order of the input rather than in
    the order of the sort key.
    """
    size = len(key)
    if size == 0:
        return np.zeros(0, dtype='i8')
    position = np.arange(size)
    last = position if tie is None else -np.asarray(tie)
    order = np.lexsort(
        (last,) + ((position,) if sort is None else (sort,)) + (key,))
    sorted_key = np.asarray(key)[order]
    is_last = np.empty(size, dtype='?')
    is_last[-1] = True
    is_last[:-1] = sorted_key[1:] != sorted_key[:-1]
    return np.sort(order[is_last])


def group_fraction(key, weights):
    """
    Return, for each row, the mean of ``weights`` over the rows sharing its ``key``.

    The survey pipeline builds this as a dictionary in a Python loop over the distinct keys and
    then reads it back in a second loop over the rows; at the tens of millions of fiber
    locations of a full survey that is most of the cost of the stage.
    """
    _, dense, counts = np.unique(key, return_inverse=True, return_counts=True)
    return np.bincount(dense, weights=np.asarray(weights, dtype='f8'),
                       minlength=len(counts))[dense] / counts[dense]


def match(left, right):
    """
    Return, for each entry of ``left``, the index of the entry of ``right`` holding the same
    key, and -1 where there is none. ``right`` must hold each key at most once.
    """
    left, right = np.asarray(left), np.asarray(right)
    if not len(right):
        return np.full(len(left), -1, dtype='i8')
    order = np.argsort(right, kind='stable')
    # The queries are sorted before they are probed, and the answers scattered back. A binary
    # search over a sorted array of tens of millions is a cache miss at nearly every level, and
    # probing it in a random order pays that for every query; probing it in order walks the same
    # memory the array is laid out in. Measured on a join of 89 million rows against 34 million,
    # 93 s against 17 s, for the same indices.
    argsort = np.argsort(left, kind='stable')
    position = np.empty(len(left), dtype='i8')
    position[argsort] = np.searchsorted(right[order], left[argsort])
    index = order[np.clip(position, 0, len(order) - 1)]
    return np.where(right[index] == left, index, -1)


def join_left(left, right, keys, columns=None, fill=None, rename=None):
    """
    Add ``columns`` of ``right`` to ``left``, matched on ``keys``, filling the rows of ``left``
    that ``right`` has no entry for.

    This stands in for :func:`astropy.table.join` with ``join_type='left'``, for the case every
    join of the pipeline is in: the right-hand side holds each key once, so the result has
    exactly the rows of ``left``, in their order. It returns plain arrays rather than the masked
    columns astropy produces, so a filled value has to be given for each added column.

    Parameters
    ----------
    left, right : Table, array
        Tables, ``mpytools`` catalogs or structured arrays. ``right`` must hold each key
        combination at most once.
    keys : str, list
        Name, or names, of the columns to match on.
    columns : list, default=None
        Columns of ``right`` to add. Defaults to all of them but the keys.
    fill : dict, default=None
        Value to give a column where the key is absent from ``right``. Defaults to ``nan`` for
        a floating column, :data:`NULL` for an integer one, and to the type's zero otherwise,
        which is what the survey pipeline's masked columns come out as when written.
    rename : dict, default=None
        New name for a column of ``right``, for the cases where it would clash.

    Returns
    -------
    table : Table
        A new table; the columns of ``left`` are shared with it, not copied.
    """
    left, right = as_table(left), as_table(right)
    keys = [keys] if isinstance(keys, str) else list(keys)
    if columns is None:
        columns = [name for name in right.colnames if name not in keys]
    fill, rename = dict(fill or {}), dict(rename or {})
    if len(keys) == 1:
        index = match(left[keys[0]], right[keys[0]])
    else:
        # Encoding the two sides apart would give them unrelated labels; encode them together.
        size = len(left)
        code = encode_keys(
            *[np.concatenate([left[key], right[key]]) for key in keys])
        index = match(code[:size], code[size:])
    absent = index < 0
    logger.info('joined on {}: {:d} of {:d} rows unmatched'.format(
        keys, absent.sum(), len(left)))
    toret = left.copy(copy_data=False)
    for name in columns:
        column = right[name].value[np.where(absent, 0, index)]
        if absent.any():
            value = fill.get(name, None)
            if value is None:
                value = np.nan if column.dtype.kind == 'f' \
                    else NULL if column.dtype.kind in 'iu' else column.dtype.type()
            # Back to the type of ``right``: the fill value may have promoted it.
            column = np.where(absent.reshape((-1,) + (1,) * (column.ndim - 1)), value,
                              column).astype(right[name].dtype, copy=False)
        # A column ``left`` already has keeps its type, as a structured array's field would.
        set_column(toret, rename.get(name, name), column)
    return toret


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


def _group_code(targetid, value):
    """
    Return the distinct targets, how many distinct ``value`` each has, and a code standing for
    the sorted set of those values.

    The code is a hash of the set rather than a dense label, because the data and its randoms
    are grouped in separate calls and still have to agree: a later stage joins one to the other
    on this column, so the same set of tiles has to come out the same wherever it is met. A
    dense label would depend on which sets happened to be present in the call.

    Writing the set out as text agrees across calls too, and is what the survey pipeline does,
    but a name like ``1230-4560-7890`` takes a hundred and twenty bytes a row against eight,
    and at the thirty million rows of a random catalog that one column is four gigabytes.

    Two independent codes are accumulated in the same pass, which costs an exclusive or and a
    multiply on top of the masking that dominates the loop. Only the first is returned; the
    second is there to be checked against, so that a collision is raised rather than quietly
    merging two different overlaps of tiles and mis-weighting every target under them.
    """
    order = np.lexsort((value, targetid))
    sorted_targetid, sorted_value = np.asarray(
        targetid)[order], np.asarray(value)[order]
    distinct = np.empty(len(order), dtype='?')
    distinct[0] = True
    distinct[1:] = ((sorted_targetid[1:] != sorted_targetid[:-1])
                    | (sorted_value[1:] != sorted_value[:-1]))
    sorted_targetid, sorted_value = sorted_targetid[distinct], sorted_value[distinct]

    unique, start, counts = np.unique(
        sorted_targetid, return_index=True, return_counts=True)
    # Position of each value within its target's group, the values being already sorted.
    index = np.repeat(np.arange(len(unique)), counts)
    rank = np.arange(len(sorted_value)) - start[index]
    # Fowler-Noll-Vo over the values in sorted order, so that the code stands for the set and
    # not for the order the rows happened to arrive in.
    # the initial state for a fast, deterministic, non-cryptographic hash—likely to generate reproducible IDs, partitions, or assignments from catalog values. The corresponding FNV-64 prime is 1099511628211.
    code = np.full(len(unique), np.uint64(14695981039346656037), dtype='u8')
    # this is a different starting state for a second, independent hash stream
    other = np.full(len(unique), np.uint64(14313749767032793493), dtype='u8')
    # these are the two FNV-like odd multipliers used for the two accumulators
    prime, other_prime = np.uint64(1099511628211), np.uint64(880355133930503)
    width = counts.max() if len(counts) else 0
    for i in range(width):
        select = rank == i
        at = index[select]
        seen = sorted_value[select].astype('u8')
        code[at] = (code[at] ^ seen) * prime
        other[at] = (other[at] ^ seen) * other_prime
    _check_code_collision(code, other)
    return unique, counts.astype('i8'), code.view('i8')


def _check_code_collision(code, other):
    """
    Raise if two distinct sets share a code, judged by a second, independent code.

    Distinct sets that agree on both codes would have to collide in a hundred and twenty eight
    bits at once, so counting the codes and counting the pairs is as good as comparing the sets
    themselves, and costs one sort of the targets rather than the sets written out.
    """
    if not len(code):
        return
    order = np.lexsort((other, code))
    first, second = code[order], other[order]
    changed = first[1:] != first[:-1]
    ncode = 1 + int(np.count_nonzero(changed))
    npair = 1 + int(np.count_nonzero(changed | (second[1:] != second[:-1])))
    if npair != ncode:
        raise ValueError('{:d} distinct sets of tiles share {:d} codes: the 64 bit code has '
                         'collided, which would merge unrelated overlaps of tiles'
                         .format(npair, ncode))


def _group_names(targetid, value):
    """
    Return the distinct targets, how many distinct ``value`` each has, and those values sorted
    and joined by ``-``.
    """
    # Sorted on the two columns and cut to distinct pairs. Deliberately not a unique over a
    # structured array: numpy falls back to a generic element comparison for those, which at
    # the hundred million rows of a random catalog costs minutes rather than seconds.
    order = np.lexsort((value, targetid))
    sorted_targetid, sorted_value = np.asarray(
        targetid)[order], np.asarray(value)[order]
    distinct = np.empty(len(order), dtype='?')
    distinct[0] = True
    distinct[1:] = ((sorted_targetid[1:] != sorted_targetid[:-1])
                    | (sorted_value[1:] != sorted_value[:-1]))
    sorted_targetid, sorted_value = sorted_targetid[distinct], sorted_value[distinct]

    unique, start, counts = np.unique(
        sorted_targetid, return_index=True, return_counts=True)
    # Position of each value within its target's group, the values being already sorted.
    index = np.repeat(np.arange(len(unique)), counts)
    rank = np.arange(len(sorted_value)) - start[index]
    strings = np.char.mod('%d', sorted_value)
    width = counts.max() if len(counts) else 0
    toret = np.zeros(len(unique), dtype='U{:d}'.format(
        max(1, width * (strings.dtype.itemsize // 4 + 1))))
    for i in range(width):
        select = rank == i
        at = index[select]
        toret[at] = strings[select] if i == 0 else np.char.add(
            np.char.add(toret[at], '-'), strings[select])
    return unique, counts.astype('i8'), toret
