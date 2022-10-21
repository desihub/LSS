'''
Taken from desiutil/bitmask.py
'''

import yaml

class _MaskBit(int):
    """A single mask bit.
    Subclasses :class:`int` to act like an :class:`int`, but allows the
    ability to extend with blat.name, blat.comment, blat.mask, blat.bitnum.
    Attributes
    ----------
    name : :class:`str`
        The name of the bit.
    bitnum : :class:`int`
        The number of the bit.  The value of the bit is ``2**bitnum``.
    mask : :class:`int`
        The value of the bit, ``2**bitnum``.
    comment : :class:`str`
        A comment explaining the meaning of the bit.
    """
    def __new__(cls, name, bitnum, comment, extra=dict()):
        self = super(_MaskBit, cls).__new__(cls, 2**bitnum)
        self.name = name
        self.bitnum = bitnum
        self.mask = 2**bitnum
        self.comment = comment
        self._extra = extra
        for key, value in extra.items():
            if hasattr(self, key):
                raise AttributeError(
                    "Bit {0} extra key '{1}' is already in use by int objects.".format(name, key))
            self.__dict__[key] = value
        return self

    def __str__(self):
        return ('{0.name:16s} bit {0.bitnum} mask 0x{0.mask:X} - ' +
                '{0.comment}').format(self)

    # def __repr__(self):
    #     return "_MaskBit(name='{0.name}', bitnum={0.bitnum:d}, comment='{0.comment}')".format(self)


#  Class to provide mask bit utility functions
class BitMask(object):
    """BitMask object to represent bit names, masks, and comments.
    Typical users are not expected to create BitMask objects directly;
    other packages like desispec and desitarget will have used this
    to pre-create the bitmasks for them using definition files in those
    packages.
    Parameters
    ----------
    name : :class:`str`
        Name of this mask, must be key in `bitdefs`.
    bitdefs : :class:`dict`
        Dictionary of different mask bit definitions;
        each value is a list of ``[bitname, bitnum, comment]``.
        A 4th entry is optional, which must be a dictionary.
    """

    def __init__(self, name, bitdefs):
        """Init.
        """
        self._bits = dict()
        self._name = name
        for x in bitdefs[name]:
            bitname, bitnum, comment = x[0:3]
            if len(x) == 4:
                extra = x[3]
                if not isinstance(extra, dict):
                    raise ValueError(
                        '{} extra values should be a dict'.format(bitname))
            else:
                extra = dict()
            self._bits[bitname] = _MaskBit(bitname, bitnum, comment, extra)
            self._bits[bitnum] = self._bits[bitname]

    def __getitem__(self, bitname):
        """Return mask for individual bitname.
        """
        return self._bits[bitname]

    def bitnum(self, bitname):
        """Return bit number (int) for this `bitname` (string).
        Parameters
        ----------
        bitname : :class:`str`
            The bit name.
        Returns
        -------
        :class:`int`
            The bit value.
        """
        return self._bits[bitname].bitnum

    def bitname(self, bitnum):
        """Return bit name (string) for this `bitnum` (integer).
        Parameters
        ----------
        bitnum : :class:`int`
            The number of the bit.
        Returns
        -------
        :class:`str`
            The name of the bit.
        """
        return self._bits[bitnum].name

    def comment(self, bitname_or_num):
        """Return comment for this bit name or bit number.
        Parameters
        ----------
        bitname_or_num : :class:`int` or :class:`str`
            Name of number of the mask.
        Returns
        -------
        :class:`str`
            The comment string.
        """
        return self._bits[bitname_or_num].comment

    def mask(self, name_or_num):
        """Return mask value.
        Parameters
        ----------
        name_or_num : :class:`int` or :class:`str`
            Name of number of the mask.
        Returns
        -------
        :class:`int`
            The value of the mask.
        Examples
        --------
        >>> bitmask.mask(3)         #  2**3
        8
        >>> bitmask.mask('BLAT')
        >>> bitmask.mask('BLAT|FOO')
        """
        if isinstance(name_or_num, int):
            return self._bits[name_or_num].mask
        else:
            mask = 0
            for name in name_or_num.split('|'):
                mask |= self._bits[name].mask
            return mask

    def names(self, mask=None):
        """Return list of names of masked bits.
        Parameters
        ----------
        mask : :class:`int`, optional
            The mask integer to convert to names. If not supplied,
            return names of all known bits.
        Returns
        -------
        :class:`list`
            The list of names contained in the mask.
        """
        names = list()
        if mask is None:
            # return names in sorted order of bitnum
            bitnums = [x for x in self._bits.keys() if isinstance(x, int)]
            for bitnum in sorted(bitnums):
                names.append(self._bits[bitnum].name)
        else:
            mask = int(mask)  # workaround numpy issue #2955 for uint64
            bitnum = 0
            while 2**bitnum <= mask:
                if (2**bitnum & mask):
                    if bitnum in self._bits.keys():
                        names.append(self._bits[bitnum].name)
                    else:
                        names.append('UNKNOWN' + str(bitnum))
                bitnum += 1

        return names

    def __getattr__(self, name):
        """Enable ``mask.BITNAME`` equivalent to ``mask['BITNAME']``.
        """
        if name in self._bits:
            return self._bits[name]
        else:
            raise AttributeError('Unknown mask bit name ' + name)

    def __repr__(self):
        '''Return yaml representation defining the bits of this bitmask.
        '''
        result = list()
        result.append(self._name + ':')
        # return names in sorted order of bitnum
        bitnums = [x for x in self._bits.keys() if isinstance(x, int)]
        for bitnum in sorted(bitnums):
            bit = self._bits[bitnum]
            # format the line for single bit, with or without extra keys
            line = '  - [{:16s} {:2d}, "{}"'.format(
                bit.name+',', bit.bitnum, bit.comment)
            if len(bit._extra) > 0:
                line = line + ', '+str(bit._extra)+']'
            else:
                line = line + ']'
            result.append(line)

        return "\n".join(result)
    
    
# TODO: Wrap in main?
_bitdefs = yaml.safe_load('''
    lumfn_mask:
     - [DDP1ZLIM,     0, "Galaxy not in DDP limits"]
     - [FILLFACTOR,   1, "Fillfactor < 0.8"]
     - [INBGSBRIGHT,  2, "Galaxy not in BGS Bright"]
     - [CONSERVATIVE, 3, "Galaxy not conserved; see CONSERVATIVE mask"]
     - [DESI_HICOMP,  4, "High completeness region of DESI, e.g. 0.5 - 1.5 deg. of rosette."]
''')

_cbitdefs = yaml.safe_load('''
    consv_mask:
     - [DDP1ZLIM,     0, "Galaxy not in conservative DDP limits"]   
     - [BOUNDDIST,    1, "Boundary distance < 8"]
''')

lumfn_mask = BitMask('lumfn_mask', _bitdefs)
consv_mask = BitMask('consv_mask', _cbitdefs)
