#!/usr/bin/python
"""Reading and querying the TRANSP namelist (TR.DAT) as a dictionary of keys to values.

Functions over an explicit ``nmlist`` dict (and a logger), so the converter keeps the dictionary and a
test can pass its own.  ``nc2imas`` exposes these under their historical method names.
"""

import re
from typing import List

import numpy as np

from input_provenance import strip_fortran_comment
from name_patterns import make_namelist_re_pattern
from utilities import check_file_exists


# A Fortran real literal: ``1.0D-3``, ``2.5d0``, ``1D3``, ``.5``, ``-4.``
FORTRAN_REAL = re.compile(r'[+-]?(?:\d+\.?\d*|\.\d+)(?:[dDeE][+-]?\d+)?\Z')
_D_EXPONENT = re.compile(r'(?<=[0-9.])[dD](?=[+-]?\d)')


def normalize_exponents(text):
    """``1.0D-3`` -> ``1.0E-3``: Fortran D exponents written as the E that Python reads (numeric tokens only)."""
    return _D_EXPONENT.sub('E', text)


def fortran_number(text):
    """The value of a single Fortran numeric literal as int or float, or ``None`` if ``text`` is not one."""
    token = text.strip()
    try:
        return int(token)
    except ValueError:
        pass
    if FORTRAN_REAL.match(token):
        return float(normalize_exponents(token))
    return None


def parse_bool_value(value):
    """Parse TRANSP namelist logical values into Python bool."""
    if isinstance(value, (bool, np.bool_)):
        return bool(value)
    if isinstance(value, (int, np.integer)):
        return bool(value)
    if value is None:
        return None
    text = str(value).strip().lower().replace('.', '')
    if text in ('true', 't', 'yes', 'y', '1'):
        return True
    if text in ('false', 'f', 'no', 'n', '0'):
        return False
    return None


def get_bool(nmlist, name:str, default=None):
    """Return a scalar TRANSP logical from the loaded TR.DAT namelist.

    ``load_namelist`` preserves the spelling/case of keys.  Field-orientation
    controls are important enough that lookup here is deliberately
    case-insensitive.
    """
    if not nmlist:
        return default
    key = next((k for k in nmlist if k.strip().upper() == name.upper()), None)
    if key is None:
        return default
    value = parse_bool_value(nmlist[key])
    return default if value is None else value


def load_namelist(path, logger):
    """The namelist file as ``{key: value}`` (int or float where the text casts), or ``None`` without a path."""
    if path is None:
        print('**Error: namelist file not set')
        return None
    nmlist = {}
    if not check_file_exists(path):
        print(f'**WARNING: namelist file {path} not found')
        return nmlist

    with open(path, encoding = 'ISO-8859-1') as nmf:
        logger.info(f'loading {path}')
        for line in nmf:
            lns = line.strip()
            if lns and (lns[0] not in ['!', '&', '/']):
                lns = strip_fortran_comment(lns)           # '!' inside a quoted string is not a comment
                kv = lns.split('=',1)
                if len(kv) == 2:
                    val = kv[1].strip()
                    number = fortran_number(val)           # int or float where the text is one numeric literal
                    nmlist[kv[0].strip()] = val if number is None else number
    return nmlist


def parse_number(ss:str, tp):
    if not isinstance(ss, str):
        return tp(ss)

    #replace exponential 10** if any
    ss0 = ss
    ss = normalize_exponents(ss)
    res= ss.find('10**')
    if res >=0:
        ss = ss.replace('10**', 'E')
    opls = ['*', '/', '+', '-']
    nn = []
    ops = []
    p = 0

    #splitting into separate number strings
    for i,a in enumerate(ss):
        if a in opls and i>0 and ss[i-1].upper() !='E':
            nn.append(ss[p:i])
            ops.append(a)
            p = i+1

    if p<len(ss):
        nn.append(ss[p:])
    if len(ops)!=len(nn)-1:
        raise Exception(f'Error in parsing string to number {ss0}')
    n = 0
    try:
        n = tp(nn[0]) if len(nn)>0 else None
        for i,s in enumerate(ops):
            d = tp(nn[i+1])
            if s=='*':
                n = n * d
            elif s=='/':
                n = n/d
            elif s=='+':
                n = n+d
            elif s=='-':
                n = n-d
    except (ValueError, TypeError, IndexError):
        # an unparseable value must not become a plausible number: the caller gets an explicit unknown
        raise ValueError(f'cannot parse {ss0!r} as {tp.__name__} (parts {nn})') from None
    return n


def count_vars(nmlist, logger, name:str, ndig:int=2):
    if name[-1] != '?':
        raise Exception(f'Error: {name} needs "?" as last char')
    name = name[:-1]
    pattern = rf'^{re.escape(name)}(?:\((\d{{1,{ndig}}})\))?$'
    if nmlist is None:
        logger.info(f'**Error: {name}: namelist not found')
        return 0
    try:
        matches = [nmlist[k] for k in nmlist
                   if re.fullmatch(pattern, k)]
        return len(matches)
    except (AttributeError, TypeError):
        logger.info(f'**Error: {name}: namelist not found')
        return 0


def keys_re(nmlist, name:str, cut:str='', nd1:int=0, nd2:int=0, nch1:int=0,
      nch2:int=0, nd:int=0, nch:int=0) -> List:
    # name: regex need to be escaped \(
    pat = make_namelist_re_pattern(name, nd1, nd2, nch1, nch2, nd, nch)

    if pat is None:
        return []
    matches = [k for k in nmlist.keys() if re.search(pat, k)]
    if cut:
        matches = [s.split(cut,1)[1] for s in matches if len(s.split(cut,1))>1]
    return matches


#returns scalar or list
def get_var_(nmlist, logger, name:str, tp=None, scale=1., ind:int=-1):
    try:
        val = nmlist[name]
        if isinstance(val, str):
           val = re.split(r'[\s,;]+', val)
        if ind>=0:
           val = val[ind]
    except (KeyError, TypeError, IndexError, AttributeError):
        #logger.info(f'** {name} NOT found in name list')
        if nmlist is None:
            logger.info(f'** {nmlist} is not set')
        return None

    narr = 0
    try:
        if ind==-1:
           narr = len(val)
    except TypeError:
        pass

    if tp in [int, float]:
        try:
            if narr>0:
               val =  [parse_number(val[i], tp) for i in range(narr)]
               val = np.array(val)
            else:
               val = parse_number(val, tp)
        except ValueError as exc:
            # present but not a number: unknown, said aloud, never a plausible 0
            logger.warning(f'namelist {name}: {exc}; treated as unknown')
            return None
        if scale != 1.0:
           val = scale*val
    return val


# Some scalar can stand for arrays, e.g: FREQECH, which causes error when
#  an array is expected. So, call for array should handle single variable also.
# Returns scalar for single, list for array on RHS
def get_vars(nmlist, logger, name:str, tp=None, scale=1.0, ind:int=-1,
      size:int=1, resize:bool=False):
    #name: don't escape (); only regex needs \(
    # search for single or without () for array on rhs
    res = get_var_(nmlist, logger, name, tp, scale, ind=ind)
    n = 0
    # search for index
    if ind>=0 and res is None:
       namei = f'{name}({ind+1})'
       res = get_var_(nmlist, logger, namei, tp, scale)
       n = 1
    # search for array
    elif res is None:
       n = count_vars(nmlist, logger, f'{name}?')
       res = []
       for i in range(n):
          namei = f'{name}({i+1})'
          val = get_var_(nmlist, logger, namei, tp, scale)
          if val is not None:
             res.append(val)
       res = np.array(res)

    # resize only if scalar, not an array to different size
    if res is not None and resize:
       try:
           n = len(res)
       except TypeError:
           pass
       if n<=1:
           res = np.full(size, res)

    if size>1 and res is not None and len(np.array(res))!= size:
        raise Exception(f'resize failed: size of {name} != {size}')
    return res


def extract_str_between(logger, name:str, pat1:str, pat2:str, default:str, tp=None):
    res = None
    try:
        li = name.find(pat1) + len(pat1) + 1
        #can have multiple occurences. eg: ','
        ris = [s.start() for s in re.finditer(pat2, name)]
        ri = [i for i in ris if i>li][0]
        res = name[li:ri] if ri>li+1 else default
        if tp is not None:
            res = tp(res)
    except (IndexError, ValueError, TypeError):
        logger.debug(f'extract error: {name} doesnot have {pat1} and {pat2}')
    return res
